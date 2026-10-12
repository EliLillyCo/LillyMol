#!/usr/bin/env ruby

require 'csv'
require 'English'
require 'fileutils'
require 'json'
require 'open3'
require 'optparse'
require 'rbconfig'
require 'shellwords'
require 'tmpdir'

# Group activity changes by transformation, keeping measured and model-predicted
# changes separate. Run mode generates and looks up transformed molecules first;
# analyse mode consumes the same intermediate files without invoking those tools.
NAME_SEPARATOR = '__MCTA__'.freeze

Stats = Struct.new(:n, :min, :max, :mean, :median, keyword_init: true)

# Retain individual changes so the median can be calculated alongside the mean.
# Empty groups have a count of zero and no numerical summary values.
class Values
  def initialize
    @values = []
  end

  def add(value)
    @values << value
  end

  def stats
    return Stats.new(n: 0, min: nil, max: nil, mean: nil, median: nil) if @values.empty?

    sorted = @values.sort
    nvalues = sorted.length
    median = if nvalues.odd?
               sorted[nvalues / 2]
             else
               0.5 * (sorted[nvalues / 2 - 1] + sorted[nvalues / 2])
             end

    Stats.new(
      n: nvalues,
      min: sorted.first,
      max: sorted.last,
      mean: @values.sum / nvalues,
      median:
    )
  end
end

def die(message)
  warn message
  exit 1
end

def quote_command(argv)
  argv.map(&:to_s).shelljoin
end

def run_command(argv, verbose:)
  warn quote_command(argv) if verbose
  system(*argv)
  die "Command failed: #{quote_command(argv)}" unless $CHILD_STATUS.success?
end

def run_pipeline(commands, verbose:)
  warn commands.map { |cmd| quote_command(cmd) }.join(' | ') if verbose

  # Check every process: a successful database lookup must not hide a failure
  # in the upstream molecule generator.
  wait_threads = Open3.pipeline_start(*commands)
  wait_threads.each_with_index do |thread, ndx|
    status = thread.value
    next if status.success?

    die "Pipeline command #{ndx + 1} failed: #{quote_command(commands[ndx])}"
  end
end

def each_non_blank_line(fname)
  return enum_for(__method__, fname) unless block_given?

  File.foreach(fname).with_index(1) do |line, line_number|
    line.chomp!
    next if line.empty?

    yield line, line_number
  end
end

def read_smiles_ids(fname)
  ids = {}
  each_non_blank_line(fname) do |line, line_number|
    tokens = line.split
    die "#{fname}:#{line_number}: expected 2 tokens, got #{tokens.size}" unless tokens.size == 2

    _smiles, id = tokens
    die "#{fname}:#{line_number}: id contains reserved separator #{NAME_SEPARATOR}" if id.include?(NAME_SEPARATOR)
    die "#{fname}:#{line_number}: duplicate id '#{id}'" if ids.key?(id)

    ids[id] = true
  end

  die "#{fname}: no molecules read" if ids.empty?
  ids
end

def read_csv_activity(fname)
  activity = {}

  CSV.foreach(fname, headers: true).with_index(2) do |row, line_number|
    die "#{fname}:#{line_number}: expected 2 columns, got #{row.fields.size}" unless row.fields.size == 2

    id = row[0].to_s
    die "#{fname}:#{line_number}: duplicate id '#{id}'" if activity.key?(id)

    activity[id] = Float(row[1])
  rescue ArgumentError
    die "#{fname}:#{line_number}: invalid activity '#{row[1]}'"
  end

  activity
end

def read_space_activity(fname)
  activity = {}
  saw_header = false

  each_non_blank_line(fname) do |line, line_number|
    unless saw_header
      saw_header = true
      next
    end

    tokens = line.split
    die "#{fname}:#{line_number}: expected 2 columns, got #{tokens.size}" unless tokens.size == 2

    id, value = tokens
    die "#{fname}:#{line_number}: duplicate id '#{id}'" if activity.key?(id)

    activity[id] = Float(value)
  rescue ArgumentError
    die "#{fname}:#{line_number}: invalid activity '#{value}'"
  end

  activity
end

def read_activity(fname, ids = nil)
  activity = if File.extname(fname).casecmp?('.csv')
               read_csv_activity(fname)
             else
               read_space_activity(fname)
             end

  die "#{fname}: no activity values read" if activity.empty?

  # Run mode supplies the input molecule ids and requires activity for each.
  # Additional activity records are allowed, including database lookup matches.
  if ids
    ids.each_key do |id|
      die "#{fname}: missing activity for smiles id '#{id}'" unless activity.key?(id)
    end
  end

  activity
end


def extract_reaction_name(fname)
  contents = File.read(fname)

  # Reaction names can come from either textproto or legacy MSI reaction files.
  if (match = contents.match(/^\s*name:\s*"([^"]+)"/))
    return match[1]
  end

  if (match = contents.match(/\(A\s+C\s+Comment\s+"([^"]+)"\)/))
    return match[1]
  end

  nil
end

def default_reactions_file
  return nil unless ENV['LILLYMOL_HOME']

  fname = File.join(ENV.fetch('LILLYMOL_HOME'), 'data', 'MedchemWizard', 'REACTIONS')
  File.file?(fname) ? fname : nil
end

def read_reaction_file_map(fname)
  mapping = {}
  root = File.dirname(fname)

  each_non_blank_line(fname) do |line, line_number|
    next if line.start_with?('#')

    reaction_file = line.start_with?('PROTO:') ? line.sub(/\APROTO:/, '') : line
    path = File.join(root, reaction_file)
    die "#{fname}:#{line_number}: missing reaction file '#{reaction_file}'" unless File.file?(path)

    name = extract_reaction_name(path)
    die "#{fname}:#{line_number}: cannot determine reaction name in '#{reaction_file}'" unless name

    if name.match?(/\s/) || name.include?(NAME_SEPARATOR)
      die "#{fname}:#{line_number}: reaction name must contain no whitespace or #{NAME_SEPARATOR}: '#{name}'"
    end
    if mapping.key?(name)
      die "#{fname}:#{line_number}: duplicate transformation name '#{name}' in '#{reaction_file}' and '#{mapping[name]}'"
    end

    # Store the original REACTIONS token. For textproto reactions this preserves
    # the PROTO: prefix needed when writing a medchem_wizard reaction list.
    mapping[name] = line
  end

  mapping
end

# Database matches contain: smiles combined_name found_id. Legacy files with
# separate starting_id and transformation columns are also accepted.
# Compare the matched molecule's measured activity with the starting molecule's
# measured activity. Positive deltas mean the numerical activity increased;
# whether that is desirable depends on the activity scale supplied by the user.
def split_product_name(name, context)
  parts = name.split(NAME_SEPARATOR, -1)
  unless parts.size == 2 && parts.all? { |part| !part.empty? && !part.match?(/\s/) }
    die "#{context}: invalid product name '#{name}'; expected id#{NAME_SEPARATOR}reaction"
  end
  parts
end

# Check wizard output before lookup or scoring can truncate a malformed name.
def validate_generated_names
  STDIN.each_line.with_index(1) do |line, number|
    tokens = line.split
    next if tokens.empty?
    die "Wizard output:#{number}: expected SMILES and one product-name token" unless tokens.size == 2
    split_product_name(tokens[1], "Wizard output:#{number}")
    $stdout.write(line)
  end
end

def accumulate_found(fname, activity, observed)
  each_non_blank_line(fname) do |line, line_number|
    tokens = line.split
    if tokens.size == 3
      _smiles, name, found_id = tokens
      starting_id, transformation = split_product_name(name, "#{fname}:#{line_number}")
    elsif tokens.size == 4
      _smiles, starting_id, transformation, found_id = tokens
    else
      die "#{fname}:#{line_number}: expected 3 or 4 tokens, got #{tokens.size}"
    end
    die "#{fname}:#{line_number}: unknown starting id '#{starting_id}'" unless activity.key?(starting_id)
    die "#{fname}:#{line_number}: unknown found id '#{found_id}'" unless activity.key?(found_id)

    observed[transformation].add(activity[found_id] - activity[starting_id])
  end
end

# Optional predicted-value input. The first non-blank line is a header.
# Subsequent records contain:
#   starting_id__MCTA__transformation predicted_activity
# Legacy starting_id transformation predicted_activity records are also accepted.
# Predictions are absolute activities for transformed molecules. Subtract the
# starting molecule's measured activity to put them on the observed delta scale.
def accumulate_predictions(fname, activity, predicted)
  skip_header = true
  each_non_blank_line(fname) do |line, line_number|
    if skip_header
      skip_header = false
      next
    end

    tokens = line.split
    if tokens.size == 2
      name, value = tokens
      starting_id, transformation = split_product_name(name, "#{fname}:#{line_number}")
    elsif tokens.size == 3
      starting_id, transformation, value = tokens
    else
      die "#{fname}:#{line_number}: expected 2 or 3 tokens, got #{tokens.size}"
    end
    die "#{fname}:#{line_number}: unknown starting id '#{starting_id}'" unless activity.key?(starting_id)

    predicted[transformation].add(Float(value) - activity[starting_id])
  rescue ArgumentError
    die "#{fname}:#{line_number}: invalid predicted activity '#{value}'"
  end
end


def command_from_template(command, input, output)
  die 'Empty model command' if command.strip.empty?

  # Preserve shell syntax such as pipes and stdout redirection. Escape inserted
  # filenames, replacing any quotes around a standalone placeholder as well.
  # Example: score_model {input} > {output}
  used_input = false
  used_output = false
  command = command.gsub(/'\{input\}'|"\{input\}"|\{input\}/) do
    used_input = true
    Shellwords.escape(input)
  end
  command = command.gsub(/'\{output\}'|"\{output\}"|\{output\}/) do
    used_output = true
    Shellwords.escape(output)
  end

  # Templates without placeholders retain the positional input/output convention.
  command += " #{Shellwords.escape(input)}" unless used_input
  command += " #{Shellwords.escape(output)}" unless used_output
  ['sh', '-c', command]
end

def format_stat(value)
  return '.' if value.nil?

  format('%.6g', value)
end

def write_table(output, observed, predicted)
  output.puts %w[
    transformation
    observed_n observed_min observed_max observed_mean observed_median
    predicted_n predicted_min predicted_max predicted_mean predicted_median
  ].join(' ')

  # Include transformations seen in either source. The missing source gets an
  # empty group, displayed as count 0 and dots for unavailable statistics.
  transformations = (observed.keys + predicted.keys).uniq.sort
  transformations.each do |transformation|
    observed_stats = observed[transformation].stats
    predicted_stats = predicted[transformation].stats

    output.puts [
      transformation,
      observed_stats.n,
      format_stat(observed_stats.min),
      format_stat(observed_stats.max),
      format_stat(observed_stats.mean),
      format_stat(observed_stats.median),
      predicted_stats.n,
      format_stat(predicted_stats.min),
      format_stat(predicted_stats.max),
      format_stat(predicted_stats.mean),
      format_stat(predicted_stats.median)
    ].join(' ')
  end
end


def summary_json(stats)
  # Omit absent statistics rather than writing nulls in the profile JSON.
  result = { 'n' => stats.n }
  result['min'] = stats.min unless stats.min.nil?
  result['max'] = stats.max unless stats.max.nil?
  result['mean'] = stats.mean unless stats.mean.nil?
  result['median'] = stats.median unless stats.median.nil?
  result
end

def activity_effect(property, source, stats)
  {
    'property' => property,
    'source' => source,
    'delta' => summary_json(stats)
  }
end

def transformation_profiles(observed, predicted, reaction_files)
  # Resolve display names back to REACTIONS entries so each exported profile
  # identifies the reaction file as well as its observed and predicted effects.
  transformations = (observed.keys + predicted.keys).uniq.sort
  transformations.map do |transformation|
    reaction_file = reaction_files[transformation]
    die "No reaction file for transformation '#{transformation}'" unless reaction_file

    {
      'transformation' => transformation,
      'reactionFile' => reaction_file,
      'effect' => [
        activity_effect('activity', 'OBSERVED_ACTIVITY', observed[transformation].stats),
        activity_effect('activity', 'PREDICTED_ACTIVITY', predicted[transformation].stats)
      ]
    }
  end
end

def write_profile_json(output, observed, predicted, reaction_files)
  output.puts JSON.pretty_generate({ 'profile' => transformation_profiles(observed, predicted, reaction_files) })
end

def usage(parser, mode = nil)
  command = mode ? " #{mode}" : ' [run|generate|analyse]'
  warn <<~USAGE
    Accumulate activity changes by medchem_wizard transformation.

    Usage: #{File.basename($PROGRAM_NAME)}#{command} [options] input.smi

    Modes:
      run       Generate transformed molecules, optionally score not-found molecules, and analyse.
      generate  Generate transformed molecules and split them into found/notfound files.
      analyse   Analyse existing found and prediction files.

    If no mode is specified, run mode is used for backward compatibility.

    #{parser}

    The smiles input must contain exactly two tokens per non-blank line:
      smiles id

    The activity file must contain a header and two columns:
      id activity

    If the activity file name ends in .csv, it is parsed as CSV.
    Otherwise it is parsed as whitespace separated text.

    The --model-command option might be something like
      --model-command 'xgbd_evaluate.sh -mdir /path/to/model {input} > {output}'
  USAGE
  exit 1
end

def default_options
  {
    buildsmidb: 'buildsmidb_bdb',
    in_database: 'in_database_bdb',
    medchem_wizard: 'medchem_wizard.sh',
    output: '-',
    profile_json: nil,
    reactions: nil,
    verbose: false,
    keep: false
  }
end

def add_common_options(opts, options)
  opts.on('-A', '--activity FILE', 'Activity file') { |value| options[:activity] = value }
  opts.on('-d', '--database NAME', 'BerkeleyDB database name') { |value| options[:database] = value }
  opts.on('-o', '--output FILE', 'Output table, default stdout') { |value| options[:output] = value }
  opts.on('--profile-json FILE', 'Write TransformationProfileCollection JSON') do |value|
    options[:profile_json] = value
  end
  opts.on('--reactions FILE', 'REACTIONS file used to map transformation names to reaction files') do |value|
    options[:reactions] = value
  end
  opts.on('--found FILE', 'Found structures file') { |value| options[:found] = value }
  opts.on('--notfound FILE', 'Not-found structures file') { |value| options[:notfound] = value }
  opts.on('--predictions FILE', 'Predictions: combined_name value (or legacy id transformation value)') do |value|
    options[:predictions] = value
  end
  opts.on('--model-command COMMAND', 'Shell command producing predictions; use unquoted {input}/{output} placeholders') do |value|
    options[:model_command] = value
  end
  opts.on('--buildsmidb EXE', 'buildsmidb executable, default buildsmidb_bdb') do |value|
    options[:buildsmidb] = value
  end
  opts.on('--in-database EXE', 'in_database executable, default in_database_bdb') do |value|
    options[:in_database] = value
  end
  opts.on('--medchem-wizard EXE', 'medchem_wizard executable, default medchem_wizard.sh') do |value|
    options[:medchem_wizard] = value
  end
  opts.on('--keep', 'Keep temporary found/notfound files') { options[:keep] = true }
  opts.on('-v', '--verbose', 'Verbose execution') { options[:verbose] = true }
end

def parse_command_line(argv)
  mode = if %w[run generate analyse].include?(argv.first)
           argv.shift
         else
           'run'
         end

  options = default_options
  parser = OptionParser.new do |opts|
    opts.banner = "Usage: #{File.basename($PROGRAM_NAME)} #{mode} [options] input.smi"
    add_common_options(opts, options)
  end

  begin
    parser.parse!(argv)
  rescue OptionParser::ParseError => e
    die e.message
  end

  [mode, options, parser, argv]
end

def generate_transformed_molecules(smiles, options)
  # Index the starting molecules, then stream wizard products into the lookup.
  # Lookup writes matching products to found and unmatched products to notfound;
  # the latter are candidates for model scoring.
  run_command(
    [options.fetch(:buildsmidb), '-d', options.fetch(:database), '-g', 'all', '-l', '-c', smiles],
    verbose: options.fetch(:verbose)
  )

  medchem = [options.fetch(:medchem_wizard), '-W', NAME_SEPARATOR, '-V', '.', smiles]
  lookup = [
    options.fetch(:in_database), '-d', options.fetch(:database), '-l', '-c', '-p',
    '-F', options.fetch(:found), '-U', options.fetch(:notfound), '-i', 'smi', '-'
  ]
  validator = [RbConfig.ruby, File.expand_path(__FILE__), '--validate-generated-names']
  run_pipeline([medchem, validator, lookup], verbose: options.fetch(:verbose))
end

def analyse_transformations(options, activity_ids = nil)
  activity = read_activity(options.fetch(:activity), activity_ids)

  # Each transformation accumulates its own samples in each activity source.
  observed = Hash.new { |hash, key| hash[key] = Values.new }
  predicted = Hash.new { |hash, key| hash[key] = Values.new }

  accumulate_found(options.fetch(:found), activity, observed)
  accumulate_predictions(options.fetch(:predictions), activity, predicted) if options[:predictions]

  if options.fetch(:output) == '-'
    write_table($stdout, observed, predicted)
  else
    File.open(options.fetch(:output), 'w') { |file| write_table(file, observed, predicted) }
  end

  return unless options[:profile_json]

  reactions = options[:reactions] || default_reactions_file
  die 'Must specify --reactions or set LILLYMOL_HOME for --profile-json' unless reactions

  reaction_files = read_reaction_file_map(reactions)
  File.open(options.fetch(:profile_json), 'w') do |file|
    write_profile_json(file, observed, predicted, reaction_files)
  end
end

def run_mode(options, argv, parser)
  usage(parser, 'run') unless argv.size == 1
  usage(parser, 'run') unless options[:activity] && options[:database]
  die 'Specify only one of --predictions and --model-command' if options[:predictions] && options[:model_command]

  smiles = argv.fetch(0)
  ids = read_smiles_ids(smiles)

  # Allocate temporary paths only for intermediate files the caller did not name.
  # Explicitly named files survive cleanup; --keep also retains temporary files.
  tmpdir = nil
  unless options[:found] && options[:notfound] && (options[:predictions] || !options[:model_command])
    tmpdir = Dir.mktmpdir('medchem_transform_activity')
    options[:found] ||= File.join(tmpdir, 'found.smi')
    options[:notfound] ||= File.join(tmpdir, 'notfound.smi')
    options[:predictions] ||= File.join(tmpdir, 'predictions.txt') if options[:model_command]
  end

  begin
    generate_transformed_molecules(smiles, options)

    # Scorers return a header and combined product-name/value records. Legacy
    # three-column predictions remain accepted by the analysis parser.
    if options[:model_command]
      model = command_from_template(options.fetch(:model_command), options.fetch(:notfound),
                                    options.fetch(:predictions))
      run_command(model, verbose: options.fetch(:verbose))
    end

    analyse_transformations(options, ids)
  ensure
    FileUtils.remove_entry(tmpdir) if tmpdir && !options.fetch(:keep)
  end
end

def generate_mode(options, argv, parser)
  usage(parser, 'generate') unless argv.size == 1
  usage(parser, 'generate') unless options[:database]

  options[:found] ||= 'found.smi'
  options[:notfound] ||= 'notfound.smi'

  generate_transformed_molecules(argv.fetch(0), options)
end

def analyse_mode(options, argv, parser)
  usage(parser, 'analyse') unless argv.empty?
  usage(parser, 'analyse') unless options[:activity] && options[:found]
  analyse_transformations(options)
end

if ARGV == ['--validate-generated-names']
  validate_generated_names
  exit 0
end

mode, options, parser, remaining_args = parse_command_line(ARGV)

case mode
when 'run'
  run_mode(options, remaining_args, parser)
when 'generate'
  generate_mode(options, remaining_args, parser)
when 'analyse'
  analyse_mode(options, remaining_args, parser)
else
  die "Unrecognised mode '#{mode}'"
end
