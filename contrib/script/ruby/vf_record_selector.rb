# selects records for vf.
# rubocop:disable Layout/EmptyLineAfterGuardClause

require_relative 'lib/iwcmdline'

def usage
  $stderr << "Selects records from a file\n"
  $stderr << " -f <n>        skip the first <n> records of the file\n"
  exit 1
end

# Information about how records are filtered.
class Filters
  def initialize
    @min_smiles_length = 0
    @max_smiles_length = 0
    @rx_must_have = nil
    @rx_exclude = nil
    @tests = 0
    # The filters need to know how to find the smiles
    @input_separator = ' '
    @smiles_column = 1
  end

  def build(cl) # rubocop:disable Naming/UncommunicativeMethodParamName
    if cl.option_present('slen')
      @min_smiles_length = cl.value('slen')
      @tests += 1
    end
    if cl.option_present('SLEN')
      @max_smiles_length = cl.value('SLEN')
      @tests += 1
    end
    if cl.option_present('t')
      @rx_must_have = Regexp.new(cl.value('t'))
      @tests += 1
    end
    if cl.option_present('T')
      @rx_exclude = Regexp.new(cl.value('T'))
      @tests += 1
    end
    internally_consistent
  end

  def internally_consistent
    if @min_smiles_length > @max_smiles_length
      $stderr << "min #{@min_smiles_length} and max #{@max_smiles_length} inconsistent\n"
      return false
    end

    true
  end

  def active
    @tests > 0
  end

  def ok(line)
    return true if @tests.zero?

    f = line.split(@input_separator)
    return false if @min_smiles_length > 0 &&
                    f[@smiles_column].length < @min_smiles_length
    return false if @max_smiles_length > 0 &&
                    f[@smiles_column].length > @max_smiles_length
    return false if @rx_must_have &&
                    !@rx_must_have.match(line)
    true
  end
end

# Information about how records are selected.
class SelectionCriteria
  def initialize
    @skip_first = 0
    @select_all = false
    @nsel = 0
    @fraction_to_select = 0.0
    @z_sample = 0
    @select_every = 0
    @random_sample = 0
    @following_records = 0
    @columns_to_write = []
    # Only considered if columns_to_write is specified. Maybe fix sometime.
    @output_separator = ' '

    # During selection, we keep track of the number of records selected.
    @selected = 0
  end

  def build(cl) # rubocop:disable Naming/UncommunicativeMethodParamName
    @skip_first = cl.value('f') if cl.option_present('f')
    @select_all = cl.option_present('a')
    @nsel = cl.value('n') if cl.option_present('n')
    @fraction_to_select = cl.value('p') if cl.option_present('p')
    @z_sample = cl.value('z') if cl.option_present('z')
    @select_every = cl.value('e') if cl.option_present('e')
    @random_sample = cl.option_present('Z')
    @following_records = cl.value('j') if cl.option_present('j')
    @columns_to_write = cl.values('col').map { |col| col.to_i - 1}
    internally_consistent
  end

  def internally_consistent
    if @select_all
      if @fraction_to_select > 0.0 ||
         @nsel > 0 ||
         @random_sample ||
         @z_sample > 0 ||
         @select_every > 0
        $stderr << "The -a option is not compatible with -n -p -z -e -Z options\n"
        $stderr << "#{@nsel} #{@random_sample} #{@z_sample} #{@select_every}\n"
        return false
      end
      return true
    end
    if @random_sample
      if @fraction_to_select > 0.0 ||
         @z_sample > 0 ||
         @select_every > 0
        $stderr << "The -Z option is incompatible with the -e and -z options\n"
        return false
      end
    end

    true
  end

  # Return true if we have selected our required number of records
  def nsel_satisfied
    @selected >= @nsel
  end

  # If all tokens are being written, just write `line` to `output`.
  # Otherwise, write the columns in @columns_to_write
  def do_write(line, output)
    if @columns_to_write.empty?
      output << line
      return true
    end
    f = line.split(@input_separator)
    buffer = []
    @columns_to_write.each do |col|
      buffer << f[col]
    end
    output << buffer.join(@output_separator) << "\n"
  end

  def maybe_write(line, filters, output)
    return false unless filters.ok(line)
    do_write(line, output)
  end

  # Return whether or not our criteria are consistent with reading from stdin.
  def ok_stdin
    if @random_sample && @random_sample > 0
      $stderr << "Random sample not consistent with stdin\n"
      return false
    end
    true
  end

  def echo_following_records(input, output)
    @following_records.times do
      return true if input.eof
      do_write(input.readline, output)
    end
    true
  end

  def do_select_all(input, filters, output)
    input.each do |line|
      maybe_write(line, filters, output)
    end
    true
  end

  # Simple case of just writing @nsel items
  def do_select_n(input, filters, output)
    while @selected < @nsel && !input.eof
      next unless maybe_write(input.readline, filters, output)
      @selected += 1
    end
    true
  end

  def do_select_every(input, filters, output)
    line_number = @skip_first - 1
    while @selected < @nsel && !input.eof
      line_number += 1
      next unless (line_number % @select_every).zero?
      next unless maybe_write(line, filters, output)
      @selected += 1
      echo_following_records(input, output)
    end
    true
  end

  def random_sample_inner(input, file_size, filters, output)
    start_pos = input.tell
    processed = {}
    attempts = 0
    while attempts < 3 * @nsel && @selected < @nsel
      attempts += 1
      # 10 bytes is just a guess, covers previous line as well...
      offset = rand(start_pos..file_size - 10)
      input.seek(offset, IO::SEEK_SET)
      input.readline
      next if input.eof
      pos = input.tell
      next if processed.key?(pos)
      processed[pos] = true
      next unless maybe_write(input.readline, filters, output)
      @selected += 1
      echo_following_records(input, output)
    end
    true
  end

  def do_random_sample(input, filters, output)
    file_size = input.size
    random_sample_inner(input, file_size, filters, output)
  end

  # Helps implement z_sample. Return a set of line numbers.
  # Note that this is kind of wrong in the case of following
  # records. It is kind of ambiguous what it would mean.
  # Not sure anyone uses this functionality...
  def random_lines
    result = {}
    attempts = 0
    while attempts < 3 * @nsel && @selected < @nsel
      attempts += 1
      line_number = rand(@z_sample)
      next if result.key?(line_number)
      result[line_number] = true
      @selected += 1
      @following_records.times do |i|
        result[line_number + i + 1] = true
      end
    end
    @selected = 0
    result
  end

  # A random sample where the number of records int the file
  # is specified - the -z option.
  # This is broken in the case of filters, which are applied here
  # after the lines to be written have been chosen. Don't care...
  def do_subset_random_sample(input, filters, output)
    lines_to_write = random_lines
    line_number = -1

    while @selected < lines_to_write.size && !input.eof
      line = input.readline
      line_number += 1
      next unless lines_to_write.key?(line_number)
      next unless maybe_write(line, filters, output)
      @selected += 1
    end
    true
  end

  def do_select_via_probability(input, filters, output)
    while @selected < @nsel && !input.eof
      line = input.readline
      next unless rand(0.0..1.0) < @fraction_to_select
      next unless maybe_write(line, filters, output)
      @selected += 1
      echo_following_records(input, output)
    end
  end

  def do_skip_first(input)
    @skip_first.times do |i|
      input.readline
      if input.eof
        $stderr << "Skipping first eof #{i}\n"
        return false
      end
    end
    true
  end

  def select(input, filters, output)
    return false unless do_skip_first(input)
    # This is not correct, need more exclusion conditions
    @nsel = 30 unless (@nsel > 0 || @select_all)

    return do_select_all(input, filters, output) if @select_all
    return do_random_sample(input, filters, output) if @random_sample
    return do_subset_random_sample(input, filters, output) if @z_sample > 0
    return do_select_every(input, filters, output) if @select_every > 0
    return do_select_via_probability(input, filters, output) if @fraction_to_select > 0.0
    do_select_n(input, filters, output)
  end
end

def main
  cl = IWCmdline.new('-v-a-n=ipos-p=float-z=ipos-Z-m=ipos-j=ipos-col=list' \
                     '-e=ipos-slen=ipos-SLEN=ipos-f=ipos-grep=s-grepv=s')
  if cl.unrecognised_options_encountered
    $stderr << "unrecognised_options_encountered\n"
    usage
  end
  selection_criteria = SelectionCriteria.new
  usage unless selection_criteria.build(cl)

  filters = Filters.new
  usage unless filters.build(cl)

  if ARGV.empty?
    $stderr << "No input\n"
    usage
  end

  if ARGV[0] == '-'
    exit 1 unless selection_criteria.ok_stdin
    if selection_criteria.select($stdin, filters, $stdout)
      exit 0
    else
      exit 1
    end
  end

  ARGV.each do |fname|
    inp = File.open(fname, 'r')
    exit 1 unless selection_criteria.select(inp, filters, $stdout)
    inp.close
    break if selection_criteria.nsel_satisfied
  end
  exit 0
end

# rubocop:enable Layout/EmptyLineAfterGuardClause

main
