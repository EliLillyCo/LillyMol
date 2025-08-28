#!/usr/bin/env ruby

# Read a QSAR dataset consisting of chiral variants and construct statistics
# on the degree to which chiral variants are different.

require_relative "lib/iwcmdline"

def usage(rc)
  $stderr << "Within a QSAR dataset identifies molecules that are chiral variants and creates\n"
  $stderr << "statistics about the degree of difference between chiral variants\n"
  $stderr << "Writes a histogram of chiral differences and counts to stdout.\n"
  $stderr << "Take that to your favourite plotting package\n"
  $stderr << " -E <fname>       acivity file\n"
  $stderr << " -xlabel <label>  header for the output file - likely the name of the target\n"
  $stderr << " -v               verbose output\n"
  exit(rc)
end

class StructureGroup attr_accessor :usmi, :smiles, :activity
  def initialize(usmi, smiles, activity)
    @usmi = usmi
    @smiles = []
    @smiles << smiles
    @activity = []
    @activity << activity
  end

  def extra(smiles, activity)
    @smiles << smiles
    @activity << activity
  end

  def number_smiles
    return @smiles.size
  end

  def min_activity
    return @activity.min
  end

  def max_activity
    return @activity.max
  end

  def chirality_present
    return @smiles.any? { |smiles| smiles.match(/@/) }
  end

  def all_smiles_same
    return @smiles.all? { |smi| smi == @smiles[0] }
  end

  def update_delta(deltas, expt_range)
    m1,m2 = @activity.minmax

    ndx = ((m2 - m1) / expt_range * 100.0).to_int

    # $stderr << "range #{m1} to #{m2} range #{expt_range} ndx #{ndx}\n"

    deltas[ndx] += 1
  end
end

# Given a sorted array of values, 'x', and a count associated with
# each value, return the approximate `percentile` value.
# We are working on the assumption that the number of values is
# significantly larger than 100. This is pretty rough but not
# worth improving.
def percentile(percentile, x, count)
  n = count.sum()
  npct = (n * percentile) / 100

  $stderr << "At #{percentile} and n #{n} npct #{npct}\n"

  return x[0] if npct == 0  # Hopefully not

  sum = 0
  count.each_with_index do |c, i|
    sum += c
    return x[i] if sum >= npct
  end

  # This should not happen
  return x.end
end

# Write Julia commands for plotting this data.
#   xs: x coordinates
#   count: y coordinates - the number of items at each X value.
#   label: name of response
#   activity_range: range of differences - used for normalisation.
#   ave: average activity - in the non normalised range.
#   fname: Julia file that will be created.
def write_julia_plot(xs, count, xlabel, activity_range, ave, fname)
  File.open(fname, "w") do |output|
    output << "using Plots\n"
    output << "x = ["
    xs.each_with_index do |x, ndx|
      output << ',' if ndx > 0
      output << x / activity_range
    end
    output << "]\n"
    output << "y = ["
    (0...xs.size).each do |i|
      output << ',' if i > 0
      output << count[i]
    end
    output << "]\n"
    output << "n = sum(y)\n"
    pct90 = "%.3f" % (percentile(90, xs, count) / activity_range)
    ave_norm = "%.3f" % (ave / activity_range)
    label = "ave #{ave_norm} 90% #{pct90}"
    output << "bar(x, y, ylabel=\"Count\", xlabel=\"#{xlabel} diff\", label=\"#{label}\","
    output << " legend=:right,"

    tmp = "%.3f" % activity_range
    title = "Stereoisomer Differences\n#{xlabel}\\n N $(n) range #{tmp}"
    output << " color=:blue, title=\"#{title}\")\n"

    output << "savefig(\"/tmp/#{xlabel}.png\")\n"
  end

  cmd = "julia #{fname}"
  system(cmd)
end

def main
  cl = IWCmdline.new("-v-E=sfile-xlabel=s-keep-Julia=s")

  if cl.unrecognised_options_encountered
    $stderr << "unrecognised_options_encountered\n"
    usage(1)
  end

  verbose = cl.option_present('v')

  unless cl.option_present('E')
    $stderr << "Must specify activity file name via the -E option\n";
    usage(1)
  end

  activity_fname = cl.value('E')

  if ARGV.empty?
    $stderr << "Insufficient arugments, must specify training set smiles as an argument\n"
    usage(1);
  end

  tmpstem = "/tmp/cstereo#{Process.pid}"
  tmpsmi = "#{tmpstem}x.smi"

  files_to_be_removed = []
  files_to_be_removed << tmpsmi

  train_smi = ARGV[0]

  cmd = "fetch_smiles_quick.sh -C 2 -c 1  #{activity_fname} #{train_smi} |"
  cmd << "sed -E -e 's/(\\S+) /\\1 \\1 /' |"
  cmd << "fileconv.sh -s 0 -o usmi -S - - |"
  cmd << "iwcut.sh -f 2,1,3,4 -|"
  cmd << "fileconv.sh -o usmi -S - - |"
  cmd << "iwcut.sh -f 2,1,3,4 - > #{tmpsmi}"
  $stderr << "Executing #{cmd}\n" if verbose
  system(cmd)

  unless File.size?(tmpsmi)
    $stderr << "#{cmd} failed\n"
    return 1
  end

  structure_group = {}
  lines_read = 0
  File.foreach(tmpsmi) do |line|
    lines_read += 1
    f = line.chomp.split
    usmi = f[0]
    if structure_group.key?(usmi)
      structure_group[usmi].extra(f[1], f[3].to_f)
    else
      structure_group[usmi] = StructureGroup.new(usmi, f[1], f[3].to_f)
    end
  end

  $stderr << "From #{lines_read} lines, read #{structure_group.size} structure groups\n" if verbose

  min_activity = 0
  max_activity = 0
  first_call = true
  structure_group.each do |k, v|
    if first_call
      min_activity = v.min_activity
      max_activity = v.max_activity
      first_call = false
    elsif v.min_activity < min_activity
      min_activity = v.min_activity
    elsif v.max_activity > max_activity
      max_activity = v.max_activity
    end
  end

  $stderr << "Activity btw #{min_activity} and #{max_activity}\n" if verbose

  structure_group.delete_if { |k, v|
    v.number_smiles == 1
  }
  
  $stderr << "After removing singletons #{structure_group.size} groups\n" if verbose

  structure_group.delete_if { |k, v|
    ! v.chirality_present
  }

  $stderr << "After removal of non chirality related dupes #{structure_group.size} groups\n" if verbose

  structure_group.delete_if { |k, v| v.all_smiles_same } 

  $stderr << "After removal of actual duplicates #{structure_group.size}\n"

  if structure_group.empty?
    $stderr << "No chirality related dupes\n"
    return 1
  end

  deltas = Array.new(100, 0)
  range = max_activity - min_activity
  structure_group.each do |k, v|
    v.update_delta(deltas, range)
  end

  last_nonzero = deltas.rindex { |v| v > 0}

  xlabel = if cl.option_present('xlabel')
             cl.value('xlabel')
           else
             'Diff'
           end
  $stdout << "#{xlabel} Count\n"
  dx = (max_activity - min_activity) / 100.0
  sum = 0.0
  n = 0

  xs = Array.new
  (0..last_nonzero).each do |i|
    x = (i * dx).round(4)

    sum += x * deltas[i]
    n += deltas[i]

    $stdout << "#{x} #{deltas[i]}\n"

    xs << x
  end

  ave = (sum / n.to_f).round(4);
  $stderr << "#{xlabel} weighted average diff #{ave}\n"

  if cl.option_present('Julia')
    fname = cl.value('Julia')
    write_julia_plot(xs, deltas, xlabel, max_activity - min_activity, ave, fname)
  end

  unless cl.option_present('keep')
    files_to_be_removed.each do |fname|
      File.unlink(fname)
    end
  end
end

main
