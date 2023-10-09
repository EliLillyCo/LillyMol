# Collate individually generated replacement rings into a single set.
# Files that have come from ring_extraction.
#  ring_extraction.v2.sh -S rings -v -R 7 -k -x . -c -l $C/cas/casgO.smi > rings.log 2>&1 &

require 'find'

c3tk_bin = ENV['C3TK_BIN']
raise "No C3TK_BIN" unless File.directory?(c3tk_bin)
require "#{c3tk_bin}/ruby/lib/iwcmdline"

def usage(rc)
  $stderr << "Collates outputs from ring_extraction.\n"
  $stderr << "Arguments must be directories in which ring_extraction has been run\n"
  $stderr << "creating one or more files with names 'rings_*smi'\n"
  $stderr << " -destdir <dir>     directory where new files are created (def .)\n"
  $stderr << " -v                 verbose output\n"
  exit(rc)
end

cmdline = IWCmdline.new("-v-destdir=dir")
if cmdline.unrecognised_options_encountered
  $stderr << "unrecognised_options_encountered\n"
  usage(1)
end

verbose = cmdline.option_present('v')

destination = if cmdline.option_present('destdir')
                destination = cmdline.value('destdir')
              else
                destination = '.'
              end

# A class that is used to accumulate results for a given ring pattern.
# We just hold the textproto string in its entirety, and parse the n:
# attribute.
class Ring
  attr_reader :number_examples
  def initialize(line)
    @line = line
    m = / n: ([0-9]+) /.match(line)
    raise "Invalid data #{line}" unless m
    @number_examples = m[1].to_i
  end

  # Parse a textproto for the n: field and increment @number_examples
  def increment(line)
    m = / n: ([0-9]+) /.match(line)
    raise "Invalid data #{line}" unless m
    @number_examples += m[1].to_i
  end

  def to_string()
    return @line.gsub(/ :n [0-9]+ /, " n: #{@number_examples} ")
  end
end

if ARGV.empty?
  $stderr << "Insuficcient arguments\n"
  usage(1)
end

# For each kind of ring, 5a6A, a mapping from unique smiles
# to smiles/smarts/id
result = {}

name_rx = Regexp.new("rings.*_(\\S+)\\.smi")
ARGV.each do |d|
  Find.find(d) do |fname|
    m = name_rx.match(fname)
    next unless m
    rtype = m[1]
    $stderr << "#{fname} #{rtype}\n" if verbose
    if ! result.key?(rtype)
      result[rtype] = {}
    end
    File.open(fname, 'r').each do |line|
      line.chomp!
      f = line.split
      usmi = f[-1]
      if result[rtype].key?(usmi)
        result[rtype][usmi].increment(line)
      else
        result[rtype][usmi] = Ring.new(line)
      end
    end
  end
end

if verbose
  $stderr << "Found data on #{result.keys.length} ring types\n"
end

result.each do |rtype, examples|
  fname = File.join(destination, "#{rtype}.smi")
  $stderr << "Creating #{fname}\n" if verbose
  File.open(fname, 'w') do |output|
    examples.sort_by {|k, v| -v.number_examples}.each do |usmi, ringdata|
      output << ringdata.to_string() << "\n"
    end
  end
end
