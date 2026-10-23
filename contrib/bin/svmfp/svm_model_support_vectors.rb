# Generate a .gfp file containing the support vectors of a svm-lite
# model.

# frozen_string_literal: true

require_relative 'lib/iwcmdline'

def usage(retcode)
  $stderr << "Generate a .gfp file with support vector weights\n"
  $stderr << "Takes two arguments: svm_lite model file, train.gfp\n"
  $stderr << " -o <fname>    output .gfp file name\n"
  $stderr << " -v            verbose output\n"
  exit(retcode)
end

# Extract the relationship between id and svm_lite weight from `fname`
def support_vectors_and_weights(input)
  threshold_b_rx = Regexp.new('(\S+) # threshold b, each following line')

  # First read through the file till we get to threshold_b.
  threshold_b = nil
  loop do
    line = input.readline
    break unless line

    m = threshold_b_rx.match(line)
    next unless m

    threshold_b = m[1]
    break
  end

  unless threshold_b
    $stderr << "Did not find threshold_b\n"
    return nil
  end

  # The rest of the records are support vectors.
  weights = {}
  sv_rx = Regexp.new('^(\\S+) .*# (\S+)')
  input.each do |line|
    m = sv_rx.match(line.chomp)
    raise "Invalid SV record #{line}" unless m

    weights[m[2]] = m[1]
  end

  $stderr << "threshold_b #{threshold_b} #{weights.size} sv->weight values\n"

  weights
end

def support_vectors_and_weights_file(fname)
  input = File.open(fname, 'r')
  rc = support_vectors_and_weights(input)
  input.close
  rc
end

# Extract PCN<...> from `tdt`
def id_from_tdt(tdt)
  tdt.split("\n").each do |line|
    next unless line.start_with?('PCN<')

    return line[4...-1]
  end
  nil
end

def gfp_subset(id_to_weight, input, output)
  tdts_read = 0
  tdts_written = 0
  while (tdt = input.gets("|\n"))
    tdts_read += 1
    id = id_from_tdt(tdt)
    unless id
      $stderr << "No identifier in #{tdt}\n"
      exit(1)
    end
    id = id.split[0]
    next unless id_to_weight.key?(id)

    output << tdt[0...-3] + "\nWEIGHT<#{id_to_weight[id]}>\n|\n"
    tdts_written += 1
  end

  $stderr << "Read #{tdts_read} tdts, wrote #{tdts_written}\n"
  if tdts_written != id_to_weight.size
    $stderr << "Have #{id_to_weight.size} support vectors, but wrote #{tdts_written} fingerprints\n"
    exit(1)
  end
  true
end

# Handle file opening and closing for `gfp_subset`.
def gfp_subset_file(id_to_weight, input_fname, output_fname)
  input = File.open(input_fname, 'r')
  output = File.open(output_fname, 'w')
  rc = gfp_subset(id_to_weight, input, output)
  input.close
  output.close
  rc
end

cmdline = IWCmdline.new('-v-o=s')

if cmdline.unrecognised_options_encountered
  $stderr << "unrecognised_options_encountered\n"
  usage(1)
end

unless cmdline.option_present('o')
  $stderr << "Must specify name of output .gfp file (-o)\n"
  usage(1)
end

if ARGV.size != 2
  $stderr << "Must specity two arguments, model file and gfp file\n"
  usage(1)
end

output_gfp_file = cmdline.value('o')

support_vectors = support_vectors_and_weights_file(ARGV[0])
gfp_subset_file(support_vectors, ARGV[1], output_gfp_file)
