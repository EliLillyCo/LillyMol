#! /usr/bin/env ruby

ianhome = ENV['C3TK_BIN'] #"/home/rx87851"

require "#{ianhome}/ruby/lib/iwcmdline.rb"

$expert = false

def usage (rc)
  $stderr.print "Summarises results in the output files from svmfp_calibrate.rb\n"
  $stderr.print " -rx <regexp>   process files that match <rx>, -rx 'A80' for example\n"
  $stderr.print "                useful if too many files to process via cmd line\n"
  $stderr.print " -F <fname>     read files to be processed from <fname>\n"
  $stderr.print " -R B2          result of interest is Bsquared\n"
  $stderr.print " -R AE          result of interest is average absolute error\n"
  $stderr.print " -R RMS         result of interest is root mean square error\n"
  $stderr.print " -R R2          result of interest is R squared\n"
  $stderr.print " -R Q2          result of interest is Q squared\n"
  $stderr.print " -R AE50        result of interest is AE50\n"
  $stderr.print " -R AE95        result of interest is AE95\n"
# GH Removed. RMS50 and RMS95 are not calculated in svmfp_calibrate.sh based on the comment
#             from Suntara Cahya
#  $stderr.print " -R RMS50       result of interest is RMS50\n"
#  $stderr.print " -R RMS95       result of interest is RMS95\n"
  $stderr.print " -R cMSR90_mad    result of interest is cMSR (median average deviation) (Suntara)\n"
  $stderr.print " -R cMSR90_sd     result of interest is cMSR (standard deviation) (Suntara)\n"
  $stderr.print " -R cMSR95_mad    result of interest is cMSR (median average deviation) (Suntara)\n"
  $stderr.print " -R cMSR95_sd     result of interest is cMSR (standard deviation) (Suntara)\n"
  $stderr.print " -R rx=<rx>     for extracting the model type from the file name\n" if ($expert)
  $stderr.print " -C .           classification model (min fraction correct)\n" if ($expert)
  $stderr.print " -C mfc         classification (min fraction correct)\n" if ($expert)
  $stderr.print " -C afc         classification (ave fraction correct)\n" if ($expert)
# $stderr.print " -p precision   precision for outputs (default 3)\n"
  $stderr.print " -T <nt>        perform t test on the <nt> best fingerprints compared to best\n" if ($expert)
  $stderr.print " -D <fname>     file of raw values\n" if ($expert)
  $stderr.print " -B <fname>     write best fingerprint model to <fname>\n"
  $stderr.print " -x             normally the -B file contains the best fingerprint only, -x allows other forms\n" if $expert
  $stderr.print " -v             verbose output\n"
  exit(rc)
end

cl = IWCmdline.new("-v-R=s-C=s-p=ipos-r-col=ipos-RX=s-rx=s-T=ipos-D=s-expert-std-q-S=s-F=sfile-tar=sfile-B=s-x")

$expert = cl.option_present('expert')

if cl.unrecognised_options_encountered()
  $stderr.print "Unrecognised options encountered\n"
  usage(1)
end

verbose = cl.option_present('v')

if cl.option_present('tar')
  require 'archive/tar/minitar'
end

$display_mesage_for_missing_or_empty_files = ! cl.option_present('q')

$missing_or_empty_files_skipped = 0

class Next_File
  def initialize
    @argv = false
    @rx = false
    @dir = false
    @file = false
    @tarfile = false
  end

  def initialise(cl, a)
    if a.size > 0
      @argv = a
    elsif cl.option_present('F')
      fname = cl.value('F')
      @file = File.open(fname, mode='r')
      raise "Cannot open '#{fname}'" unless @file
    elsif cl.option_present('rx')
      @rx = Regexp.new(cl.value('rx'))
      @dir = Dir.open('.')
      raise "Cannot open current directory for reading" unless @dir
    elsif cl.option_present('tar')
      fname = cl.value('tar')
      @tarfile = Archive::Tar::Minitar.open(fname)
      raise "Cannot open tar file '#{fname}'" unless @tarfile
    else
      $stderr << "Must specify files to process via command line args, -rx, -F or -tar options\n"
      return false
    end

    return true
  end
  
  def file_exists(fname)
    if File.size?(fname)
      return true
    else
      $stderr << "Missing or empty file '#{fname}'\n" if $display_mesage_for_missing_or_empty_files
      $missing_or_empty_files_skipped += 1
      return false
    end
  end

  def each
    if @argv
      @argv.each do |fname|
        next unless file_exists(fname)
        yield fname, IO.read(fname)
      end
    elsif @dir
      @dir.each do |fname|
        next unless @rx.match(fname)
        next unless file_exists(fname)
        yield fname, IO.read(fname)
      end
    elsif @file
      @file.each do |fn|
        fn.chomp!
        fname = fn.gsub(/\$(\w+)/){ENV[$1]}
        next unless file_exists(fname.to_str)
        yield fname, IO.read(fname)
      end
    elsif @tarfile
      @tarfile.each do |entry|
        yield entry.name, entry.read
      end
    end
  end
end

next_file = Next_File.new

if ! next_file.initialise(cl, ARGV)
  $stderr << "Cannot initialise next file iterator\n"
  usage(1)
end

$output_precision = 3
if (cl.option_present('p'))
  $output_precision = cl.value('p')
end

column_to_process = false
if (cl.option_present('col'))
  column_to_process = cl.value('col')
  column_to_process -= 1
end

$show_std_dev = cl.option_present('std')

rx = Array.new     # an array of regular expressions
rxname = Array.new
reverse = Array.new

#Average absolute error 1.22571
#SVM_Predicted Unbiased Bsquared 0.674414

rxrx = Regexp.new('^rx=(\\S+)')

header_for_Dfile = "B2"

if (cl.option_present('R'))

  cl.values('R').each do |r|

    if ('B2' == r)
      rx.push(Regexp.new(' Bsquared '))
      rxname.push('B2')
      reverse.push(false)
    elsif ('R2' == r)
      rx.push(Regexp.new(' R2 '))
      rxname.push('R2')
      reverse.push(false)
    elsif ('Q2' == r)
      rx.push(Regexp.new(' Q2 '))
      rxname.push('Q2')
      reverse.push(false)
    elsif ('AE' == r)
      rx.push(Regexp.new('^Average absolute error '))
      rxname.push('AE')
      reverse.push(true)
    elsif ('RMSE' == r || 'RMS' == r)
      rx.push(Regexp.new('^RMS error '))
      rxname.push('RMS')
      reverse.push(true)
    elsif ('D' == r)
      rx.push(Regexp.new('^Distributional similarity'))
      rxname.push('DS')
      reverse.push(false)
    elsif ('PPV' == r)
      rx.push(Regexp.new(' PPV '))
      rxname.push('PPV')
      reverse.push(false)
    elsif ('AE50' == r)
      rx.push(Regexp.new(' AE50 '))
      rxname.push('AE50')
      reverse.push(true)
    elsif ('AE95' == r)
      rx.push(Regexp.new(' AE95 '))
      rxname.push('AE95')
      reverse.push(true)
    elsif ('cMSR90_mad' == r)
      rx.push(Regexp.new(' cMSR90_mad '))
      rxname.push('cMSR90_mad')
      reverse.push(true)
    elsif ('cMSR95_mad' == r)
      rx.push(Regexp.new(' cMSR95_mad '))
      rxname.push('cMSR95_mad')
      reverse.push(true)
    elsif ('cMSR90_sd' == r)
      rx.push(Regexp.new(' cMSR90_sd '))
      rxname.push('cMSR90_sd')
      reverse.push(true)
    elsif ('cMSR95_sd' == r)
      rx.push(Regexp.new(' cMSR95_sd '))
      rxname.push('cMSR95_sd')
      reverse.push(true)
    elsif ("bias" == r.downcase)
      rx.push(Regexp.new(' bias ', Regexp::IGNORECASE))
      rxname.push('bias')
      reverse.push(true)
    elsif ("all" == r.downcase)
      rx.push(Regexp.new(' Bsquared '))
      rxname.push('B2')
      reverse.push(false)
      rx.push(Regexp.new(' R2 '))
      rxname.push('R2')
      reverse.push(false)
      rx.push(Regexp.new(' Q2 '))
      rxname.push('Q2')
      reverse.push(false)
      rx.push(Regexp.new('^Average absolute error '))
      rxname.push('AE')
      reverse.push(true)
      rx.push(Regexp.new('^RMS error '))
      rxname.push('RMS')
      reverse.push(true)
      rx.push(Regexp.new(' AE50 '))
      rxname.push('AE50')
      reverse.push(true)
      rx.push(Regexp.new(' AE95 '))
      rxname.push('AE95')
      reverse.push(true)
      rx.push(Regexp.new(' cMSR90_sd '))
      rxname.push('cMSR90_sd')
      reverse.push(true)
      rx.push(Regexp.new(' cMSR95_sd '))
      rxname.push('cMSR95_sd')
      reverse.push(true)
      rx.push(Regexp.new(' cMSR90_mad '))
      rxname.push('cMSR90_mad')
      reverse.push(true)
      rx.push(Regexp.new(' cMSR95_mad '))
      rxname.push('cMSR95_mad')
      reverse.push(true)
    elsif (m = rxrx.match(r))
      rx.push(Regexp.new(m[0]))
      rxname.push(m[0])
      reverse.push(false)   # we really do not know
    else
      $stderr.print "Unrecognised -R qualifier '#{r}'\n"
      usage(5)
    end
  end
elsif (cl.option_present('C'))
  c = cl.value('C')

# m = MatchData.new
  class_rx = Regexp.new('CLASS:(\\S+)')

  if ('.' == c || 'mfc' == c)
    rx.push(Regexp.new('Min fraction correct'))
    rxname.push('MinFractionCorrect')
    reverse.push(false)
  elsif ('afc' == c)
    rx.push(Regexp.new('Average fraction correct'))
    rxname.push('AveFractionCorrect')
    reverse.push(false)
  elsif ('fc' == c)
    rx.push(Regexp.new('Fraction correct'))
    rxname.push('FractionCorrect')
    reverse.push(false)
  elsif (m = class_rx.match(c))
    c = m[1]
    rx.push(Regexp.new("instances of class '#{c}'"))
    rxname.push(c)
    reverse.push(false)
  else
    $stderr.print "Unrecognised -C qualifier '#{c}'\n"
    usage(5)
  end
else
  rx.push(Regexp.new(' Bsquared '))
  rxname.push('B2')
  reverse.push(false)
end

$stderr << "Defined #{rx.size} regular expressions\n" if verbose

if rx.size > 1 && ! cl.option_present('S')
  $stderr << "More than one sort criteria specified, must specify output file name stem\n"
  usage(1)
end

output_file_name_stem = cl.value('S')

class Set_of_Data
  def initialize
    @minval = false
    @maxval = false
    @sum = false
    @sum2 = false
    @values = Array.new
  end

  def n
    @values.size
  end

  def extra(v)
    if (0 == @values.size)
      @minval = v
      @maxval = v
      @sum = v
      @sum2 = v * v
    else
      @minval = v if (v < @minval)
      @maxval = v if (v > @maxval)
      @sum += v
      @sum2 += v * v
    end

    @values.push(v)
  end

  def average
    @sum / @values.size.to_f
  end

  def report (name, longest_name, os)
    os.printf("%-#{longest_name}s had #{@values.size} values between %.3f and %.3f", name, @minval, @maxval)

    n = @values.size.to_f
    ave = @sum / n
    if ($show_std_dev)
      var = (@sum2 - n * ave * ave) / (n - 1.0)
      os.printf(" ave %.#{$output_precision}f std %.#{$output_precision}f\n", ave, Math.sqrt(var))
    else
      os.printf(" ave %.#{$output_precision}f \n", ave)
    end
  end

  def write_prediction (ndx, os)
    os << @values[ndx]
  end
  def write_all_values (os)
    @values.each do |v|
      os << " #{v}"
    end
    os << "\n"
  end
end

class Model
  def initialize(s)
    @name = s
    @measures = Hash.new
    @reverse = false
  end

  def n(f)
    if ! @measures.has_key?(f)
      $stderr << "#@name has no key #{f}\n"
      return 0
    end
    @measures[f].n
  end
  def name
    @name
  end

  def extra(f, v)
    $st
    if (! @measures.has_key?(f))
      @measures[f] = Set_of_Data.new
#     $stderr << "Initialised data for #{f}\n"
    end

    @measures[f].extra(v)
  end

  def average(f)
    @measures[f].average
  end

  def report(f, longest_name, os)
    @measures[f].report(@name, longest_name, os)
  end

  def write_prediction(f, ndx, s)
    @measures[f].write_prediction(ndx, s)
  end

  def write_all_values(f, os)
    @measures[f].write_all_values(os)
  end
end

def which_rx_matches (rx, line)
  rx.each_index do |i|
    return i if rx[i].match(line)
  end

  return false
end

# This regexp is used to extract the actual name of the method from the file name

fname_rx = Regexp.new('[A-Z]\d+\\.(\\S+)\\.\d+')

if (cl.option_present('RX'))
  fname_rx = Regexp.new(cl.value('RX'))

  $stderr << "File name regular expression '#{fname_rx.source}'\n" if (verbose)
end

predictors = Hash.new

longest_name = 0


next_file.each do |fname, file_contents|
# $stderr << "Processing '#{fname}' containing #{contents.size} bytes\n"

  m = fname_rx.match(fname)

  if (! m)
    $stderr.print "Skipping non valid file name '#{fname}' #{fname_rx.source}\n"
    next
  end

  predictor = m[1]

  if ! predictor
    $stderr << "No captured value from '#{fname}' matched by '${fname_rx.source}', do you have a capture (..)\n"
    exit 1
  end

  longest_name = predictor.length if (predictor.length > longest_name)

  if (! predictors.has_key?(predictor))
    predictors[predictor] = Model.new(predictor)
  end

  lines = file_contents.split("\n")

  zvalue = false    # scope here for efficiency

  lines.each do |line|

    ndx = which_rx_matches(rx, line)
    next unless (ndx)

    f = line.split

    if (column_to_process)
      zvalue = f[column_to_process]
    else
      zvalue = f.pop
    end

    begin
      t = Float(zvalue)
    rescue
      $stderr << "Invalid float '#{zvalue}', for '#{rx[ndx].source}' set to zero\n"
      t = 0.0;
    end

#   $stderr.print "From '#{line.chomp}' get #{t}\n"

    predictors[predictor].extra(ndx, t)
  end
end

if ! $display_mesage_for_missing_or_empty_files && $missing_or_empty_files_skipped > 0
  $stderr << "Skipped #{$missing_or_empty_files_skipped} missing or empty files\n"
end

if 0 == predictors.size
  $stderr << "NO data\n"
  exit 1
end

def do_report(name_stem, m, rxname, i, longest_name)
  fname = "#{name_stem}.#{rxname}.dat"
  outp = File.open(fname, mode='w')
  raise "Cannot open '#{fname}'" unless outp
  m.report(i, longest_name, outp)
  outp.close
end

if (cl.option_present('r'))
  reverse.each_index do |i|
    reverse[i] = ! reverse[i]
  end
end

dfile_stem = false
dfile_stem = cl.value('D') if cl.option_present('D')

bfile_stem = false
bfile_stem = cl.value('B') if cl.option_present('B')

best_fp_only = true

best_fp_only = false if cl.option_present('x')

# We need a means of identifying those models that are fingerprint based

def looks_like_fingerprint s
  return false unless s
  return false unless s.length > 0
  return false if s.index("ARB")
#GH Only check the string starts with those substring
  return false if (0 == s.index("CUBIST")) 
  return false if (0 == s.index("SVMD"))
  return false if (0 == s.index("KNN"))

#  $stderr << "Valid fingerprint '#{s}'\n"
  return true
end

nt = cl.value('T')

if (1 == nt)
  $stderr << "Must specify more than one model to test with the -T option\n"
  usage(3);
end

(0...rx.size).each do |i|
  $stderr << "Processing predictor #{i} #{rxname[i]}\n" if verbose
  array_of_predictors = Array.new

  nsplit = 0

  predictors.each do |k, v|
    n = v.n(i)
    if (n > 0)
      array_of_predictors.push(v)
      nsplit = n if (n > nsplit)
    else
      $stderr.print "No #{rxname[i]} values for '#{v.name}', ignored\n"
    end
  end

#$stderr.print "Have #{array_of_predictors.size} different predictions\n"

  rvs = reverse[i]
  array_of_predictors.sort! {|a, b| 
    va = a.average(i)
    vb = b.average(i)

    if (rvs)
      tmp = va
      va = vb
      vb = tmp
    end

    va <=> vb
  }

  $stderr.print "After sort #{array_of_predictors.size} values\n" if verbose

  os = $stdout
  if rx.size > 1
    fname = "#{output_file_name_stem}.#{rxname[i]}.dat"
    os = File.open(fname, mode='w')
    raise "Cannot open '#{fname}'" unless os
  end

  array_of_predictors.each do |m|
    m.report(i, longest_name, os)
  end

  os.close if rx.size > 1

  if dfile_stem
    dfile = "#{dfile_stem}"
    dfile << i.to_s if rx.size > 1
    outp = File.open(dfile, mode='w')
    raise "Cannot open raw results file '#{dfile}'" unless outp

    outp << rx[i].source << "\n"
    predictors.each do |k, v|
      next unless (v.n(i) == nsplit)

      outp << k 
      v.write_all_values(i, outp)
    end
  end

  if bfile_stem
    fname = "#{bfile_stem}.#{rxname[i]}.dat"
    os = File.open(fname, mode='w')
    raise "Cannot open '#{fname}'" unless os

    if best_fp_only
#     $stdout << " have #{array_of_predictors.size} predictors\n"
#     $stdout << array_of_predictors.compact.select { |p| looks_like_fingerprint(p.name) }.size << "\n"
      os << array_of_predictors.compact.select { |p| looks_like_fingerprint(p.name) }.last.name().gsub(/-/, ' -').gsub(/^ /, "") << "\n"
    else
      os << array_of_predictors.last.name().gsub(/-/, ' -').gsub(/^ /, "") << "\n"
    end
  end

  next unless nt

  tmpdir = "."

  tmpdir = "/node/scratch" if (! File.writable?(tmpdir))
  tmpdir = "/tmp" if (! File.writable?(tmpdir))

  tmpfile = "#{tmpdir}/svmfp_summarise_T_test#{Process.pid}.dat"

  outp = File.open(tmpfile, mode='w')
  raise "Cannot open '#{tmpfile}'" unless (outp)

  nsplits = array_of_predictors[0].n(i)

  istart = array_of_predictors.size - nt - 1;

  outp << "N"
  (array_of_predictors.size - 1).downto(istart).each do |i|
    outp << " #{array_of_predictors[i].name}"
  end
  outp << "\n"

  (0...nsplits).each do |x|
    outp << "#{x}"
    (array_of_predictors.size - 1).downto(istart).each do |j|
      outp << ' '
      array_of_predictors[j].write_prediction(x, outp)
    end
    outp << "\n"
  end

  outp.close

  if (! File.size?(tmpfile))
    $stderr << "Did not create temporary file '#{tmpfile}'\n"
    exit 2
  end

  cmd = "test_t_test.sh -s 1 -c 2"
  (nt+1).downto(3).each do |i|
    cmd << " -c #{i}"
  end
  cmd << " #{tmpfile}"

  $stdout << "\n"    # just a line to separate results

  $stderr << "Executing '#{cmd}'\n" if (verbose)

  system(cmd)

  File.unlink(tmpfile)
end
