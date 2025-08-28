#!/usr/bin/env ruby

ianhome = ENV['LILLYMOL_HOME']


require_relative('lib/iwcmdline')
require_relative('lib/accumulator')

$expert = false

def usage (rc)
  $stderr.print "Records the average line number of a pattern in a file\n"
  $stderr.print " -f <records>   skip the first <records> records in the input\n"
  $stderr.print " -rx <regexp>   regexp to process\n"
  $stderr.print " -RX <fname>    read regular expressions from <fname>\n"
  $stderr.print " -d <dname>     name of descriptor to process\n"
  $stderr.print " -e <median>    request various percentile measurements (50 == median)\n"
  $stderr.print " -expert        more options\n" unless ($expert)
  $stderr.print " -v             verbose output\n"
  exit(rc)
end

cl = IWCmdline.new("-v-expert-rx=s-e=ipos-RX=sfile-f-d=s-f=ipos")

if (cl.option_present('expert'))
  $expert = true
end

if cl.unrecognised_options_encountered()
  $stderr.print "Unrecognised options encountered\n"
  usage(1)
end

verbose = cl.option_present('v')

multiple_rx_within_a_file = true

if cl.option_present('f')
  multiple_rx_within_a_file = false
end

dname = false
if cl.option_present('d')
  dname = cl.value('d')
end

if (cl.option_present('rx'))
  true
elsif cl.option_present('RX')
  true
else
  $stderr << "Must specify one or more regular expressions via the -rx option\n"
  usage(2)
end

$medians = Array.new
if cl.option_present('e')
  cl.values('e').each do |e|
    $medians.push(e)
  end
else
  $medians.push(50)
end

if (0 == ARGV.size)
  $stderr.print "Insufficient arguments\n"
  usage(2)
end

class Regexp_and_Count
  def initialize(s)
    @rx = Regexp.new(s)
    @acc = Accumulator.new
    @lines = Array.new
    @dname = false
    @descriptor_column = -1
  end

  def set_descriptor_name(s)
    @dname = s
  end

  def extra(line, lineno)
    if @dname
      f = line.chomp.split
      if @descriptor_column < 0
        @descriptor_column = f.index(@dname)
        raise "No #{@dname} in #{line.chomp}" unless @descriptor_column
        $stderr << "Descriptor '#{@dname}' found in column #{@descriptor_column}\n"
      end

      if @rx.match(f[@descriptor_column])
        @acc.extra(lineno)
        @lines.push(lineno)
      end
    elsif @rx.match(line)
      @acc.extra(lineno)
      @lines.push(lineno)
    end
  end

  def report(s, lines_in_file)
    s << "Pattern '#{@rx.source}' "
    if (0 == @acc.n)
      s << "not matched\n"
    elsif 1 == @lines.size
      median_fraction = @acc.minval.to_f / lines_in_file.to_f
      s.printf("matched 1 records, first #{@acc.minval} last #{@acc.maxval} of #{lines_in_file}, ave %.2f (%f), median #{@acc.maxval} fraction %.2f\n", @acc.average, 1, median_fraction)
    else
      s.printf("matched #{@acc.n} records, first #{@acc.minval} last #{@acc.maxval} of #{lines_in_file}, ave %.2f (%f),", @acc.average, @acc.average/lines_in_file)
      @lines.sort!
      $medians.each do |e|
        median = @lines[e * @lines.size / 100]
        median_fraction = median.to_f / lines_in_file.to_f
        s.printf(" #{e}%s #{median} fraction %.2f", "%", median_fraction)
      end
      s.printf("\n")
    end
  end

  def average
    return 0.0 if (0 == @acc.n)
    return @acc.average
  end
end

rc = Array.new

if cl.option_present('rx')
  cl.values('rx').each do |r|
    t = Regexp_and_Count.new(r)
    rc.push(t)
  end
elsif cl.option_present('RX')
  fname = cl.value('RX')
  rxs = IO.readlines(fname)
  raise "Cannot read regular expressions from '#{fname}'" unless rxs.size > 0
  rxs.each do |line|
    rc.push(Regexp_and_Count.new(line.chomp))
  end
end

if dname
  rc.each do |x|
    x.set_descriptor_name(dname)
  end
end

skip_first_lines = if cl.option_present('f')
                     cl.value('f')
                   else
                     0
                   end


# not safe to use ARGF.lineno - breaks on sshfs

lineno = 1

if multiple_rx_within_a_file
  lineno_prev_file = -1

  ARGF.each do |line|
    rc.each do |r|
      lineno += 1
      next if lineno <= skip_first_lines
      r.extra(line, lineno)
    end
  end

  if lineno_prev_file > 0 && lineno_prev_file != lineno
    $stderr << "Warning, line number mismatch #{lineno_prev_file} vs #{lineno} after processing '#{fname}'\n"
  end

  lineno_prev_file = lineno
else
  ARGV.each do |fname|
    inp = File.open(fname, mode='r')
    raise "Cannot open '#{fname}'" unless inp
    lineno = 0
    inp.each do |line|
      lineno += 1
      next if lineno <= skip_first_lines
      rc.each do |r|
        r.extra(line, lineno)
      end
    end
  end
end

rc.sort!{|a, b| a.average <=> b.average}

rc.each do |r|
  r.report($stdout, lineno)
end
