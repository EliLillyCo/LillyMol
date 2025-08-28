#!/usr/bin/env ruby
ianhome = ENV['LILLYMOL_HOME']

require_relative('lib/iwcmdline_v2')

$expert = false

def usage (rc)
  $stderr.print "Multi-threaded version of tsubstructure\n"
  $stderr.print " -thr <nthreads>  number of threads to use\n"
  $stderr.print " -nj              do NOT join the output files, leave in split form\n"
  $stderr.print " -tmpdir <dir>    directory for temporary files\n" if ($expert)
  $stderr.print " -expert          more options\n" unless ($expert)
  $stderr.print " -v               verbose output\n"
  exit(rc)
end

cl = IWCmdlineV2.new("-v-expert-thr=ipos-tmpdir=s-m=s-n=s-tsub=xfile-nj-R=s-a")

if (cl.option_present('expert'))
  $expert = true
end

verbose = cl.option_present('v')

tsubstructure = 'tsubstructure.sh'

tsubstructure = cl.value('tsub') if (cl.option_present('tsub'))

if ARGV.empty?
  $stderr.print "Insufficient arguments\n"
  usage(2)
end

mfile = cl.value('m')
nfile = cl.value('n')
rfile = cl.value('R')

afile = cl.option_present('a') ? "afile" : false

if (mfile || nfile || rfile)
  true
else
  $stderr << "No output specifiers (-m or -n), cannot use threaded version\n"
  usage(2)
end

common_args = ""

if (cl.option_count('m') > 1)
  cl.values('m').each do |m|
    if ('QDT' == m)
      common_args << " -m QDT"
    else
      mfile = m
    end
  end
end

nthreads = if cl.option_present('thr')
             cl.value('thr')
           else
             2
           end

input_file = ARGV.pop

unless File.size?(input_file)
  $stderr << "Missing or empty input file '#{input_file}'\n"
  exit 2
end

# We prefer /node/scratch to /tmp as it is usually larger
tmpdir = '/tmp'
if FileTest.directory?('/node/scratch')
  tmpdir = '/node/scratch'
end

if (cl.option_present('tmpdir'))
  tmpdir = cl.value('tmpdir')

  Dir.mkdir(tmpdir) unless (FileTest.directory?(tmpdir))
end

file_size = File.size?(input_file)

inp = File.open(input_file, mode='r')
raise "Cannot open '#{input_file}'" unless inp

chunk_size = file_size / nthreads

if chunk_size < 10
  $stderr << "File '#{input_file}' too small for #{nthreads} threads\n"
  exit 2
end

offset = Array.new
offset.push(0)

(1...nthreads).each do |i|
  o = i * chunk_size

  if ( ! inp.seek(o, IO::SEEK_SET))
      raise "Cannot seek to #{o} in '#{input_file}'\n"
  end
  line = inp.gets
  offset.push(inp.pos)
end

if (verbose)
  $stderr << "'#{input_file}', size #{file_size} divided into #{nthreads} chunks, each #{chunk_size} bytes\n"
  offset.each do |o|
    $stderr << o << "\n"
  end
end

ntmp = false
mtmp = false
rtmp = false
atype = false

rejoin = true

rejoin = false if (cl.option_present('nj'))

if rejoin
  ntmp = "#{tmpdir}/#{File.basename(nfile)}#{Process.pid}" if (nfile)
  mtmp = "#{tmpdir}/#{File.basename(mfile)}#{Process.pid}" if (mfile)
  rtmp = "#{tmpdir}/#{File.basename(rfile)}#{Process.pid}" if (rfile)
  atmp = "#{tmpdir}/#{File.basename(afile)}#{Process.pid}" if (afile)
else
  ntmp = "#{tmpdir}/#{File.basename(nfile)}" if (nfile)
  mtmp = "#{tmpdir}/#{File.basename(mfile)}" if (mfile)
  rtmp = "#{tmpdir}/#{File.basename(rfile)}" if (rfile)
  atmp = "#{tmpdir}/#{File.basename(afile)}" if (afile)
end

special_chars = Regexp.new('[$\[!><()@~]')

ARGV.each do |a|
  if (special_chars.match(a))
    common_args << " '#{a}'"
  else
    common_args << " #{a}"
  end
end

afiles = Array.new

pids = Array.new

(1...nthreads).each do |i|
  cmd = "#{tsubstructure} -A D #{common_args}"

  cmd << " -m #{mtmp}#{i}" if (mtmp)
  cmd << " -n #{ntmp}#{i}" if (ntmp)
  cmd << " -R #{rtmp}#{i}" if (rtmp)
  cmd << " -a " if (atmp)

  cmd << " -i seek=#{offset[i]}"
  cmd << " -i stop=#{offset[i+1]}" if offset[i+1]
    
  cmd << " #{input_file}"

  if (atmp)
    afile = "#{atmp}#{i}"
    afiles.push(afile)
    cmd << " > #{afile}"
  end

  $stderr << "Executing '#{cmd}'\n" if (verbose)

  p = fork {
    exec(cmd)
  }

  pids.push(p)
end

# Now do thread zero here

rfiles = Array.new

cmd = "#{tsubstructure} -A D #{common_args}"

if (! rejoin)
  cmd << " -m #{mtmp}0" if (mtmp)
  cmd << " -n #{ntmp}0" if (ntmp)
  cmd << " -R #{rtmp}0" if (rtmp)
  cmd << " -a" if atmp
else
  cmd << " -m #{mfile}" if (mtmp)
  cmd << " -n #{nfile}" if (ntmp)
  cmd << " -a" if atmp
  if (rtmp)
    cmd << " -R #{rfile}0"
    rfiles.push("#{rfile}0")
  end
end

cmd << " -i stop=#{offset[1]} #{input_file}"

if (atmp)
  afile = "#{atmp}0"
  cmd << " > #{afile}"
  afiles.push(afile)
end

$stderr << "Thread zero '#{cmd}'\n" if verbose
system(cmd)

cmdm = "cat" if (mfile)
cmdn = "cat" if (nfile)

mfiles_present = false
nfiles_present = false

(1...nthreads).each do |i|

  p = pids[i - 1]
# $stderr << "Waiting for #{p}\n"

  Process.wait(p)

  if (mtmp)
    mm = "#{mtmp}#{i}.smi"

    if (File.exist?(mm) && File.size(mm) > 0)
      cmdm << " #{mm}"
      mfiles_present  = true
    end
  end

  if (ntmp)
    nn = "#{ntmp}#{i}.smi"
    if (File.exist?(nn) && File.size(nn) > 0)
      cmdn << " #{nn}"
      nfiles_present  = true
    end
  end

  if (rtmp)
    rr = "#{rtmp}#{i}"
    if (File.exist?(rr) && File.size?(rr) > 0)
      rfiles.push(rr)
    end
  end
end

exit 0 if (! rejoin)

if (mfiles_present)
  cmdm << " >> #{mfile}.smi"
  system(cmdm)

  (1..nthreads).each do |i|
    File.unlink("#{mtmp}#{i}.smi") if (File.exist?("#{mtmp}#{i}.smi"))
  end
end

if (nfiles_present)
  cmdn << " >> #{nfile}.smi"
  system(cmdn)

  (1..nthreads).each do |i|
    File.unlink("#{ntmp}#{i}.smi") if (File.exist?("#{ntmp}#{i}.smi"))
  end
end

if (rfiles.size > 0)
  outp = File.open(rfile, mode='w')
  raise "Cannot open -R output file '#{rfile}'" unless outp
  count = Hash.new(0)
  rfiles.each do |fname|
    lines = File.readlines(fname)
    lines.each do |line|
      f = line.chomp.split
      count[f[1]] += f[0].to_i
#     $stderr << "Count for '#{f[1]} updated to #{count[f[1]]} from '#{f[0]}'\n"
    end
    File.unlink(fname)
  end
  count.each do |k, v|
    outp << "#{v} #{k}\n"
  end
end

if (afiles.size > 0)
  descriptor_file_cat = 'descriptor_file_cat.sh'
  cmd = "#{descriptor_file_cat}"
  afiles.each do |afile|
    if ( File.size?(afile))
      cmd << " #{afile}"
    else
      $stderr << "Missing or empty -a file '#{afile}', skipping\n"
    end
  end
  system(cmd)
  afiles.each do |afile|
    File.unlink(afile)
  end
end
