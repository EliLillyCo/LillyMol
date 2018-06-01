#! /usr/bin/env ruby
# $Id$

#tool_home = ENV['C3TK_BIN']
tool_home = ".."

require "#{tool_home}/ruby/lib/iwcmdline"

#GH Original
#~ cl = IWCmdline.new("-v-a=i-o=i-e=i-w=ipos-k-ks-dry-start=i-stop=i-do=s-stem=s-rx=s-suffix=s-dsc=sfile-qsub-cluster-cluster.seq-sync-qsubopt=close-file=sfile-col=s-echo-echon-f=s-rxdir=dir-subdir-s-array=s-sleep=ipos-sortrx-expert-parallel=ipos-l=s-V-j-q-submit=xfile-noeval")
cl = IWCmdline.new("-v-a=i-o=i-e=i-w=ipos-k-ks-dry-start=i-stop=i-do=s-stem=s-rx=s-suffix=s-dsc=sfile-qsub-cluster-cluster.seq-sync-qsubopt=close-file=sfile-col=s-echo-echon-f=s-rxdir=dir-subdir-s-array=s-sleep=ipos-sortrx-expert-l=s-V-j-q-submit=xfile-noeval")

verbose = cl.option_present('v')

$expert = cl.option_present('expert')

dry_run = cl.option_present('dry')

qsub = cl.option_present('qsub') || cl.option_present('cluster') || cl.option_present('qsubopt') || cl.option_present('sync') || cl.option_present('l') || cl.option_present('V')

doitfilename = "DOPATTERNTMP#{Process.pid}"
doitfile = false

#GH Parallel is removed
#~ if qsub || cl.option_present('parallel')
  #~ dry_run = true
  #~ doitfile = File.open(doitfilename, mode='w')
  #~ unless doitfile
    #~ $stderr.print "Cannot create temporary file for qsub submission\n"
    #~ exit 3
  #~ end
#~ end

fsub = '%'

sleeptime = 0    # sleep between issuing commands

iwqb = "#{tool_home}/ruby/iwqb.rb"

if cl.option_present('submit')
  iwqb = cl.value('submit')
end

def usage (rc)
  clusterdir = '/opt/ge'
  $stderr.print "Executes a series of commands based on a range\n"
  $stderr.print " -a <start>      range start (default 1)\n"
  $stderr.print " -o <stop>       range stop\n"
  $stderr.print " -e <step>       range step (default 1)\n"
  $stderr.print " -w <width>      insert leading 0's to make each number <width> wide\n"
  $stderr.print " -k              keep going after failure\n"
  $stderr.print " -ks             keep going after failure, and be silent about failures\n"
  $stderr.print " -dry            dry run, don't actually execute any commands\n"
  $stderr.print " -do n,n,n       rather than a sequence, just to specific numbers\n"
  $stderr.print " -stem iwsplit%.smi determine numbers by looking at sequence of files\n"
  $stderr.print " -array x,x,x    string items available in \#{array[i]} during numeric\n" if $expert
  $stderr.print " -rx <rx>        process all files in current dir that match <rx>\n"
  $stderr.print "                 if <rx> contains (grouping) that will be used for %\n"
  $stderr.print " -rxdir <dir>    scan for matching files in <dir> rather than '.'\n" if $expert
  $stderr.print " -s              when doing -rx processing, only do files of non zero size\n" if $expert
  $stderr.print " -dsc <fname>    process all the descriptor names in <fname>\n"
  $stderr.print " -file <fname>   process all the tokens in <fname> - one per line\n"
  $stderr.print " -col ...        column specification for -file, use '-col all' for whole record\n" if $expert
  $stderr.print " -subdir         process all subdirectories\n"
  $stderr.print " -f <char>       character to use instead of %\n"
  $stderr.print " -noeval         do NOT evaluate \#{..} constructs (useful for embedded invocations)\n" if $expert
  $stderr.print " -cluster        collect all commands in a file and submit to qsub (Cluster setup is required)\n" if (FileTest.directory?(clusterdir))
  $stderr.print " -V              with the -cluster option: pass the -V option to qsub (Cluster setup is required)\n" if (FileTest.directory?(clusterdir) && $expert)
  $stderr.print " -cluster.seq    collect all commands in a file and submit to qsub - if numeric sequence (Cluster setup is required)\n" if (FileTest.directory?(clusterdir) && $expert)
  $stderr.print " -sync           submit to qsub, wait for results (Cluster setup is required)\n" if (FileTest.directory?(clusterdir))
  $stderr.print " -qsubopt ... -qsubopt options passed to qsub (Cluster setup is required)\n" if $expert
  $stderr.print " -l ...          cluster resource specification(s) (Cluster setup is required)\n" if $expert
#~ $stderr.print " -parallel <j>   use gnu parallel to process <j> tasks at a time\n" if $expert
  $stderr.print " -submit <fname> specify your own job submitter (default iwqb.rb) (Cluster setup is required)\n" if $expert
  $stderr.print " -echo           echo % before each command\n"
  $stderr.print " -sleep <sec>    sleep between commands\n" if $expert
  $stderr.print " -expert         more options\n" unless ($expert)
  $stderr.print " -v              verbose output\n"
  exit(rc)
end

if cl.option_present('o')
  istop = cl.value('o')
elsif cl.option_present('stop')
  istop = cl.value('stop')
elsif cl.option_present('do')
  true
elsif cl.option_present('array')
  true
elsif cl.option_present('stem')
  true
elsif cl.option_present('rx')
  true
elsif cl.option_present('suffix')
  true
elsif cl.option_present('dsc')
  true
elsif cl.option_present('file')
  true
elsif cl.option_present('subdir')
  true
else
  $stderr.print "Must specify stop value via the -o option\n"
  usage (3)
end

sleeptime = cl.value('sleep') if (cl.option_present('sleep'))

istart = 1

if cl.option_present('a')
  istart = cl.value('a')
end

if cl.option_present('start')
  istart = cl.value('start')
end

istep = 1

if cl.option_present('e')
  istep = cl.value('e')
end

if cl.option_present('w')
  width = cl.value('w')
else
  width = 0
end

if cl.option_present('f')
  fsub = cl.value('f')
end

items_to_do = Array.new

numeric_sequence = false

suppress_stderr = cl.option_present('q')

array = Array.new

# When doing the -rx option, if the regexp comes in
# with (groupings) then this array will be populated and
# will be what is substituted

matchdata = Array.new

if cl.option_present('do')
  if 1 == cl.values('do').size && '-' == cl.value('do')
    $stdin.each do |d|
      items_to_do.push(d.chomp)
    end
  else
    cl.values('do').each do |d|
      f = d.split(/,/)
      f.each do |i|
        items_to_do.push(i)
      end
    end
  end
elsif cl.option_present('array')
  ndx = 0
  cl.values('array').each do |a|
    f = a.split(/,/)
    f.each do |i|
      array.push(i)
      items_to_do.push(ndx)
      ndx += 1
    end
  end
elsif cl.option_present('stem')
  stem = cl.value('stem')

  unless FileTest.size?(stem.gsub(fsub, '1'))
    $stderr.print "No files of form '#{stem}' present\n"
    exit 3
  end

  istart = 1
  istop = 1
  while FileTest.size?(stem.gsub(fsub, istop.to_s))
    items_to_do.push(istop)
    istop += 1
  end
elsif cl.option_present('rx')
  s = cl.value('rx')
  rx = Regexp.new(s)

  only_match_non_zero_file_size = cl.option_present('s')

  dir = '.'
  dir = cl.value('rxdir') if (cl.option_present('rxdir'))

  $stderr.print "Scanning directory '#{dir}' for '#{s}' files\n" if (verbose)

  d = Dir.open(dir)
  raise "Cannot open current directory '#{dir}'" unless(d)
  d.each do |fname|
    m = (rx.match(fname))
    next unless (m)
    next if (only_match_non_zero_file_size && ! File.size?(fname))
    items_to_do.push(fname)
    if m.size > 1
      matchdata.push(m[1])
    end
  end
  d.close

  matchdata.sort! {|a, b| a <=> b} if cl.option_present('sortrx')
elsif cl.option_present('suffix')
  s = cl.value('suffix')
  rx = Regexp.new("(\\S+)\\.#{s}$")

  only_match_non_zero_file_size = cl.option_present('s')

  dir = '.'
  dir = cl.value('rxdir') if (cl.option_present('rxdir'))

  $stderr.print "Scanning directory '#{dir}' for '#{s}' files\n" if (verbose)

  d = Dir.open(dir)
  raise "Cannot open current directory '#{dir}'" unless(d)
  d.each do |fname|
    m = (rx.match(fname))
    next unless (m)
    next if (only_match_non_zero_file_size && ! File.size?(fname))
    items_to_do.push(fname)
    if m.size > 1
      matchdata.push(m[1])
    end
  end
  d.close

  matchdata.sort! {|a, b| a <=> b} if cl.option_present('sortrx')
elsif cl.option_present('dsc')
  f = cl.value('dsc')
  inp = File.open(f, mode='r')
  raise "Cannot open -dsc file '#{f}'" unless inp
  hdr = inp.gets
  inp.close
  f = hdr.split
  raise "Cannot read header from '#{f}'" unless (f.size > 1)
  f.shift
  f.each do |d|
    items_to_do.push(d)
  end
elsif cl.option_present('subdir')
  dname = '.'
  dname = cl.value('rxdir') if (cl.option_present('rxdir'))

  d = Dir.open(dname)
  raise "Cannot open directory '#{dname}'" unless (d)

  d.each do |fname|
    next unless (File.directory?(fname))
    next if (/^\./ =~ fname)
    items_to_do.push(fname)
  end
elsif cl.option_present('file')
  f = cl.value('file')
  inp = File.open(f, mode='r')
  raise "Cannot open -file file '#{f}'" unless inp
  if ! cl.option_present('col')
    inp.each do |line|
      ff = line.split
      items_to_do.push(ff[0])
    end
  elsif 'ALL' == cl.value('col') || 'all' == cl.value('col')
    inp.each do |line|
      items_to_do.push(line.chomp)
    end
  else
    columns_to_process = Hash.new
    cl.values('col').each do |sc|
      c = sc.to_i
      columns_to_process[c-1] = true
    end

    inp.each do |line|
      f = line.split
      todo = ''
      f.each_index do |i|
        next unless columns_to_process[i]
        todo << ' ' if (todo.length > 0)
        todo << f[i]
      end

      items_to_do.push(todo) if (todo.length > 0)
    end
  end
  inp.close
else
  if istop < istart && istep > 0
    $stderr.print "Stop value '#{istop}' less than start value '#{istart}'\n"
  end

  istart.step(istop, istep) do |i|
    items_to_do.push(i)
  end

  numeric_sequence = true if (1 == istart && 1 == istep)
end

keep_going_after_failure = cl.option_present('k')
report_failures = true

if cl.option_present('ks')
  keep_going_after_failure = true
  report_failures = false
end

if 0 == ARGV.size
  $stderr.print "No command to issue\n"
  usage(1)
end

cmd = ARGV.join(' ')

if (! numeric_sequence) && (0.class != items_to_do[0].class)
  as_numeric = Array.new
  digit_rx = Regexp.new('^\\d+$')
  numeric_sequence = true
  items_to_do.each do |i|
    unless digit_rx.match(i)
      numeric_sequence = false
      break
    end

    j = i.to_i

    if j <= 0
      numeric_sequence = false
      break
    end

    as_numeric[j] = true
  end

  if numeric_sequence   # make sure no gaps
    as_numeric.each_with_index do |i, v|
      unless v
        numeric_sequence = false
        break
      end
    end
  end
end

qsub_options = ''

if cl.option_present('qsubopt')
  qsub_options = cl.value('qsubopt')
end

if cl.option_present('V')
  qsub_options << ' -V'
end

if cl.option_present('l')
  cl.values('l').each do |l|
    qsub_options << " -l #{l}"
  end
end

if numeric_sequence && cl.option_present('cluster.seq')
  cmd.gsub!('%', '${SGE_TASK_ID}')
  qsub_file = "dptrntmp#{Process.pid}"
  outp = File.open(qsub_file, mode='w')
  raise "Cannot open '#{qsub_file}'" unless outp
  outp.print "#! /bin/bash\n"
  outp.print "#\$ -cwd\n"
  outp.print "#\$ -S /bin/bash\n"
  outp.print("#{cmd}\n")
  outp.close
  if cl.option_present('dry')
    $stderr << "Dry run specified, job not submitted, commands in '#{qsub_file}'\n"
    exit(0)
  end
  cmd = "/etc/cluster-setup.sh && qsub -t 1-#{items_to_do.size}"
  cmd << " #{qsub_options}" if qsub_options.length > 0
  cmd << " #{qsub_file}"
  $stderr << "Executing '#{cmd}'\n" if verbose
  system(cmd)
  exit
end

number_rx = Regexp.new('^\\d+$')

ruby_rx = Regexp.new('#\{.+\}')

if cl.option_present('noeval')
  ruby_rx = Regexp.new("^AALKJFDHIUYRJNBASDFUASDASD$")     # something very unlikely
end

echo_percent = cl.option_present('echo')
echo_percent_newline = cl.option_present('echon')

failures = 0

(0...items_to_do.size).each do |ndx|
  i = items_to_do[ndx]

  mycmd = "#{cmd}"

# $stderr << "Processing '#{mycmd}', item '#{i}'\n"

  percent = ''

  if ruby_rx.match(mycmd)
    t = "f = \"#{mycmd}\""
#   $stderr << "T is '#{t}', i = #{i}\n"
    if matchdata.size > 0
      t.gsub!(fsub, matchdata[ndx]) if (matchdata.size > 0)
    else
      t.gsub!(fsub, i.to_s)
    end
#   $stderr << "Before eval '#{t}'\n"
    mycmd = eval(t)
  elsif matchdata.size > 0
    mycmd.gsub!(fsub, matchdata[ndx])
    percent = matchdata[ndx]
  elsif number_rx.match(i.to_s)
#   $stderr.print "'#{i}' is numeric\n"
    tmp = i.to_s
    if width > 0
      tmp = sprintf("%0#{width}d", i)
    end

    mycmd.gsub!(fsub, tmp)
    percent = tmp
  else
    mycmd.gsub!(fsub, i)
    percent = "#{i}"
  end

# $stderr << "Before substitutions '#{mycmd}'\n"
  mycmd.sub!(/dquote(.+)dquote/,'"\1"')
  mycmd.sub!(/squote(.+)squote/,"'\\1'")
  mycmd.gsub!(/redirect/, '>')

# $stderr << "After quote replacement '#{mycmd}'\n"

  if doitfile
    doitfile.print "#{mycmd}\n"
    next
  end

  if dry_run
    $stdout.print "#{mycmd}\n"
  elsif verbose
    $stderr.print "#{mycmd}\n"
  end

  next if dry_run

  $stdout.print "#{percent} " if (echo_percent)
  $stdout.print "#{percent}\n" if (echo_percent_newline)

  if suppress_stderr
    mycmd << " 2> /dev/null"
  end

  unless system(mycmd)
    $stderr.print "'#{mycmd}' failed\n" if (report_failures)
    exit(1) unless keep_going_after_failure
    failures += 1
  end

  sleep(sleeptime) if (sleeptime > 0)
end

$stderr << "#{failures} of #{items_to_do.size} commands failed\n" if ((verbose || failures > 0) && report_failures)

if doitfile
  doitfile.close
  unless File.size?(doitfilename)
    $stderr << "No command file created '#{doitfilename}', nothing to do\n"
    File.unlink(doitfilename) if (File.exist?(doitfilename))
    exit
  end

#GH Original
  #~ if cl.option_present('parallel')
    #~ j = cl.value('parallel')
    #~ cmd = "#{ianhome}/p/parallel -j #{j} < #{doitfilename}" # TODO: we may want to relocate parallel to a lib/p path instead
  #~ else
    #~ cmd = "#{iwqb} "
    #~ cmd << "-qsub #{qsub_options} -qsub " if qsub_options.length > 0
    #~ if cl.option_present('sync')
      #~ cmd << '-sync '
    #~ end
    #~ if cl.option_present('j')
      #~ cmd << ' -j'
    #~ end
    #~ cmd << "#{doitfilename}"
  #~ end
    cmd = "#{iwqb} "
    cmd << "-qsub #{qsub_options} -qsub " if qsub_options.length > 0
    if cl.option_present('sync')
      cmd << '-sync '
    end
    if cl.option_present('j')
      cmd << ' -j'
    end
    cmd << "#{doitfilename}"

  if cl.option_present('dry')
    $stderr << "Dry run specified, job not submitted, would have executed '#{cmd}'\n"
    exit(0)
  end

  $stderr.print "Executing '#{cmd}'\n" if (verbose)
  system(cmd)

#GH Original
  #~ if cl.option_present('sync') || cl.option_present('parallel')
    #~ File.unlink(doitfilename)
  #~ end
  if cl.option_present('sync')
    File.unlink(doitfilename)
  end
  
end
