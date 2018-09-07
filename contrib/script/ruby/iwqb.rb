#! /usr/bin/env ruby

# $Id$


lillymol_home = ".."
require "#{lillymol_home}/ruby/lib/iwcmdline.rb"

stem = "iwqb#{Process.pid}"

$expert = false

def usage (rc)
  $stderr.print "Vastly simplified version of qblastem - no error checking, no nothing!\n"
  $stderr.print " -stem <stem>   stem to use for temporary files\n"
  $stderr.print " -m             allow multi line scripts\n"
  $stderr.print " -M <n>         group input into multi-job chunks, size <n>\n"
  $stderr.print " -qsub ... -qsub options passed to qsub\n"
  $stderr.print " -dir ...  -dir  directive(s) to be inserted into the shell script\n"
  $stderr.print " -j              join standard out and standard error\n"
  $stderr.print " -sync           wait for jobs to complete before exiting\n"
  $stderr.print " -N <name>       job name passed to qsub\n" if ($expert)
  $stderr.print " -expert         more options\n" unless ($expert)
  $stderr.print " -v             verbose output\n"
  exit(rc)
end

cl = IWCmdline.new("-v-expert-stem=s-m-M=ipos-cluster=s-qsub=close-dir=close-sync-N=s-j")

$expert = cl.option_present('expert')

if cl.unrecognised_options_encountered
  $stderr.print "Unrecognised options encountered\n"
  usage(1)
end

verbose = cl.option_present('v')

if 0 == ARGV.size
  $stderr.print "Insufficient arguments\n"
  usage(2)
end

if cl.option_present('stem')
  stem = cl.value('stem')
end

multi_line_scripts = cl.option_present('m')

user_specified_items_per_chunk = false

user_specified_items_per_chunk = cl.value('M') if (cl.option_present('M'))

command_file = ARGV[0]

unless FileTest.size?(command_file)
  $stderr.print "Missing or empty command file '#{command_file}'\n"
  exit(3)
end

jobs_created = 0

starts_with_hash = Regexp.new('^#')

master = "iwqb.master#{Process.pid}.sh"
m = File.open(master, mode='w')
raise "Cannot open master script '#{master}'" unless (m)

m.print "#! /bin/bash\n"
m.print "#\$ -cwd\n"
m.print "#\$ -S /bin/bash\n"

if cl.option_present('dir')
  cl.values('dir').each do |d|
    m.print "#\$ #{d}\n"
  end
end

m.print "hostname >&2\n"
m.print "uname=`uname`\n"
m.print "export LILLYMOL_HOME=#{lillymol_home}\n"
m.print "export PATH=$PATH:$LILLYMOL_HOME/bin/sh\n"

m.print "\n"

inp = File.open(command_file, mode="r")
raise "Cannot open '#{command_file}'" unless (inp)

scripts = Array.new

if multi_line_scripts
  rawdata = inp.readlines.join

  scripts = rawdata.split(/\n\|\n/)
elsif user_specified_items_per_chunk
  current_chunk = ''
  items_current_chunk = 0
  inp.each do |line|
    current_chunk << line
    items_current_chunk += 1
    if items_current_chunk >= user_specified_items_per_chunk
      scripts.push(current_chunk)
      current_chunk = ''
      items_current_chunk = 0
    end
  end
  scripts.push(current_chunk) if (current_chunk.length > 0)
else
  scripts = inp.readlines
end

scripts.each do |script|
  if 0 == script.chomp.length
    next
  end

  next if (starts_with_hash.match(script))    # skip comment lines
  next if (1 == script.length)    # just the newline

  jobs_created += 1

  m.print "function iwqb#{jobs_created}\n"
  m.print "{\n"
  m.print script
  m.print "\n" if (multi_line_scripts)
  m.print "}\n"
end

$stderr.print "Created #{jobs_created} jobs for submission\n"

if 0 == jobs_created
  $stderr.print "No jobs created!!\n"
  exit 3
end

m.print "eval iwqb${SGE_TASK_ID}\n"
m.close

system("chmod +x #{master}")

cmd = " qsub -cwd -b y -t 1-#{jobs_created} "

if cl.option_present('N')
  cmd << ' -N ' << cl.value('N')
end

if cl.option_present('qsub')
  cmd << cl.value('qsub')
end

if cl.option_present('j')
  cmd << ' -j yes'
end

if cl.option_present('sync')
 cmd << ' -sync yes'
end

cmd << " ./#{master}"

$stderr.print "Executing '#{cmd}'\n" if (verbose)
system (cmd)
