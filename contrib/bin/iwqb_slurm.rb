#! /usr/bin/env ruby
# $Id$
#
# iwqb with optional Slurm support.
#
# Default behaviour remains Grid Engine/qsub, matching the historic iwqb.rb.
# Use -slurm to generate and submit a Slurm array job with sbatch.

require 'shellwords'

lillymol_home = ".."

require_relative "lib/iwcmdline.rb"

stem = "iwqb#{Process.pid}"

$expert = false

def usage(rc)
  $stderr.print "Vastly simplified version of qblastem - no error checking, no nothing!\n"
  $stderr.print " -stem <stem> stem to use for temporary files\n"
  $stderr.print " -m allow multi line scripts\n"
  $stderr.print " -M <n> group input into multi-job chunks, size <n>\n"
  $stderr.print "\n"
  $stderr.print "Grid Engine, default scheduler:\n"
  $stderr.print " -qsub ... -qsub options passed to qsub\n"
  $stderr.print " -dir ... -dir directive(s) to be inserted into the shell script\n"
  $stderr.print "\n"
  $stderr.print "Slurm:\n"
  $stderr.print " -slurm use sbatch/Slurm rather than qsub/Grid Engine\n"
  $stderr.print " -sbatch ... -sbatch options passed to sbatch\n"
  $stderr.print " -sdir ... -sdir directive(s) to be inserted as #SBATCH lines\n"
  $stderr.print " -throttle <n> limit simultaneously running Slurm array tasks, --array=1-N%<n>\n"
  $stderr.print "\n"
  $stderr.print " -j join standard out and standard error\n"
  $stderr.print " -sync wait for jobs to complete before exiting\n"
  $stderr.print " -N <name> job name passed to qsub/sbatch\n" if ($expert)
  $stderr.print " -expert more options\n" unless ($expert)
  $stderr.print " -v verbose output\n"
  exit(rc)
end

# Keep all original options. Add Slurm-specific options.
# The -tp ... -tp and similar quote-hell constructs are intentionally not part
# of iwqb and are not implemented here.
cl = IWCmdline.new("-v-expert-stem=s-m-M=ipos-cluster=s-qsub=close-dir=close-sync-N=s-j-slurm-sbatch=close-sdir=close-throttle=ipos")
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

use_slurm = cl.option_present('slurm')

if use_slurm && cl.option_present('qsub')
  $stderr.print "Warning: -qsub options ignored in -slurm mode\n"
end

if use_slurm && cl.option_present('dir')
  $stderr.print "Warning: -dir directives ignored in -slurm mode; use -sdir instead\n"
end

if ! use_slurm && cl.option_present('sbatch')
  $stderr.print "Warning: -sbatch options ignored unless -slurm is specified\n"
end

if ! use_slurm && cl.option_present('sdir')
  $stderr.print "Warning: -sdir directives ignored unless -slurm is specified\n"
end

if ! use_slurm && cl.option_present('throttle')
  $stderr.print "Warning: -throttle ignored unless -slurm is specified\n"
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

# Preserve the historic master script naming convention. The -stem option was
# present in the original script, but the original master file name did not use
# it; keep that behaviour for compatibility.
master = "iwqb.master#{Process.pid}.sh"

m = File.open(master, mode='w')
raise "Cannot open master script '#{master}'" unless (m)

m.print "#! /bin/bash\n"

if use_slurm
  # Slurm directives must appear before any non-comment/non-blank shell code.
  m.print "#SBATCH --chdir=.\n"
  m.print "#SBATCH --output=slurm-%A_%a.out\n"
  if cl.option_present('j')
    m.print "#SBATCH --error=slurm-%A_%a.out\n"
  else
    m.print "#SBATCH --error=slurm-%A_%a.err\n"
  end
  if cl.option_present('sdir')
    cl.values('sdir').each do |d|
      m.print "#SBATCH #{d}\n"
    end
  end
else
  m.print "#$ -cwd\n"
  m.print "#$ -S /bin/bash\n"
  if cl.option_present('dir')
    cl.values('dir').each do |d|
      m.print "#$ #{d}\n"
    end
  end
end

m.print "hostname >&2\n"
m.print "uname=`uname`\n"
m.print "export LILLYMOL_HOME=#{lillymol_home}\n"
m.print "export PATH=$PATH:$LILLYMOL_HOME/bin\n"
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
inp.close

scripts.each do |script|
  if 0 == script.chomp.length
    next
  end

  next if (starts_with_hash.match(script)) # skip comment lines
  next if (1 == script.length) # just the newline

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
  m.close
  exit 3
end

# Modernised dispatch: no eval needed. Expanding the function name into a word
# and invoking it is sufficient in bash.
if use_slurm
  m.print "task_id=\"${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID not set}\"\n"
else
  m.print "task_id=\"${SGE_TASK_ID:?SGE_TASK_ID not set}\"\n"
end
m.print "fn=\"iwqb${task_id}\"\n"
m.print '"$fn"' + "\n"

m.close

system("chmod +x #{Shellwords.escape(master)}")

if use_slurm
  array_spec = "1-#{jobs_created}"
  if cl.option_present('throttle')
    array_spec << "%#{cl.value('throttle')}"
  end

  cmd = "sbatch --array=#{Shellwords.escape(array_spec)}"

  if cl.option_present('N')
    cmd << ' --job-name=' << Shellwords.escape(cl.value('N'))
  end

  if cl.option_present('sbatch')
    cmd << ' ' << cl.value('sbatch')
  end

  if cl.option_present('sync')
    cmd << ' --wait'
  end

  cmd << " ./#{Shellwords.escape(master)}"
else
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
end

$stderr.print "Executing '#{cmd}'\n" if (verbose)

system(cmd)
