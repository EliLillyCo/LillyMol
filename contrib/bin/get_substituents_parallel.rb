#!/usr/bin/env ruby

require_relative "lib/iwcmdline_v2.rb"
require_relative "lib/parallel_support.rb"

$expert = false

def usage(rc)
  $stderr << "Multi-threaded version of get_substituents\n"
  $stderr << "get_substituents_parallel -thr 16 -S out ... file1\n"
  $stderr << "Note that only one input file can be processed\n";
  $stderr << " -S <stem>        output stem - mandatory\n"
  $stderr << " -thr <nthreads>  number of threads to use\n"
  $stderr << " -nj              do NOT join the output files, leave in split form\n"
  $stderr << " -get_substituents <exe>  get_substituents executable to use (default get_substituents.sh)\n"
  $stderr << " -tmpdir <dir>    directory for temporary files\n" if ($expert)
  $stderr << " -log <stem>      redirect stderr to individual log files - which are removed\n"
  $stderr << "                  you probably want this to avoid screens full of warning messages\n"
  $stderr << "                  but make sure your reaction is well behaved and that warnings can actually be ignored\n";
  $stderr << " -keeplog         do NOT remove the logfiles, passes -v to the underlying executable\n"
  $stderr << " -v               verbose output\n"
  $stderr << "All other options are passed directly to get_substituents\n"

  exit(rc)
end

def get_substituents_parallel

  cl = IWCmdlineV2.new("-v-thr=ipos-tmpdir=s-S=s-get_substituents=xfile-nj-log=s-keeplog")
  
  verbose = cl.option_present('v')
  
  get_substituents = if cl.option_present('get_substituents') 
            cl.value('get_substituents')
          else
            'get_substituents.sh'
          end
  
  if ARGV.empty?
    $stderr << "Insufficient arguments\n"
    usage(2)
  end
  
  unless cl.option_present('S')
    $stderr << "Must specify output file name via the -S option\n"
    usage(1)
  end

  nthreads = if cl.option_present('thr')
               cl.value('thr')
             else
               2
             end
  
  # We prefer /node/scratch to /tmp as it is usually larger
  tmpdir = '/tmp'
  if FileTest.directory?('/node/scratch')
    tmpdir = '/node/scratch'
  end
  
  if cl.option_present('tmpdir')
    tmpdir = cl.value('tmpdir')
  
    Dir.mkdir(tmpdir) unless (FileTest.directory?(tmpdir))
  end
  
  PSupport::looks_like_multiple_files(ARGV) && exit

  input_file = ARGV.pop

  offset = PSupport::get_offsets(input_file, nthreads, verbose)
  
  sfile = cl.value('S')
  
  rejoin = if cl.option_present('nj')
             false
           else
             true
           end
  
  stem = "#{tmpdir}/#{File.basename(sfile)}"
  # If we are rejoining, a file that will be removed.
  stem << "#{Process.pid}" if rejoin
  
  common_args = PSupport::quote_special_characters(ARGV)

  common_args << ' -n'
  
  log_stem = if cl.option_present('log')
              cl.value('log')
            else
              false
            end

  keeplog = cl.option_present('keeplog')
  
  summary_files = (0...nthreads).to_a.map { |i| "#{stem}#{i}.textproto" }
  log_files = if log_stem
                (0...nthreads).to_a.map { |i| "#{log_stem}#{i}.log" }
              else
                []
              end

  pids = Array.new
  
  (1...nthreads).each do |i|
    cmd = "#{get_substituents} #{common_args} -S #{summary_files[i]}"
    cmd << ' -v' if keeplog
  
    cmd << " -i seek=#{offset[i]}"
    cmd << " -i stop=#{offset[i+1]}" if offset[i+1]
      
    cmd << " #{input_file}"
  
    cmd << " 2> #{log_files[i]}" if log_stem
  
    $stderr << "Executing '#{cmd}'\n" if verbose
  
    p = fork {
      exec(cmd)
    }
  
    pids.push(p)
  end
  
  # Now do thread zero here
  
  cmd = "#{get_substituents} #{common_args} -S #{summary_files[0]}"
  cmd << ' -v' if keeplog
  
  cmd << " -i stop=#{offset[1]} #{input_file}"

  cmd << " 2>#{log_files[0]}" if log_stem
  
  $stderr << "Thread zero '#{cmd}'\n" if verbose
  system(cmd)

  pids.each do |p|
    Process.wait(p)
  end
  
  # Append suffixes to chunk files.

  exit 1 unless PSupport::all_files_exist(summary_files)
  
  PSupport::remove_files(log_files) unless keeplog

  exit 0 unless rejoin

  cmd = "dicer_fragments_collate.sh " << summary_files.join(' ') << " > #{sfile}"
  
  $stderr << "Executing #{cmd}\n" if verbose
  system(cmd)
  
  PSupport::remove_files(summary_files)
end

get_substituents_parallel
