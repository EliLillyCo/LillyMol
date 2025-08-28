#!/usr/bin/env ruby

require_relative "lib/iwcmdline_v2.rb"
require_relative "lib/parallel_support.rb"

$expert = false

def usage(rc)
  $stderr << "Multi-threaded version of molecule_filter\n"
  $stderr << "fileconv_parallel -thr 16 ... file1\n"
  $stderr << "Note that only one input file can be processed\n";
  $stderr << " -thr <nthreads>  number of threads to use\n"
  $stderr << " -molecule_filter <exe>  molecule_filter executable to use (default molecule_filter.sh)\n"
  $stderr << " -tmpdir <dir>    directory for temporary files\n" if ($expert)
  $stderr << " -log <stem>      redirect stderr to individual log files - which are removed\n"
  $stderr << "                  you probably want this to avoid screens full of warning messages\n"
  $stderr << "                  but make sure your reaction is well behaved and that warnings can actually be ignored\n";
  $stderr << " -keeplog         do NOT remove the logfiles, passes -v to the underlying executable\n"
  $stderr << " -v               verbose output\n"
  $stderr << "All other options are passed directly to molecule_filter\n"

  exit(rc)
end

def molecule_filter_parallel

  cl = IWCmdlineV2.new("-v-thr=ipos-tmpdir=s-molecule_filter=xfile-log=s-keeplog-F=sfile")
  
  verbose = cl.option_present('v')
  
  molecule_filter = if cl.option_present('molecule_filter') 
            cl.value('molecule_filter')
          else
            'molecule_filter.sh'
          end
  
  if ARGV.empty?
    $stderr << "Insufficient arguments\n"
    usage(2)
  end

  unless cl.option_present('F')
    $stderr << "Must specify config file via the -F option\n"
    usage(1)
  end

  config_file = cl.value('F')
  
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
  
  stem = "#{tmpdir}/mfparallel#{Process.pid}"
  
  common_args = PSupport::quote_special_characters(ARGV)
  
  log_stem = if cl.option_present('log')
              cl.value('log')
            else
              false
            end

  keeplog = cl.option_present('keeplog')
  
  smiles_files = (0...nthreads).to_a.map { |i| "#{stem}#{i}" }
  log_files = if log_stem
                (0...nthreads).to_a.map { |i| "#{log_stem}#{i}.log" }
              else
                []
              end

  # The default command.
  base_cmd = "#{molecule_filter} -F #{config_file} #{common_args}"
  base_cmd << ' -v' if keeplog
  pids = Array.new
  
  (1...nthreads).each do |i|
    cmd = "#{base_cmd}"
  
    cmd << " -i seek=#{offset[i]}"
    cmd << " -i stop=#{offset[i+1]}" if offset[i+1]
      
    cmd << " #{input_file}"

    cmd << " > #{smiles_files[i]}"
  
    cmd << " 2> #{log_files[i]}" if log_stem
  
    $stderr << "Executing '#{cmd}'\n" if verbose
  
    p = fork {
      exec(cmd)
    }
  
    pids.push(p)
  end
  
  # Now do thread zero here
  
  cmd = "#{base_cmd}"
  
  cmd << " -i stop=#{offset[1]} #{input_file}"

  cmd << " 2>#{log_files[0]}" if log_stem
  
  $stderr << "Thread zero '#{cmd}'\n" if verbose
  system(cmd)

  pids.each do |p|
    Process.wait(p)
  end
  
  # We did not create the first temporary file.
  smiles_files.shift
  exit 1 unless PSupport::all_files_exist(smiles_files)
  
  PSupport::remove_files(log_files) unless keeplog

  cmd = "/bin/cat " << smiles_files.join(' ')
  
  $stderr << "Executing #{cmd}\n" if verbose
  system(cmd)
  
  PSupport::remove_files(smiles_files)
end

molecule_filter_parallel
