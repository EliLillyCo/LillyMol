#!/usr/bin/env ruby

require_relative "lib/iwcmdline_v2.rb"
require_relative "lib/parallel_support.rb"

$expert = false

def usage(rc)
  $stderr << "Multi-threaded version of gfp_make\n"
  $stderr << "gfp_make_parallel -thr 16 ... file.smi > file.gfp\n"
  $stderr << "Note that only one input file can be processed\n";
  $stderr << " -thr <nthreads>  number of threads to use\n"
  $stderr << " -gfp_make <exe>  gfp_make executable to use (default gfp_make.sh)\n"
  $stderr << " -tmpdir <dir>    directory for temporary files\n" if ($expert)
  $stderr << " -log <stem>      redirect stderr to individual log files - which are removed\n"
  $stderr << "                  you probably want this to avoid screens full of warning messages\n"
  $stderr << "                  but make sure your reaction is well behaved and that warnings can actually be ignored\n";
  $stderr << " -keeplog         do NOT remove the logfiles\n"
  $stderr << " -v               verbose output\n"
  $stderr << "All other options are passed directly to gfp_make\n"

  exit(rc)
end

def gfp_make_parallel

  cl = IWCmdlineV2.new("-v-thr=ipos-tmpdir=s-gfp_make=xfile-log=s-keeplog")
  
  verbose = cl.option_present('v')
  
  gfp_make = if cl.option_present('gfp_make') 
            cl.value('gfp_make')
          else
            'gfp_make.sh'
          end
  
  if ARGV.empty?
    $stderr << "Insufficient arguments\n"
    usage(2)
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
  
  stem = File.join(tmpdir, "gfp_make_p#{Process.pid}")
  
  common_args = PSupport::quote_special_characters(ARGV)
  
  log_stem = if cl.option_present('log')
              cl.value('log')
            else
              false
            end

  keeplog = cl.option_present('keeplog')
  
  output_files = (0...nthreads).to_a.map { |i| "#{stem}#{i}.gfp" }
  log_files = if log_stem
                (0...nthreads).to_a.map { |i| "#{log_stem}#{i}.log" }
              else
                []
              end

  pids = Array.new
  
  (1...nthreads).each do |i|
    cmd = "#{gfp_make} #{common_args}"
  
    cmd << " -i seek=#{offset[i]}"
    cmd << " -i stop=#{offset[i+1]}" if offset[i+1]
      
    cmd << " #{input_file} > #{output_files[i]}"
  
    cmd << " 2> #{log_files[i]}" if log_stem
  
    $stderr << "Executing '#{cmd}'\n" if verbose
  
    p = fork {
      exec(cmd)
    }
  
    pids.push(p)
  end
  
  # Now do thread zero here
  
  cmd = "#{gfp_make} #{common_args}"
  
  cmd << " -i stop=#{offset[1]} #{input_file}"

  cmd << " > #{output_files[0]}"

  cmd << " 2>#{log_files[0]}" if log_stem
  
  $stderr << "Thread zero '#{cmd}'\n" if verbose
  system(cmd)

  pids.each do |p|
    Process.wait(p)
  end
  
  exit 1 unless PSupport::all_files_exist(output_files)
  
  PSupport::remove_files(log_files) unless keeplog

  cmd = "cat " << output_files.join(' ')
  
  $stderr << "Executing #{cmd}\n" if verbose
  system(cmd)
  
  PSupport::remove_files(output_files)
end

gfp_make_parallel
