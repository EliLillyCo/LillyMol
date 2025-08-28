#!/usr/bin/env ruby

require_relative "lib/iwcmdline_v2.rb"
require_relative "lib/parallel_support.rb"

$expert = false

def usage(rc)
  $stderr << "Multi-threaded version of in_database_bdb\n"
  $stderr << "Note that only the -F and -U options are aggregated across parallel lookups\n"
  $stderr << "in_database_parallel -thr 16 -d out.bdb ... file1\n"
  $stderr << "Note that only one input file can be processed\n";
  $stderr << " -d <dbname>      the database in which to perform lookups - mandatory\n"
  $stderr << " -F <fname>       write molecules     found in the database to <fname>\n"
  $stderr << " -U <fname>       write molecules not found in the database to <fname>\n"
  $stderr << "                  must specify one or both of -F and -U\n"
  $stderr << " -thr <nthreads>  number of threads to use\n"
  $stderr << " -nj              do NOT join the output files, leave in split form\n"
  $stderr << " -in_database <exe>  in_database_bdb executable to use (default in_database_bdb.sh)\n"
  $stderr << " -tmpdir <dir>    directory for temporary files\n" if ($expert)
  $stderr << " -log <stem>      redirect stderr to individual log files - which are removed\n"
  $stderr << "                  you probably want this to avoid screens full of warning messages\n"
  $stderr << "                  but make sure your reaction is well behaved and that warnings can actually be ignored\n";
  $stderr << " -keeplog         do NOT remove the logfiles, passes -v to the underlying executable\n"
  $stderr << " -v               verbose output\n"
  $stderr << "All other options are passed directly to in_database_bdb\n"

  exit(rc)
end

def append_f_and_u(cmd, thr, f_files, u_files)
  if f_files
    cmd << " -F #{f_files[thr]}"
  end
  if u_files
    cmd << " -U #{u_files[thr]}"
  end
end

def in_database_parallel

  cl = IWCmdlineV2.new("-v-thr=ipos-tmpdir=s-d=s-in_database=xfile-nj-log=s-keeplog-F=s-U=s")
  
  verbose = cl.option_present('v')
  
  in_database_bdb = if cl.option_present('in_database') 
            cl.value('in_database')
          else
            'in_database_bdb.sh'
          end
  
  if ARGV.empty?
    $stderr << "Insufficient arguments\n"
    usage(2)
  end
  
  unless cl.option_present('d')
    $stderr << "Must specify database file name via the -d option\n"
    usage(1)
  end

  fstem = if cl.option_present('F')
            cl.value('F')
          else
            ""
          end
  ustem = if cl.option_present('U')
            cl.value('U')
          else
            ""
          end

  if fstem.empty? && ustem.empty?
    $stderr << "Must specify one of both -F and -U options\n"
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
  
  dbname = cl.value('d')
  unless File.size?(dbname)
    $stderr << "Missing or empty database #{dbname}\n"
    return 1
  end

  rejoin = if cl.option_present('nj')
             false
           else
             true
           end
  
  ftmp = "#{tmpdir}/#{File.basename(fstem)}"
  utmp = "#{tmpdir}/#{File.basename(fstem)}"
  # If we are rejoining, a file that will be removed.
  if rejoin
    ftmp << "#{Process.pid}" unless fstem.empty?
    ftmp << "#{Process.pid}" unless ustem.empty?
  end
  
  common_args = PSupport::quote_special_characters(ARGV)

  f_files = (0...nthreads).to_a.map { |i| "#{ftmp}#{i}" } if fstem.length > 0
  u_files = (0...nthreads).to_a.map { |i| "#{utmp}#{i}" } if ustem.length > 0
  
  log_stem = if cl.option_present('log')
              cl.value('log')
            else
              false
            end

  keeplog = cl.option_present('keeplog')
  
  f_files = (0...nthreads).to_a.map { |i| "#{ftmp}#{i}" } if fstem.length > 0
  u_files = (0...nthreads).to_a.map { |i| "#{utmp}#{i}" } if ustem.length > 0

  log_files = if log_stem
                (0...nthreads).to_a.map { |i| "#{log_stem}#{i}.log" }
              else
                []
              end

  pids = Array.new
  
  (1...nthreads).each do |i|
    cmd = "#{in_database_bdb} #{common_args} -d #{dbname}"
    append_f_and_u(cmd, i, f_files, u_files)
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
  
  cmd = "#{in_database_bdb} #{common_args} -d #{dbname}"
  append_f_and_u(cmd, 0, f_files, u_files)
  cmd << ' -v' if keeplog
  
  cmd << " -i stop=#{offset[1]} #{input_file}"

  cmd << " 2>#{log_files[0]}" if log_stem
  
  $stderr << "Thread zero '#{cmd}'\n" if verbose
  system(cmd)

  pids.each do |p|
    Process.wait(p)
  end
  
  f_files = f_files.map { |s| "#{s}.smi" } if f_files
  u_files = u_files.map { |s| "#{s}.smi" } if u_files

  exit 1 unless PSupport::all_files_exist(f_files) if f_files
  exit 1 unless PSupport::all_files_exist(u_files) if u_files
  
  PSupport::remove_files(log_files) unless keeplog

  exit 0 unless rejoin

  if f_files
    cmd = "cat " << f_files.join(' ') << " > #{fstem}.smi"
  
    $stderr << "Executing #{cmd}\n" if verbose
    system(cmd)
  end

  if u_files
    cmd = "cat " << u_files.join(' ') << " > #{ustem}.smi"
  
    $stderr << "Executing #{cmd}\n" if verbose
    system(cmd)
  end

  PSupport::remove_files(f_files) if f_files
  PSupport::remove_files(u_files) if u_files
end

in_database_parallel
