#!/usr/bin/env ruby

require_relative "lib/iwcmdline_v2.rb"
require_relative "lib/parallel_support.rb"

$expert = false

def usage(rc)
  $stderr << "Multi-threaded version of unique_molecules\n"
  $stderr << "unique_molecules_parallel -thr 16 -S out ... file1\n"
  $stderr << "Note that only one input file can be processed\n";
  $stderr << " -S <stem>        output stem - mandatory\n"
  $stderr << " -thr <nthreads>  number of threads to use\n"
  $stderr << " -nj              do NOT join the output files, leave in split form\n"
  $stderr << " -unique_molecules <exe>  unique_molecules executable to use (default unique_molecules.sh)\n"
  $stderr << " -tmpdir <dir>    directory for temporary files\n" if ($expert)
  $stderr << " -log <stem>      redirect stderr to individual log files - which are removed\n"
  $stderr << "                  you probably want this to avoid screens full of warning messages\n"
  $stderr << "                  but make sure your reaction is well behaved and that warnings can actually be ignored\n";
  $stderr << " -keeplog         do NOT remove the logfiles\n"
  $stderr << " -v               verbose output\n"
  $stderr << "All other options are passed directly to unique_molecules\n"

  exit(rc)
end

def unique_molecules_parallel

  cl = IWCmdlineV2.new("-v-thr=ipos-tmpdir=s-S=s-unique_molecules=xfile-nj-log=s-keeplog")

  verbose = cl.option_present('v')

  unique_molecules = if cl.option_present('unique_molecules') 
            cl.value('unique_molecules')
          else
            'unique_molecules.sh'
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

  stem = File.join(tmpdir, "#{File.basename(sfile)}")

  dstem = File.join(tmpdir, "D#{File.basename(sfile)}")

  # First task is to sort the file
  cmd = "msort_parallel.sh -p -D #{dstem} -M nc=#{nthreads} -k formula #{input_file}"
  $stderr << "Executing #{cmd}\n" if verbose
  system(cmd)

  dfiles = []
  (0...nthreads).each do |i|
    fname = "#{dstem}#{i}.smi"
    unless File.size?(fname)
      $stderr << "msort did not create #{fname}\n"
      return 1;
    end

    dfiles << fname
  end

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

  pids = Array.new

  (0...nthreads).each do |i|
    cmd = "#{unique_molecules} #{common_args} -S #{smiles_files[i]} #{dfiles[i]}"
    cmd << ' -v' if keeplog

    cmd << " 2> #{log_files[i]}" if log_stem

    $stderr << "Executing '#{cmd}'\n" if verbose

    p = fork {
      exec(cmd)
    }

    pids.push(p)
  end

  pids.each do |p|
    Process.wait(p)
  end

  # Append suffixes to chunk files.
  smiles_files = smiles_files.map { |fname| "#{fname}.smi" }

  exit 1 unless PSupport::all_files_exist(smiles_files)

  PSupport::remove_files(log_files) unless keeplog

  exit 0 unless rejoin

  cmd = "cat " << smiles_files.join(' ') << " > #{sfile}.smi"

  $stderr << "Executing #{cmd}\n" if verbose
  system(cmd)

  PSupport::remove_files(smiles_files)
end

unique_molecules_parallel
