#!/usr/bin/env ruby

require_relative "lib/iwcmdline_v2.rb"
require_relative "lib/parallel_support.rb"

$expert = false

def usage(rc)
  $stderr << "Multi-threaded version of trxn\n"
  $stderr << "trxn_parallel -thr 16 -P rxn.textproto -S product ... scaffold.smi sidechain.smi\n"
  $stderr << " -R <rxn>         reaction file\n"
  $stderr << " -P <rxn>         textproto reaction file\n"
  $stderr << " -S <stem>        output file name - mandatory\n"
  $stderr << " -thr <nthreads>  number of threads to use\n"
  $stderr << " -nj              do NOT join the output files, leave in split form\n"
  $stderr << " -trxn <exe>      trxn executable to use (default trxn.sh)\n"
  $stderr << " -tmpdir <dir>    directory for temporary files\n" if ($expert)
  $stderr << " -log <stem>      redirect stderr to individual log files - which are removed\n"
  $stderr << "                  you probably want this to avoid screens full of warning messages\n"
  $stderr << "                  but make sure your reaction is well behaved and that warnings can actually be ignored\n";
  $stderr << " -keeplog         do NOT remove the logfiles, passes -v to the underlying executable\n"
  $stderr << " -v               verbose output\n"
  $stderr << "All other options are passed directly to trxn\n"

  exit(rc)
end

# Scan backwards through argv looking for the scaffold file.
def get_input_file(argv)
  # The last token on the command line must be a file of structures,
  # either a sidechain file or the scaffold
  # $stderr << "Examining #{argv}\n"
  unless File.size?(argv.last)
    $stderr << "Missing or empty file #{argv.last}\n"
    exit 1
  end

  # Something like
  # trxn_parallel -P rxn.textproto -thr 8 file.smi
  # is quite OK.
  return argv[0] if argv.size == 1

  prev_file = argv.last
  (argv.length - 2).downto(0) do |i|
    fname = argv[i]
    if File.size?(fname)
      prev_file = fname
    else
      return prev_file
    end
  end

  return argv[0]
end

def trxn_parallel

  cl = IWCmdlineV2.new("-v-thr=ipos-tmpdir=s-S=s-trxn=xfile-R=sfile-P=sfile-nj-log=s-keeplog")

  verbose = cl.option_present('v')

  trxn = if cl.option_present('trxn')
            cl.value('trxn')
          else
            'trxn.sh'
          end

  if ARGV.empty?
    $stderr << "Insufficient arguments\n"
    usage(2)
  end

  unless cl.option_present('S')
    $stderr << "Must specify output file name via the -S option\n"
    usage(1)
  end

  if cl.option_present('P')
    rxnfile = "-P " << cl.value('P')
  elsif cl.option_present('R')
   rxnfile = "-R " << cl.value('R')
  else
    $stderr << "must specify reactin with either the -R or -P option\n"
    usage
  end

  nthreads = if cl.option_present('thr')
               cl.value('thr')
             else
               2
             end

  # Scan backwards from the end of ARGV to identify the scaffold file

  input_file = get_input_file(ARGV)

  unless File.size?(input_file)
    $stderr << "Missing or empty input file '#{input_file}'\n"
    exit 2
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

  common_args = "#{rxnfile} " << PSupport::quote_special_characters(ARGV)

  keeplog = cl.option_present('keeplog')

  log_stem = if cl.option_present('log')
              cl.value('log')
            else
              false
            end

  smiles_files = (0...nthreads).to_a.map { |i| "#{stem}#{i}" }

  log_files = if log_stem
                (0...nthreads).to_a.map { |i| "#{log_stem}#{i}.log" }
              else
                []
              end

  files_to_delete = []

  pids = Array.new

  (1...nthreads).each do |i|
    cmd = "#{trxn} #{common_args} -S #{smiles_files[i]}"
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

  cmd = "#{trxn} #{common_args} -S #{smiles_files[0]}"
  cmd << ' -v' if keeplog

  cmd << " -i stop=#{offset[1]} #{input_file}"

  cmd << " 2>#{log_files[0]}" if log_stem

  $stderr << "Thread zero '#{cmd}'\n" if verbose
  system(cmd)

  pids.each do |p|
    Process.wait(p)
  end

  # Append suffixes to chunk files.
  smiles_files = smiles_files.map { |fname| "#{fname}.smi" }

  exit 1 unless PSupport::all_files_exist(smiles_files)

  exit 0 unless rejoin

  cmd = "cat " << smiles_files.join(' ') << " > #{stem}.smi"
  $stderr << "Executing #{cmd}\n" if verbose
  system(cmd)

  PSupport::remove_files(smiles_files)
end

trxn_parallel
