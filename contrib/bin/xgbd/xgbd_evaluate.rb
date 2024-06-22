# Front end for xgbd_evaluate.py

require_relative("../lib/iwcmdline.rb")

def usage(rc)
  exit(rc)
end

def main
  cl = IWCmdline.new("-v-mdir=dir-tmpdir=dir")

  verbose = cl.option_present("v")

  if cl.unrecognised_options_encountered()
    $stderr << "unrecongised_options_encountered\n"
    usage(1)
  end

  unless cl.option_present("mdir")
    $stderr << "Must specify model directory via the -mdir option\n"
    usage(1)
  end

  mdir = cl.value("mdir")
  $stderr << "mdir #{mdir}\n"
  traindat = File.join(mdir, "train.dat")
  unless File.size?(traindat)
    $stderr << "Missing or empty training data file #{traindat}\n"
    return 1
  end

  tmpdir = if cl.option_present("tmpdir")
             cl.value("tmpdir")
           else
             "/tmp"
           end

  if ARGV.size == 0
    $stderr << "Insufficient arguments\n"
    usage(1)
  end

  tmpfile = File.join(tmpdir, "xgbd_eval#{Process.pid}.dat")
  cmd = "iwcut -D #{traindat} -v #{ARGV[0]} > #{tmpfile}"
  $stderr << "Executing #{cmd}\n" if verbose
  rc = system(cmd)
  if rc && File.size?(tmpfile)
  else
    $stderr << "#{cmd} failed, did not generate #{tmpfile}\n"
  end

  python_script = __FILE__.gsub("\.rb", ".py")
  cmd = "python #{python_script} -mdir #{mdir} #{tmpfile}"
  $stderr << "Executing #{cmd}\n" if verbose
  system(cmd)


  File.unlink(tmpfile)
end

main()
