#!/usr/bin/env ruby

# Evaluate an RF model from either smiles or a descriptor file.

require 'set'
require 'google/protobuf'

c3tk_home = ENV['C3TK_HOME']
raise 'C3TK_HOME not defined' unless c3tk_home

require "#{c3tk_home}/bin/ruby/lib/iwcmdline"
require "#{c3tk_home}/bin/py/pytk/xgbd/random_forest_model_pb"

def usage
msg = <<-END
Scores a random_forest descriptor model built with rf_make.
Takes either a smiles file or a descriptor file,
If a smiles file is given, the descriptors required for the model will be computed by make_descriptors
  -mdir <dir>           model directory created by rf_make
  -smi                  input is a smiles file - model descriptors will be computed.
  -pipe                 use descriptor_pipeline for descriptor compuation,
                        only some descriptors are supported, but maybe faster.
  -j <threads>          send the -j option to make_descriptors.sh
  -v                    verbose output
END
  $stderr << msg

  exit(0)
end

# Return the proto mdir/model_metadata.dat
def get_model_metadata(mdir)
  fname = File.join(mdir, 'model_metadata.dat')

  unless File.size?(fname)
    $stderr << "Empty or missing model metadata #{fname}\n"
    return 1
  end

  return RandomForestModel::RandomForestModel.decode(File.read(fname))
end

def rf_evaluate_smiles(fname, mdir, proto, cl)
  descriptors = Set.new()
  proto.name_to_col.each do |name, col|
    descriptors.add(name.gsub(/_.*/, ""))
  end

  if cl.option_present('pipe')
    cmd = 'descriptor_pipeline.sh'
  else
    cmd = 'make_descriptors.sh'
  end

  if cl.option_present('j')
    j = cl.value('j')
    cmd << " -j #{j}"
  end

  tmpfile = File.join(ENV['TMPDIR'], "rf_evaluate_smiles_#{Process.uid}.#{Process.pid}.dat")

  descriptors.each do |d|
    cmd << " -#{d}"
  end
  cmd << " #{fname} > #{tmpfile}"
  $stderr << "Executing #{cmd}\n" if cl.option_present('v')

  system(cmd)
  unless File.size?(tmpfile)
    $stderr << "#{cmd} failed\n"
    return
  end

  rf_evaluate_descriptors(tmpfile, mdir, proto, cl)

  File.unlink(tmpfile)
end

def rf_evaluate_descriptors(fname, mdir, proto, cl)
  cmd = "rf_evaluate.sh -mdir #{mdir} #{fname}"

  $stderr << "Executing #{cmd}\n" if cl.option_present('v')
  system(cmd)
end

def main
  cl = IWCmdline.new("-v-mdir=dir-smi-pipe-j=ipos")

  if cl.unrecognised_options_encountered
    $stderr << "unrecognised_options_encountered\n"
    usage
  end

  unless cl.option_present('mdir')
    $stderr << "Must specify random forest model directory via the -mdir option\n"
    usage
  end

  if ARGV.empty?
    $stderr << "Insufficient arguments\n"
    usage
  end

  mdir = cl.value('mdir')

  proto = get_model_metadata(mdir)

  if cl.option_present('smi') || ARGV[0].match(/\.smi$/)
    rf_evaluate_smiles(ARGV[0], mdir, proto, cl)
  else
    rf_evaluate_descriptors(ARGV[0], mdir, proto, cl)
  end

end

main
