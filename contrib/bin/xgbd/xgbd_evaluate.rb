#!/usr/bin/env ruby

# Evaluate an xgboost model

require 'set'
require 'google/protobuf'

c3tk_home = ENV['C3TK_HOME']
raise 'C3TK_HOME not defined' unless c3tk_home

require "#{c3tk_home}/bin/ruby/lib/iwcmdline"
require "#{c3tk_home}/bin/py/pytk/xgbd/xgboost_model_pb"

def usage
msg = <<-END
Scores an xgboost descriptor model built with xgbd_make.
Takes either a smiles file or a descriptor file,
If a smiles file is given, the descriptors required for the model will be computed by make_descriptors
  -mdir <dir>           model directory created by xgbd_make
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

  return XgboostModel::XGBoostModel.decode(File.read(fname))
end

def xgbd_evaluate_smiles(fname, mdir, proto, cl)
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

  descriptors.each do |d|
    cmd << " -#{d}"
  end
  cmd << " #{fname}"
  cmd << "| xgboost_model_evaluate.sh -mdir #{mdir} -"

  $stderr << "Executing #{cmd}\n" if cl.option_present('v')

  system(cmd)
end

def xgbd_evaluate_descriptors(fname, mdir, proto, cl)
  cmd = "xgboost_model_evaluate.sh -mdir #{mdir} #{fname}"

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
    $stderr << "Must specify xgboost model directory via the -mdir option\n"
    usage
  end

  if ARGV.empty?
    $stderr << "Insufficient arguments\n"
    usage
  end

  mdir = cl.value('mdir')

  proto = get_model_metadata(mdir)

  if cl.option_present('smi') || ARGV[0].match(/\.smi$/)
    xgbd_evaluate_smiles(ARGV[0], mdir, proto, cl)
  else
    xgbd_evaluate_descriptors(ARGV[0], mdir, proto, cl)
  end

end

main
