# Evaluate an svmfp model built with svmfp_make

# frozen_string_literal: true

require 'google/protobuf'

require_relative 'lib/iwcmdline'
require_relative 'lib/gfp_model_pb'

def usage(retcod)
  $stderr << "Evaluate an svmfp model build with svmfp_make\n"
  $stderr << " -m <dir>     model directory\n"
  $stderr << " -v           verbose output\n"
  exit(retcod)
end

cmdline = IWCmdline.new('-v-mdir=dir-gfp_make=xfile-eval=xfile-eval_opts=close')
if cmdline.unrecognised_options_encountered
  $stderr << "unrecognised_options_encountered\n"
  usage(1)
end

verbose = cmdline.option_present('v')

unless cmdline.option_present('mdir')
  $stderr << "Must specify model directory via the -mdir option\n"
  usage(1)
end

mdir = cmdline.value('mdir')

# A model that will be scored. Attributes are populated by reading
# the proto file in the directory. Note that the proto file is encoded.
class Model
  attr_reader :model_data_fname

  def initialize(mdir)
    @mdir = mdir
    @model_data_fname = "#{mdir}/model.dat"
    raise "Missing #{@model_data_fname}" unless File.size?(@model_data_fname)

    @data = GfpModel::SvmfpModel.decode(File.read(@model_data_fname))
  end

  def gfp
    @data.metadata.fingerprints
  end
end

# When give a model directory, we decide if that is a single
# model directory, or a directory containing other model directories.

def identify_models(dir)
  model_proto = "#{dir}/model.dat"
  return [dir] if File.size?(model_proto)

  result = []
  Dir.open(dir).each do |fname|
    next if fname.start_with?('.')
    next unless File.directory?(fname)

    model_proto = "#{dir}/model.dat"
    result << dir if File.size?(model_proto)
  end
  result
end

model_dirs = identify_models(mdir)
if model_dirs.empty?
  $stderr << "Invalid model directory #{mdir}\n"
  exit(1)
end

models = []
model_dirs.each do |dir|
  models << Model.new(dir)
end
# $stderr << models << "\n"

if ARGV.empty?
  $stderr << "Insufficient arguments\n"
  usage(1)
end

gfp_make = if cmdline.option_present('gfp_make')
             cmdline.value('gfp_make')
           else
             'gfp_make.sh'
           end

evaluate = if cmdline.option_present('evaluate')
             cmdline.value('evaluate')
           else
             'gfp_svmfp_score_v2'
           end

if cmdline.option_present('eval_opts')
  opts = cmdline.value('eval_opts')
  evaluate = "#{evaluate} #{opts}"
elsif verbose
  evaluate = "#{evaluate} -v"
end

smiles_files = ARGV.join(' ')

# For now, multiple models are not combined, the output will contain
# a column for each model.

mnames = []
models.each do |model|
  mnames << '-M ' << model.model_data_fname
end
mnames = mnames.join(' ')
# $stderr << mnames << "\n"

models.each do |model|
  cmd = "#{gfp_make}  #{model.gfp} #{smiles_files} | #{evaluate} -M #{model.model_data_fname} -"
  $stderr << "Executing #{cmd}\n" if verbose
  system(cmd)
end
