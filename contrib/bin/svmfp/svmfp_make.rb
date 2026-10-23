# Build svmfp models
# Subsequently altered to also build lightgbm models.

# frozen_string_literal: true

require 'date'
require 'fileutils'
require 'tempfile'
require 'google/protobuf'

require_relative 'lib/iwcmdline'
require_relative 'lib/gfp_model_pb'

def usage(retcod)
  $stderr << "Builds an svmfp model from smiles and activity\n"
  $stderr << " -mdir <dir>   model directory to create\n"
  $stderr << " -A <fname>    file containing activity data\n"
  $stderr << " -C            classification model\n"
  $stderr << " -gfp ... -gfp fingerprint specification (to gfp_make)\n"
  $stderr << " -p <support>  support level for bit inclusion\n"
  $stderr << " -m <ram>          MB of RAM used for the cache. Default 500\n"
  $stderr << " -w <width>    svm_learn tube width (default 0.1)\n"
  $stderr << " -svml ... -svml  passed directly to svm_learn\n"
  $stderr << " -flatten      flatten sparse fingerprint counts to 1\n"
  $stderr << " -lightgbm ... -lightgbm build a lightgbm model\n"
  $stderr << " -catboost ... -catboost build a catboost model\n"
  $stderr << " -rescore       rescore the training set to establish a linear correction\n"
  $stderr << " -v            verbose output\n"

  exit retcod
end

def execute_cmd(cmd, verbose, files_needed = [])
  $stderr << "Executing '#{cmd}'\n" if verbose
  rc = system(cmd)
  raise "Non zero rc #{rc} from #{cmd}" unless rc

  files_needed.each do |fname|
    next if fname == nil

    next if File.size?(fname)

    raise "#{cmd} did not produce #{fname}"
  end
end

def get_threshold_b(model_file)
  threshold_b_rx = Regexp.new('(\S+) # threshold b, each following line')

  File.foreach(model_file).each do |line|
    m = threshold_b_rx.match(line)
    return m[1].to_f if m
  end
  raise "No threshold_b in #{model_file}"
end

# A model file has been created in `mdir`. return the name.
# Note that this will silently fail of there are multiple model
# files in the directory and the one we want is not the first found.
def get_xgboost_file(mdir)
  maybe = Dir.children(mdir).filter { |fname| /^\d+\.model$/.match(fname) }

  if maybe.empty?
    $stderr << "No model file in #{mdir}\n"
    return nil
  end

  # sort by time modified and take the most recent.

  maybe.sort_by! { |t| File.mtime("#{mdir}/#{t}") }
  return maybe.last
end

# `activity_file` has been identified as containing classification data.
# generate a class label translation file in `mdir`.
def perform_class_label_translation(activity_file, mdir, train_activity, verbose)
  cmd = "class_label_translation -C #{mdir}/class_label_translation.txt " \
  "-cbin #{mdir}/class_label_translation.dat #{activity_file} > #{train_activity}"
  execute_cmd(cmd, verbose, [train_activity, "#{mdir}/class_label_translation.dat"])
end

# Only difference is that an extra option is needed for class_label_translation.
def perform_class_label_translation_lightgbm(activity_file, mdir, train_activity, verbose)
  cmd = "class_label_translation -lightgbm -C #{mdir}/class_label_translation.txt " \
  "-cbin #{mdir}/class_label_translation.dat #{activity_file} > #{train_activity}"
  execute_cmd(cmd, verbose, [train_activity, "#{mdir}/class_label_translation.dat"])
end

def perform_response_scaling(activity_file, mdir, train_smi, train_activity, verbose)
  cmd = "feature_scaling -i auto -bin -C #{mdir}/response_scaling -v -subset #{train_smi} -scol 2 #{activity_file} > #{train_activity}"
  execute_cmd(cmd, verbose, [train_activity])
end

def get_response_name(activity_file)
  input = File.open(activity_file, 'r')
  header = input.readline
  input.close
  header.split[1]
end

cmdline = IWCmdline.new('-v-mdir=s-A=sfile-activity=sfile-C-gfp=close-svml=close-p=ipos-w=fraction-flatten-gfp_make=xfile' \
                        '-svm_learn=xfile-gfp_to_svm_lite=xfile-lightgbm=close-lightgbm_config=sfile' \
                        '-catboost=close' \
                        '-xgboost=close-xgboost_config=sfile-xgb_test=sfile' \
                        '-m=ipos' \
                        '-rescore')
if cmdline.unrecognised_options_encountered
  $stderr << "unrecognised_options_encountered\n"
  usage(1)
end

verbose = cmdline.option_present('v')

if cmdline.option_present('A')
  activity_file = cmdline.value('A')
elsif cmdline.option_present('activity')
  activity_file = cmdline.value('activity')
else
  $stderr << "Must specify activity file via the -A option\n"
  usage(1)
end

unless cmdline.option_present('mdir')
  $stderr << "Must specify model directory via the -mdir option\n"
  usage(1)
end

mdir = cmdline.value('mdir')

begin
  FileUtils.mkdir_p(mdir)
rescue => e # rubocop:disable Style/RescueStandardError
  $stderr << "Did not create model directory #{mdir} #{e.message}'\n"
  exit(1)
end

gfp_make = if cmdline.option_present('gfp_make')
             cmdline.value('gfp_make')
           else
             'gfp_make.sh'
           end

gfp_to_svm_lite = if cmdline.option_present('gfp_to_svm_lite')
                    cmdline.value('gfp_to_svm_lite')
                  else
                    'gfp_to_svm_lite'
                  end

svm_learn = if cmdline.option_present('svm_learn')
              cmdline.value('svm_learn')
            else
              'svm_learn'
            end

svm_learn_options = if cmdline.option_present('svml')
                      cmdline.value('svml')
                    else
                      '-t 4 -m 500'
                    end

svm_learn_options = "#{svm_learn_options} -w #{cmdline.value('w')}" if cmdline.option_present('w')
svm_learn_options = "#{svm_learn_options} -m #{cmdline.value('m')}" if cmdline.option_present('m')

# the alternate model forms will be nil if not specified.
lightgbm = cmdline.value('lightgbm')
default_lightgbm_config = cmdline.value('lightgbm_config')
catboost = cmdline.value('catboost')
xgboost = cmdline.value('xgboost')
xgboost_config = cmdline.value('xgboost_config')

if lightgbm
  if ! default_lightgbm_config
    $stderr << "When building a lightgbm model, must specify -lightgbm_config\n"
    usage(1)
  end
  lightgbm = "lightgbm config=#{default_lightgbm_config} #{lightgbm} force_row_wise=true"
  FileUtils.cp(default_lightgbm_config, mdir)
end

if xgboost
  if ! xgboost_config
    $stderr << "When building a xgboost model, must specify -xgboost_config\n"
    usage(1)
  end
  xgboost = "xgboost #{xgboost_config} model_dir=#{mdir}"
  FileUtils.cp(xgboost_config, mdir)
end

catboost = "catboost fit #{catboost} --train-dir #{mdir} --fstr-file fstr.dat " \
           "--use-best-model --min-data-in-leaf=2 " \
           "--model-format CatboostBinary,CPP " if catboost

if ARGV.empty?
  $stderr << "Insufficient arguments\n"
  usage(1)
end

fingerprints = if cmdline.option_present('gfp')
                 cmdline.value('gfp')
               else
                 '-EC3:ACHRY'
               end
flatten_sparse_fingerprints = cmdline.option_present('flatten')

smiles = ARGV[0]

train_smi = "#{mdir}/train.smi"
train_gfp = "#{mdir}/train.gfp"
train_activity = "#{mdir}/train.activity"

cmd = "#{gfp_make} #{fingerprints} #{smiles} > #{train_gfp}"
execute_cmd(cmd, verbose, [train_gfp])

FileUtils.cp(smiles, train_smi)

if cmdline.option_present('C')  # Classification.
  if lightgbm || catboost || xgboost
    perform_class_label_translation_lightgbm(activity_file, mdir, train_activity, verbose)
    lightgbm = "#{lightgbm} objective=binary" if lightgbm
    catboost = "#{catboost} --loss-function Logloss --custom-metric=MCC --auto-class-weights Balanced" if catboost
    xgboost = "#{xgboost} objective=binary:logistic" if xgboost
  else
    perform_class_label_translation(activity_file, mdir, train_activity, verbose)
    svm_learn_options = "#{svm_learn_options} -z c"
  end
else  # Regression
  perform_response_scaling(activity_file, mdir, train_smi, train_activity, verbose)
  svm_learn_options = "#{svm_learn_options} -z r"
  lightgbm = "#{lightgbm} objective=regression" if lightgbm
  catboost = "#{catboost} --loss-function RMSE" if catboost
  xgboost = "#{xgboost} objective=reg:squarederror" if xgboost
end

bit_xref = "bit"
bit_subset = "bit"

f = if flatten_sparse_fingerprints
      '-f'
    else
      ''
    end

l = if lightgbm || catboost || xgboost
      '-l'
    else
      ''
    end

cmd = "#{gfp_to_svm_lite} #{f} #{l} -C #{mdir}/#{bit_xref} -A #{train_activity} -A sep=auto -S #{mdir}/train "
if cmdline.option_present('p')
  support = cmdline.value('p')
  cmd = "#{cmd} -p #{support}"
end

train_svml = "#{mdir}/train.svml"
cmd = "#{cmd} #{train_gfp}"
execute_cmd(cmd, verbose, [train_svml, "#{mdir}/bit_xref.dat", "#{mdir}/bit_subset.dat"])

xgb_test = nil
if lightgbm
  model_file = "#{mdir}/LightGBM_model.txt"
  cmd = "#{lightgbm} data=#{mdir}/train.svml output_model=#{model_file}"
elsif catboost
  model_file = "#{mdir}/Catboost.model.bin"
  cmd = "#{catboost} --learn-set libsvm://#{mdir}/train.svml --model-file Catboost.model.bin"
elsif xgboost
  # need to interpolate the number of training rounds
  uri = File.absolute_path(train_svml)
  cmd = "#{xgboost} data=#{uri}?format=libsvm"
  if cmdline.option_present('xgb_test')
    test_fname = cmdline.value('xgb_test')
    test_svml = Tempfile.new('xgboost_test')
    cmd_test = "#{gfp_make} #{fingerprints} #{test_fname} | #{gfp_to_svm_lite} -l -X #{mdir}/bit_xref.dat -S #{test_svml.path} -"
    xgb_test = "#{test_svml.path}.svml"
    execute_cmd(cmd_test, verbose, [xgb_test])
    cmd = "#{cmd} test:data=#{xgb_test}?format=libsvm"
  end
else
  model_file = "#{mdir}/train.model"
  cmd = "#{svm_learn} #{svm_learn_options} #{train_svml} #{model_file}"
end
execute_cmd(cmd, verbose, [model_file])

# The name of the model file for xgboost is not known until after
# the model is built. We could also try to discern the number of
# boosting cycles.
xgboost_model_file = nil
if xgboost
  xgboost_model_file = get_xgboost_file(mdir)
end

# The metadata attribute is common among all model types.
def populate_metadata(model, fingerprints, response_name, classification, flatten_sparse_fingerprints)
  model.metadata = GfpModel::ModelMetadata.new
  model.metadata.date_built = Time.now.to_s
  model.metadata.fingerprints = fingerprints
  model.metadata.response_name = response_name
  if classification
    model.metadata.class_label_translation = 'class_label_translation.dat'
  else
    model.metadata.response_scaling = 'response_scaling.dat'
  end
  model.metadata.flatten_sparse_fingerprints = flatten_sparse_fingerprints
end

response_name = get_response_name(train_activity)

catboost = "#{catboost} --name=#{response_name}" if catboost

# Write the model proto
if lightgbm
  model = GfpModel::LightGbmModel.new
  populate_metadata(model, fingerprints, response_name, cmdline.option_present('C'), flatten_sparse_fingerprints)
  model.bit_xref = 'bit_xref.dat'
  model.training_cmdline = lightgbm
  File.write("#{mdir}/model.dat", GfpModel::LightGbmModel.encode(model))
  File.write("#{mdir}/model.json", GfpModel::LightGbmModel.encode_json(model))
elsif catboost
  model = GfpModel::CatboostModel.new
  populate_metadata(model, fingerprints, response_name, cmdline.option_present('C'), flatten_sparse_fingerprints)
  model.bit_xref = 'bit_xref.dat'
  model.training_cmdline = catboost
  File.write("#{mdir}/model.dat", GfpModel::CatboostModel.encode(model))
  File.write("#{mdir}/model.json", GfpModel::CatboostModel.encode_json(model))
elsif xgboost
  model = GfpModel::XGBoostModel.new
  populate_metadata(model, fingerprints, response_name, cmdline.option_present('C'), flatten_sparse_fingerprints)
  model.bit_xref = 'bit_xref.dat'
  model.config_file = xgboost_config
  model.training_cmdline = xgboost
  model.model_file = xgboost_model_file
  File.write("#{mdir}/model.dat", GfpModel::XGBoostModel.encode(model))
  File.write("#{mdir}/model.json", GfpModel::XGBoostModel.encode_json(model))
else
  support_vectors = "#{mdir}/support_vectors.gfp"
  cmd = "svm_model_support_vectors.sh -o #{support_vectors} #{model_file} #{train_gfp}"
  execute_cmd(cmd, verbose, [support_vectors])

  model = GfpModel::SvmfpModel.new
  populate_metadata(model, fingerprints, response_name, cmdline.option_present('C'), flatten_sparse_fingerprints)

  model.threshold_b = get_threshold_b(model_file)
  model.bit_subset = 'bit_subset.dat'
  model.bit_xref = bit_xref
  model.train_gfp = 'train.gfp'
  model.support_vectors = 'support_vectors.gfp'
  File.write("#{mdir}/model.dat", GfpModel::SvmfpModel.encode(model))
  File.write("#{mdir}/model.json", GfpModel::SvmfpModel.encode_json(model))
end

if cmdline.option_present('rescore')
  train_pred = File.join(mdir, 'train.pred')
  cmd = "svmfp_evaluate.sh -mdir #{mdir} #{mdir}/train.smi > #{train_pred}"
  execute_cmd(cmd, verbose, [train_pred])

  scaling_file = File.join(mdir, 'output_scaling')
  cmd = "iwstats -w -Y allequals -E #{activity_file} -p 2 -C #{scaling_file} #{train_pred}"
  execute_cmd(cmd, verbose, ["#{scaling_file}.dat"])
end
