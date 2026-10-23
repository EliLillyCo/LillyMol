#!/usr/bin/env ruby
# Build multiple svmfp models based on a set of fingerprints.
# Useful for figuring out the best fingerprint to use for a given target.

require_relative 'lib/iwcmdline.rb'

def usage
  $stderr << "Build svmfp models\n"
  $stderr << " -fp <fname>          file of fingerprints to try\n"
  $stderr << " -nsplit <nsplit>     number of splits to create\n"
  exit(1)
end

cl = IWCmdline.new('-v-trpct=ipos-nsplit=ipos-niter=ipos-A=sfile-fp=sfile-catboost=close-xgboost=close-xgboost_config=sfile-lightgbm=close-lightgbm_config=sfile-i=ipos-extrasplit=s')

if cl.unrecognised_options_encountered
  $stderr << "unrecognised_options_encountered\n"
  usage
end

verbose = cl.option_present('v')

if ARGV.empty?
  $stderr << "Must specify smiles file\n"
  usage
end

smiles = ARGV[0]

if ! cl.option_present('trpct')
  $stderr << "Must specify the training set percent via the -trpct option\n"
  usage
end

trpct = 80
prefix = 'A80'
if cl.option_present('trpct')
  x = cl.value('trpct')
  if x < 0 || x > 100
    $stderr << "The training percent (-trpct) option must be a valid percent\n"
    usage
  end
  prefix = "A#{x}"
  # trpct = x.to_i / 100.0
end

lightgbm = cl.values('lightgbm')
lightgbm_config = cl.value('lightgbm_config')
catboost = cl.values('catboost')
xgboost = cl.values('xgboost')
xgboost_config = cl.value('xgboost_config')

if ! cl.option_present('fp')
  $stderr << "Must specify file of fingerprints via the -fp option\n"
  usage
end

fingerprints = File.readlines(cl.value('fp'))

if ! cl.option_present('A')
  $stderr << "Must specify activity file via the -A option\n"
  usage
end

activity_fname = cl.value('A')

def execute_cmd(cmd, verbose, expected_files)
  $stderr << "Executing '#{cmd}'\n" if verbose
  system(cmd)
  expected_files.each do |fname|
    if ! File.size?(fname)
      "$stderr << "#{cmd} did not create #{fname}\n"
      return false
    end
  end
  return true
end

# The -ps option has been given. We need to determine the files implied by this.
# Return a list of the train and test splits.
def get_pre_split(stem)
  
end


if cl.option_present('ps')
  train_files, test_files = get_pre_split(cl.value('ps'))
  nsplit = train_files.size
else
  if cl.option_present('nsplit')
    nsplit = cl.value('nsplit')
  elsif cl.option_present('niter')
    nsplit = cl.value('niter')
  else
    nsplit = 10
  end

  cmd = "stratified_sample -A #{activity_fname} -n #{nsplit} -R #{prefix}_train -E #{prefix}_test #{smiles}"
  cmd = "stratified_samples -N #{nsplit} -R #{prefix}_TRAIN -E #{prefix}_TEST -p #{trpct} -s 1 -M #{smiles} #{activity_fname}"
  execute_cmd(cmd, verbose, [])
  train_files = (0..nsplit).map { |i| "#{prefix}_TRAIN#{i}.smi"}
  test_files = (0..nsplit).map { |i| "#{prefix}_TEST#{i}.smi"}
end

if cl.option_present('extrasplit')
  cl.values('extrasplit').each do |split|
    train_fname = "#{split}R.smi"
    test_fname = "#{split}E.smi"
    if File.size?(train_fname) && File.size?(test_fname)
      train_files << train_fname
      test_files << test_fname
      nsplit += 1
    else
      raise "Missing or empty file(s) #{train_fname} #{test_fname}"
    end
  end
end

$stderr << train_files << ' ' << test_files << "\n";

stem = if cl.option_present('stem')
    cl.value('stem')
  else
    'model'
  end

support = if cl.option_present('p')
    cl.value('p')
  else
    1
  end

fingerprints.each do |fp|
  fp.chomp!
  fps = fp.gsub(/ /, '')
  preds = []
  (0...nsplit).each do |split|
    mdir = "#{stem}_#{fps}_#{split}"
    cmd = "svmfp_make.sh -p #{support} -mdir #{mdir} -gfp #{fp} -gfp -A #{activity_fname} #{train_files[split]}"
    execute_cmd(cmd, verbose, [])
    pred = "#{stem}_#{fps}_#{split}.pred"
    cmd = "svmfp_evaluate.sh -mdir #{mdir} #{test_files[split]} > #{pred}"
    execute_cmd(cmd, verbose, [pred])
    preds << pred
    catboost.each_with_index do |cb, ndx|
      c_mdir = "#{mdir}_catboost_#{ndx}"
      cmd = "svmfp_make.sh -mdir #{c_mdir} -gfp #{fp} -gfp -catboost #{cb} -catboost -A #{activity_fname} #{train_files[split]}"
      execute_cmd(cmd, verbose, [])
      pred = "#{stem}_#{fps}_cb_#{split}_#{ndx}.pred"
      cmd = "catboost_evaluate.sh -mdir #{c_mdir} #{test_files[split]} > #{pred}"
      execute_cmd(cmd, verbose, [pred])
      preds << pred
    end
    xgboost.each_with_index do |xg, ndx|
      c_mdir = "#{mdir}_xgboost"
      cmd = "svmfp_make.sh -mdir #{c_mdir} -gfp #{fp} -gfp -xgboost -xgboost -xgboost_config #{xgboost_config} -A #{activity_fname} #{train_files[split]}"
      execute_cmd(cmd, verbose, [])
      pred = "#{stem}_#{fps}_xg_#{split}_#{ndx}.pred"
      cmd = "xgboost_evaluate.sh -v -mdir #{c_mdir} #{test_files[split]} > #{pred}"
      execute_cmd(cmd, verbose, [pred])
      preds << pred
    end
    lightgbm.each_with_index do |lg, ndx|
      c_mdir = "#{mdir}_lightgbm_#{ndx}"
      cmd = "svmfp_make.sh -mdir #{c_mdir} -gfp #{fp} -gfp -lightgbm #{lg} -lightgbm -lightgbm_config #{lightgbm_config} -A #{activity_fname} #{train_files[split]}"
      execute_cmd(cmd, verbose, [])
      pred = "#{stem}_#{fps}_lg_#{split}_#{ndx}.pred"
      cmd = "lightgbm_evaluate.sh -mdir #{c_mdir} #{test_files[split]} > #{pred}"
      execute_cmd(cmd, verbose, [pred])
      preds << pred
    end
  end
  # ifile = "I#{stem}_#{fps}"
  # cmd = "iwstats -w -Y allequals -E #{activity_fname} -p 2 -I #{ifile} #{preds.join(' ')}"
  # stats_file = "#{stem}_#{fps}.#{cycle}"
  # cmd = "iwstats -w -Y allequals -E #{activity_fname} -p 2 #{preds.join(' ')}"
  # execute_cmd(cmd, verbose, [ifile])
end
