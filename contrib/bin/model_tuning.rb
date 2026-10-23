#!/usr/bin/env ruby
# Build svmfp models based on a pre-split set of files
# and a set of fingerprints.

require_relative 'lib/iwcmdline.rb'

def usage(cl)
  expert = cl.option_present('expert')

  $stderr << "Build svmfp models with different fingerprints\n"
  $stderr << " -A <fname>           activity file.\n"
  $stderr << " -trpct <pct>         percent of dataset to use for training\n"
  $stderr << " -fp <fname>          file of fingerprints to try\n"
  $stderr << " -fpextra <fp> -fpextra   duplicate the fingerprint array and add <fp> to each\n"
  $stderr << " -nsplit <nsplit>     number of splits to create\n"
  $stderr << " -keep_models         do NOT remove model directories\n"
  $stderr << " -chrono              make an extra chronological split - uses numeric values in ids\n"
  $stderr << " -iwstats ... -iwstats passed to iwstats\n" if expert
  $stderr << " -svml ... -svml      passed to svmfp_make\n" if expert
  $stderr << " -expert              more options\n" unless expert
  $stderr << " -v                   verbose output\n"

  exit(1)
end

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
# Return arrays of the train and test smiles files.
def get_pre_split(train_stem, test_stem)
  # We return these two arrays.
  train_files = []
  test_files = []
  0.upto(1000) do |ndx|
    train_smi = "#{train_stem}#{ndx}.smi"
    test_smi = "#{test_stem}#{ndx}.smi"
    if File.size?(train_smi) && File.size?(test_smi)
      train_files << train_smi
      test_files << test_smi
      next
    elsif ndx == 0
      # Ok to miss index 0, maybe the numbering sequence starts with 1
    else
      return train_files, test_files
    end
  end

  return train_files, test_files
end

# Given a set of training set smiles and a set of test set files, return
# an array of training set descriptor files and an array of test set
# descriptor files
def fetch_split_descriptor_files(train_smiles, test_smiles, haystack, verbose)
  train_files = []
  test_files = []
  raise "Inconsistent sizes" unless train_smiles.size == test_smiles.size
  train_smiles.each_with_index do |smiles_file, ndx|
    output = smiles_file.gsub(/\.smi$/, '.dat')
    execute_cmd("descriptor_file_select_rows -c 2 #{smiles_file} #{haystack} > #{output}", verbose, [output])
    train_files << output

    smiles_file = test_smiles[ndx]
    output = smiles_file.gsub(/\.smi$/, '.dat')
    execute_cmd("descriptor_file_select_rows -c 2 #{smiles_file} #{haystack} > #{output}", verbose, [output])
    test_files << output
  end

  return train_files, test_files
end

def model_tuning
  cl = IWCmdline.new('-v-trpct=ipos-nsplit=ipos-niter=ipos-A=sfile-fp=sfile-fpextra=close' +
                     '-catboost=close-xgboost=close-lightgbm=close-lightgbm_config=sfile' +
                     '-dfile=sfile-ps-chrono-iwstats=close-expert' +
                     '-sameline' +
                     '-i=ipos-svml=close-keep_models')

  if cl.unrecognised_options_encountered
    $stderr << "unrecognised_options_encountered\n"
    usage(cl)
  end

  verbose = cl.option_present('v')

  if ARGV.empty?
    $stderr << "Must specify smiles file\n"
    usage(cl)
  end

  smiles = ARGV[0]
  raise "Missing or empty smiles file #{smiles}" unless File.size?(smiles)

  trpct = 80
  if cl.option_present('trpct')
    trpct = cl.value('trpct')
    if trpct < 0 || trpct > 100
      $stderr << "The training percent (-trpct) option must be a valid percent\n"
      usage(cl)
    end
  end

  prefix = "A#{trpct}"

  # Not fully implemented...
  classification = cl.option_present('C')

  lightgbm = cl.values('lightgbm')
  lightgbm_config = cl.value('lightgbm_config')
  catboost = cl.values('catboost')
  xgboost = cl.values('xgboost')

  fingerprint_files = []
  if cl.option_present('fp')
    fingerprint_files = cl.values('fp')
  else
    lillymol_home = ENV['LILLYMOL_HOME']
    fp = File.join(lillymol_home, 'contrib', 'data', 'default_fingerprints')
    raise "Missing or empty default_fingerprints #{fp}" unless File.size?(fp) > 0
    fingerprint_files << fp
  end

  fingerprints = []
  fingerprint_files.each do |fname|
    fingerprints.push(*File.readlines(fname).map { |fp| fp.chomp})
  end

  cl.values('fpextra').each do |fpextra|
    extrafp = fingerprints.map { |fp| "#{fp} #{fpextra}" }
    fingerprints.push(*extrafp)
  end

  $stderr << fingerprints << "\n" if verbose

  unless cl.option_present('A')
    $stderr << "Must specify activity file via the -A option\n"
    usage(cl)
  end

  activity_fname = cl.value('A')

  train_stem = 'TRAIN'
  test_stem = 'TEST'

  if cl.option_present('ps')
    train_files, test_files = get_pre_split(train_stem, test_stem)
    nsplit = train_files.size
    $stderr << "Find #{nsplit} previous splits\n" if verbose
    exit 1 if nsplit == 0
  else
    if cl.option_present('nsplit')
      nsplit = cl.value('nsplit')
    elsif cl.option_present('niter')
      nsplit = cl.value('niter')
    else
      nsplit = 10
    end

    cmd = "stratified_samples -s 1 -N #{nsplit} -p #{trpct} -R #{train_stem} -E #{test_stem} -M #{smiles}"
    if cl.option_present('chrono')
      cmd << ' -C'
      nsplit += 1
    end
    cmd << " #{activity_fname}"
    execute_cmd(cmd, verbose, ["#{train_stem}0.smi", "#{test_stem}0.smi"])
    train_files = (0..nsplit).map { |i| "#{train_stem}#{i}.smi"}
    test_files = (0..nsplit).map { |i| "#{test_stem}#{i}.smi"}
  end

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

  # Deliberate limitation, we only handle one descriptor file
  train_x = []
  train_y = []
  dfile = cl.values('dfile')
  if dfile.size > 0
    train_x, test_x = fetch_split_descriptor_files(train_files, test_files, dfile[0], verbose)
  end

  keep_models = cl.option_present('keep_models')

  # MayBe NewLine
  mbnl = if cl.option_present('sameline')
           mbnl = ";"
         else
           mbnl = "\n"
         end

  cmd_stream = File.open("model_tuning.txt", "w")

  if classification
    iwstats = 'confusion_matrix.rb'
  else
    iwstats = "iwstats -w -Y allequals -p 2 -E #{activity_fname}"
    iwstats << ' ' << cl.value('iwstats') if cl.option_present('iwstats')
  end

  svmfp_make = "time svmfp_make.sh -A #{activity_fname}"
  svmfp_make << " -p #{support}" if support > 1
  svmfp_make << " -C" if classification
  svmfp_make << ' ' << cl.value('svml') if cl.option_present('svml')

  fingerprints.each do |fp|
    fps = fp.gsub(/ /, '')
    preds = []
    (0...nsplit).each do |split|
      mdir = "#{stem}.#{fps}.#{split}"
      cmd_stream << "#{svmfp_make} -mdir #{mdir} -gfp #{fp} -gfp #{train_files[split]}#{mbnl}"
      pred = "#{stem}_#{fps}_#{split}.pred"
      cmd_stream << "svmfp_evaluate.sh -mdir #{mdir} #{test_files[split]} > #{pred}#{mbnl}"
      preds << pred
      cmd_stream << "#{iwstats} #{pred} > #{prefix}.#{fps}.#{split}#{mbnl}"
      cmd_stream << "rm -r #{mdir}#{mbnl}" unless keep_models
      cmd_stream << "\n"

      catboost.each_with_index do |cb, ndx|
        c_mdir = "#{mdir}_catboost_#{ndx}"
        cmd_stream << "svmfp_make.sh -mdir #{c_mdir} -gfp #{fp} -gfp -catboost #{cb} -catboost -A #{activity_fname} #{train_files[split]}\n"
        pred = "#{stem}_#{fps}_cb_#{split}_#{ndx}.pred"
        cmd_stream << "catboost_evaluate.sh -mdir #{c_mdir} #{test_files[split]} > #{pred}\n"
        preds << pred
        cmd_stream << "#{iwstats} #{pred} > #{prefix}.ctbfp.#{ndx}.#{split}\n"
        cmd_stream << "rm -r #{c_mdir}\n" unless keep_models
      end
      # Only one descriptor file is processed, but any number of xgboost configs.
      xgboost.each_with_index do |xg, ndx|
        c_mdir = "#{stem}_xgboost.#{ndx}.#{split}"
        cmd_stream << "xgbd_make.sh --mdir #{c_mdir} #{xg} --activity #{activity_fname} #{train_x[0]}\n"
        pred = "#{stem}_xg_#{ndx}_#{split}.pred"
        cmd_stream << "xgboost_evaluate.sh -v -mdir #{c_mdir} #{test_x[0]} > #{pred}\n"
        preds << pred
        cmd_stream << "#{iwstats} #{pred} > #{prefix}.xgbfp.#{ndx}.#{split}\n"
        cmd_stream << "rm -r #{c_mdir}\n" unless keep_models
      end
      lightgbm.each_with_index do |lg, ndx|
        c_mdir = "#{mdir}_lightgbm_#{ndx}"
        cmd_stream << "svmfp_make.sh -mdir #{c_mdir} -gfp #{fp} -gfp -lightgbm #{lg} -lightgbm -lightgbm_config #{lightgbm_config} -A #{activity_fname} #{train_files[split]}\n"
        pred = "#{stem}_#{fps}_lg_#{split}_#{ndx}.pred"
        cmd_stream << "lightgbm_evaluate.sh -mdir #{c_mdir} #{test_files[split]} > #{pred}\n"
        preds << pred
        cmd_stream << "#{iwstats} #{pred} > #{prefix}.lgbfp.#{ndx}.#{split}\n"
        cmd_stream << "rm -r #{c_mdir}\n" unless keep_models
      end
    end
  end
  #  ifile = "I#{stem}_#{fps}"
  # cmd = "iwstats -E #{activity_fname} -p 2 -I #{ifile} #{preds.join(' ')}"
  # execute_cmd(cmd, verbose, [ifile])
end

model_tuning()
