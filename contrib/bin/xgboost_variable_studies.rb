#!/usr/bin/env ruby

require 'fileutils'

require_relative('lib/iwcmdline')

def usage
  help = <<-END
Sets up jobs for analysing XGBoost feature importance.
Observe that XGBoost feature importance can be a strong function of the order in which columns are presented
to the model.
This tool performs -niter random permutations of the columns in a descriptor file.
For each permutation an XGBoost model is built, and the feature importance recorded.
A small number of random forest models are also built.
The resulting models and feature importance files can be analysed with xgboost_variable_studies_analyse
Usage:

 -A <fname>             Activity file - mandatory.
 -niter <n>             Number of XGBoost models to build (default 10).
 -output <fname>        Outut file to create - contains commands for batch computation.
 -score <fname>         If you wish to score a test set with each model, specify here. Optional.
 -submit                Submit the job to the cluster with iwqb.rb - otherwise you can examine and submit.
 -v                     verbose output.
END

  $stderr << help

  exit(0)
end

def descriptor_file_name(name_stem, ndx)
  return "#{name_stem}#{ndx + 1}.dat"
end

# Return a string for building a couple of RF models.
def build_rf_model(replicates, activity_file, name_stem)
  result = []
  replicates.times do |i|
    mdir = "#{name_stem}.RF#{i}"
    descriptor_file = descriptor_file_name(name_stem, i)
    $stderr << "Descriptor file is #{descriptor_file}\n"
    result << "rf_make.sh --mdir #{mdir} --activity #{activity_file} --feature_importance #{descriptor_file}\n"
  end

  return result
end

# Return a set of commands to self score the models built.
def self_score(name_stem, niter)
  result = []
  niter.times do |i|
    fname = descriptor_file_name(name_stem, i)
    mdir = "#{name_stem}.XGBD#{i}"
    pred = "#{mdir}.pred"
    cmd = "xgbd_evaluate.sh -mdir #{mdir} #{fname} > #{pred}"
    result << cmd
  end

  return result
end

def main
  cl = IWCmdline.new('-v-niter=ipos-A=sfile-output=s-score=sfile-submit')

  verbose = cl.option_present('v')

  unless cl.option_present('A')
    $stderr << "Must specify activity file via the -A option\n"
    usage
  end

  activity_file = cl.value('A')

  if ARGV.empty?
    $stderr << "Must specify input descriptor file\n"
    usage
  end

  niter = if cl.option_present('niter')
            cl.value('niter')
          else
            10
          end

  descriptor_file = "#{ARGV[0]}"
  unless File.size?(descriptor_file)
    $stderr << "Missing or empty input file #{descriptor_file}\n"
    exit(1)
  end

  name_stem = 'xgbd_vstudy'
  cmd = "random_column_order.sh -N #{niter} -j -S #{name_stem} -U .dat #{descriptor_file}"
  $stderr << "Executing #{cmd}\n" if verbose
  system(cmd)

  # RF models seem not to have much variability across models, so only 2 replicates.
  rf_niter = 0  # Turned off, creates too many problems.

  test_set = if cl.option_present('score')
                 cl.value('score')
               else
                 nil
               end

  cmds = []

  niter.times do |i|
    mdir = "#{name_stem}.XGBD#{i}"
    fname = "#{name_stem}#{i + 1}.dat"
    raise "No #{fname}" unless File.size?(fname)
    cmds << "/lrlhps/users/rx87690/bin.2025/xgbd_make.sh --mdir #{mdir} --activity #{activity_file} --feature_importance def --nthreads 1 #{fname}"

    if test_set
      cmds << "xgbd_evaluate.sh -mdir #{mdir} #{test_set} > #{name_stem}#{i}.pred"
    end
    cmds << "|\n"   # enable multi-line jobs in iwqb
  end

  # Rf last deliberately
  cmds << build_rf_model(rf_niter, activity_file, name_stem) if rf_niter > 0

# $stderr << cmds.join("\n") << "\n"

  cmd_file = if cl.option_present('output')
               cl.value('output')
             else
               'xgbd_vstudy.sh'
             end

  File.write(cmd_file, cmds.join("\n") << "\n")

  if cl.option_present('submit')
    cmd = "iwqb.sh -m -qsub -V -qsub #{cmd_file}"
    $stderr << "Submitting #{cmd}\n"
    system("iwqb.sh -m -qsub -V -qsub #{cmd_file}")
  else
    $stderr << "Submit job with      iwqb.sh -m -qsub -V -qsub #{cmd_file}\n"
  end

# if cl.option_present('self_score')
#   cmds << self_score(name_stem, niter)
# end
end

main
