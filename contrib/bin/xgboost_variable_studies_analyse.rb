#!/usr/bin/env ruby
# A series of models, including feature importance data, have been built, likely by
# xgboost_variable_studies.rb. Scan the feature importance files and work out a
# prioritised variable ordering.

require_relative('lib/iwcmdline')

def usage
  exit(0)
end

# Fetch the names of a series of directories whose names begin with `stem` and
# end with a sequential integer.
# Return an array of directory names.
def get_models(stem)
  ndx = 0
  result = []
  while true
    fname = "#{stem}#{ndx}"
    # $stderr << "Looking for #{fname} isdir #{File.directory?(fname)}\n"
    return result unless File.directory?(fname)
    return result unless File.size?(File.join(fname, 'feature_importance.txt'))
    result << fname
    ndx += 1
  end
end

def write_ranks(ranks, output)
  if ranks.empty?
    output << " . . . ."
    return
  end

  sum = ranks.sum
  mean = ranks.sum / ranks.size.to_f
  output << ' ' << ranks.size() << ' ' << mean.round(3)
  ranks.minmax.each do |r|
    output << ' ' << r
  end
end

# Given an array of scores, write N, mean, min and max to `output`.
def write_scores(scores, output)
  if scores.empty?
    output << " . . . ."
    return
  end

  sum = scores.sum
  mean = scores.sum / scores.size.to_f
  output << ' ' << scores.size() << ' ' << mean.round(3)
  scores.minmax.each do |v|
    output << ' ' << v.round(3)
  end
end

class Feature attr_accessor :average_combined_rank, :name, :times_seen
  def initialize(n)
    @name = n
    @times_seen = 0
    @rf_scores = []
    @rf_ranks = []
    # The average rank across RF and XGBoost
    @average_combined_rank
    @xgbd_scores = []
    @xgbd_ranks = []

    @times_top_rf_rank = 0
    @times_top_xgbd_rank = 0
  end

  def compute_average_combined_rank
    sum_ranks = @rf_ranks.sum + @xgbd_ranks.sum + @xgbd_ranks.sum
    if sum_ranks == 0
      @average_combined_rank = 0
    else
      @average_combined_rank = sum_ranks / (@rf_ranks.size + @xgbd_ranks.size).to_f
    end
  end

  def another_top_rf_score
    @times_top_rf_rank += 1
  end
  def another_top_xgbd_score
    @times_top_xgbd_rank += 1
  end

  def extra_rf(score, rank)
    @rf_scores << score
    @rf_ranks << rank
    @times_seen += 1
  end

  def extra_xgbd(score, rank)
    @xgbd_scores << score
    @xgbd_ranks << rank
    @times_seen += 1
  end

  def highest_xgbd_rank
    @xgbd_ranks.max
  end

  def write_summary(process_rf_data, output)
    output << "#{@name}"

    output << " #{@times_top_rf_rank}" if process_rf_data

    output << " #{@xgbd_ranks.min}"
    output << " #{@times_top_xgbd_rank}"

    output << " #{@average_combined_rank.round(3)}" if process_rf_data

    write_scores(@rf_scores, output) if process_rf_data
    write_scores(@xgbd_scores, output)

    write_ranks(@rf_ranks, output) if process_rf_data
    write_ranks(@xgbd_ranks, output)
    output << "\n"
  end
end

# Return the class of feature 'w_natoms' -> 'w'
def class_name(feature)
  return feature.gsub(/_.*/, "")
end

# Read a feature importance file and accumulate statistics about
# individual features, by_feature, and feature classes, by_feature.
def get_scores(dirname, is_rf, by_feature, by_class)
  highest_score = 0.0
  fname = File.join(dirname, 'feature_importance.txt')
  feature_names = []
  File.foreach(fname, chomp:true).with_index do |line, line_number|
    if line_number == 0  # skip header
      next
    end

    f = line.split
    feature_names << f[0]
    score = f[1].to_f
    highest_score = score if line_number == 1

    to = if by_feature.key?(f[0])
           by_feature[f[0]]
         else
           by_feature[f[0]] = Feature.new(f[0])
         end

    fclass = class_name(f[0])
    by_class_to = if by_class.key?(fclass)
                    by_class[fclass]
                  else
                    by_class[fclass] = Feature.new(fclass)
                  end

    if is_rf
      to.extra_rf(score/highest_score, line_number - 1)
      by_class_to.extra_rf(score/highest_score, line_number - 1)
    else
      to.extra_xgbd(score/highest_score, line_number - 1)
      by_class_to.extra_xgbd(score/highest_score, line_number - 1)
    end
  end

  # $stderr << " Have #{feature_names.size} features\n"
  # $stderr << "#{by_feature.size} by feature and #{by_class.size} by class\n"

  (0...feature_names.size / 10).each do |i|
   ftype = class_name(feature_names[i])

    if is_rf
      by_feature[feature_names[i]].another_top_rf_score
      by_class[ftype].another_top_rf_score
    else
      by_feature[feature_names[i]].another_top_xgbd_score
      by_class[ftype].another_top_xgbd_score
    end
  end
end

def write_header(process_rf_data, output)
  output << "Id"
  if process_rf_data
    output << " RFtop10"
  end
  output << " XGTopRank XGTop10"
  if process_rf_data
    output << " AveRank"
  end
  if process_rf_data
    output << " NRF RFavescore RFminscore RFmaxscore"
  end
  output << " NXG XGavescore XGminscore XGmaxscore"
  if process_rf_data
    output << " NRF RFaverank RFminrank RFmaxrank"
  end
  output << " NXG XGaverank XGminrank XGmaxrank"

  output << "\n"
end

def main
  cl = IWCmdline.new('-v-stem=s-RF-iwcut=s-support=ipos')

  verbose = cl.option_present('v')

  stem = if cl.option_present('stem')
           cl.value('stem')
         else
           'xgbd_vstudy'
         end


  support = if cl.option_present('support')
              cl.value('support')
            else
              0
            end

  process_rf_data = cl.option_present('RF')

  if process_rf_data
    rf_models = get_models("#{stem}.RF")
    raise "No RF models" unless rf_models.size > 0
  else
    rf_models = []
  end

  xgbd_models = get_models("#{stem}.XGBD")

  if xgbd_models.empty?
    $stderr << "No xgboost models with stem #{stem}\n"
    return 1;
  end

  $stderr << "Read #{rf_models.size} RF models and #{xgbd_models.size} XGBoost models\n" if verbose

  # Individual features, w_natoms
  features = {}
  # Feature class, all w_* descriptors
  by_class = {}

  rf_models.each do |m|
    get_scores(m, true, features, by_class)
  end
  xgbd_models.each do |m|
    get_scores(m, false, features, by_class)
  end

  $stderr << "Read data on #{features.size} descriptors\n" if verbose

  features.each do |k, v|
    v.compute_average_combined_rank
  end
  by_class.each do |k, v|
    v.compute_average_combined_rank
  end

  f = features.values.sort { |f1, f2| f1.average_combined_rank <=> f2.average_combined_rank }

  discarded_lack_of_support = 0
  write_header(process_rf_data, $stdout)
  f.each do |feature|
    if feature.times_seen < support
      discarded_lack_of_support += 1
      next
    end
    feature.write_summary(process_rf_data, $stdout)
  end

  if support > 0
    $stderr << "#{discarded_lack_of_support} of #{f.size} features skipped for failing support #{support}\n"
  end

  f = by_class.values.sort { |f1, f2| f1.average_combined_rank <=> f2.average_combined_rank }

  $stdout << "\n"
  write_header(process_rf_data, $stdout)
  f.each do |fclass|
    fclass.write_summary(process_rf_data, $stdout)
  end

  if cl.option_present('iwcut')
    fname = cl.value('iwcut')
    f = features.values.sort { |f1, f2| f1.average_combined_rank <=> f2.average_combined_rank }
    File.open(fname, 'w') do |writer|
      f.each do |feature|
        next if feature.times_seen < support
        writer << feature.name << "\n"
      end
    end
  end
end


main
