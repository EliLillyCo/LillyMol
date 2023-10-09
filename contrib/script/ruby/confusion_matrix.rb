#! /usr/bin/env ruby
# $Id$
ianhome = ENV['C3TK_BIN']
# ianhome = "/home/rx87851"

require "#{ianhome}/ruby/lib/iwcmdline.rb"

$expert = false

def usage (rc)
  $stderr.print "Reports classification accuracies\n"
  $stderr.print " -E <fname>     experimental values file\n"
  $stderr.print " -e <col>       experimental values in column <col>\n"
  $stderr.print " -P <fname>     predicted values file\n"
  $stderr.print " -p <col>       predicted values column <col>\n"
  $stderr.print " -pipe          read predicted values from stdin\n" if ($expert)
  $stderr.print " -pos <class>   compute PPV and NPV measures assuming 'class' is the positive class\n"
  $stderr.print " -n <c>         convert floating point data to class labels, < c and > c\n"
  $stderr.print " -n <l,h>       convert floating point data to class labels, < l and > h\n"
  $stderr.print " -En <c>        convert floating point data to class labels in expt data\n" if ($expert)
  $stderr.print " -Pn <c>        convert floating point data to class labels in pred data \n" if ($expert)
  $stderr.print " -nm            suppress display of the matrix\n"
  $stderr.print " -sum           write sum of all measures (SENS + PPV + SPEC + NPV)\n"
  $stderr.print " -v             verbose output\n"
  exit(rc)
end

cl = IWCmdline.new("-v-expert-P=sfile-p=ipos-E=sfile-e=ipos-n=s-cutoff=f-N-pos=s-En=f-Pn=f-nm-sum-pipe")

$expert = cl.option_present('expert')

if cl.unrecognised_options_encountered()
  $stderr.print "Unrecognised options encountered\n"
  usage(1)
end

verbose = cl.option_present('v')

$convert_experimental_to_two_class = false
$convert_predicted_to_two_class_lower_bound = false
$convert_predicted_to_two_class_upper_bound = false

if (cl.option_present('n'))
  n = cl.value('n')
  comma_separated = Regexp.new('(\S+),(\S+)')
  m = comma_separated.match(n)
  if (m)
    $convert_experimental_to_two_class = 0.0   # not sure what else to do
    $convert_predicted_to_two_class_lower_bound = m[1].to_f
    $convert_predicted_to_two_class_upper_bound = m[2].to_f
    if ($convert_predicted_to_two_class_lower_bound > $convert_predicted_to_two_class_upper_bound)
      $stderr << "Lower bound #{$convert_predicted_to_two_class_upper_bound} greater than upper bound #{$convert_predicted_to_two_class_lower_bound}, impossible\n"
      exit 2
    end
  else
    $convert_experimental_to_two_class = n.to_f
    $convert_predicted_to_two_class_lower_bound = n.to_f
    $convert_predicted_to_two_class_upper_bound = n.to_f
  end
end

if (cl.option_present('En'))
  $convert_experimental_to_two_class = cl.value('En')
end

if (cl.option_present('Pn'))
  $convert_predicted_to_two_class_lower_bound = cl.value('Pn')
  $convert_predicted_to_two_class_upper_bound = cl.value('Pn')
end

if ($convert_experimental_to_two_class && $convert_predicted_to_two_class_upper_bound)
  true    # both set, that's good
elsif ($convert_experimental_to_two_class || $convert_predicted_to_two_class_upper_bound)
  $stderr.print "Must use both of -En and -Pn\n"
  usage(4)
end

$convert_to_nearest_integer = cl.option_present('N')

if ($convert_experimental_to_two_class && $convert_to_nearest_integer)
  $stderr.print "Cannot combine two class processing with the -N option\n"
  usage(3)
end

$display_matrix = true

if (cl.option_present('nm'))
  $display_matrix = false
end

# By default, the input file will be 'id expt pred'

experimental_column = 2
predicted_column = 2

if (cl.option_present('e'))
  experimental_column = cl.value('e')
end

if (cl.option_present('p'))
  predicted_column = cl.value('p')
end

experimental_column -= 1
predicted_column -= 1

def do_convert_to_nearest_integer(v)
  if (v < 0.0)
    v = (v - 0.49999).to_i
  else
    v = (v  + 0.49999).to_i
  end
  
  return v
end

# Converts the experimental data - where there is only one cutoff

def convert_to_class_if_needed1(v, convert_to_two_class)
  if (convert_to_two_class)
    v = v.to_f
    if (v < convert_to_two_class)
      return '-1'
    elsif (v >= convert_to_two_class)
      return '1'
    else
      $stderr.print "What to do with experimental value '#{v}'\n";
      return v.to_s
    end
  elsif ($convert_to_nearest_integer)
    return do_convert_to_nearest_integer(v.to_f).to_s
  end

  v
# $stderr.print "Converted to '#{v}'\n";
end

# Converts the predicted data - where there might be two cutoffs

def convert_to_class_if_needed2(v, convert_to_two_class_lower_bound,
                               convert_to_two_class_upper_bound)
# $stderr << "Converting #{v}, between #{convert_to_two_class_lower_bound} and #{convert_to_two_class_upper_bound}\n"

  if (convert_to_two_class_lower_bound)
    v = v.to_f
    if (v < convert_to_two_class_lower_bound)
      return '-1'
    elsif (v >= convert_to_two_class_upper_bound)
      return '1'
    else
      return false
    end
  elsif ($convert_to_nearest_integer)
    return do_convert_to_nearest_integer(v.to_f).to_s
  end

# $stderr.print "Converted to '#{v}'\n";
  return v
end

# If there is a separate experimental file, read truth from there

experimental_value = Hash.new

if (cl.option_present('E'))
  experimental_file_name= cl.value('E')
  inp = File.open(experimental_file_name, mode="r")
  raise "Cannot open experiment data file '#{experimental_file_name}'" unless (inp)
  inp.each do |line|
    line.chomp!
    f = line.split

    id = f[0]
    expt = f[experimental_column]

    expt = convert_to_class_if_needed1(expt, $convert_experimental_to_two_class) 

    experimental_value[id] = expt

#   $stderr.print "ID '#{id}', expt '#{expt}'\n"
  end

  inp.close

  if (verbose)
    $stderr.print "Read #{experimental_value.size} experimental data values from '#{experimental_file_name}'\n"
  end
end

predicted_file = false
if (cl.option_present('P'))
  predicted_file = cl.value('P')
elsif (cl.option_present('pipe'))
  predicted_file = '-'
elsif (ARGV.size > 0)
  predicted_file = ARGV[0]
else
  $stderr.print "Insufficient arguments\n"
  usage(2)
end

header_records_to_skip = 1

inp = "-" == predicted_file ? $stdin : File.open(predicted_file, mode="r")
raise "Cannot open '#{predicted_file}'" unless (inp)

(0...header_records_to_skip).each do |i|
  line = inp.gets
# $stderr.print "Got header record\n";
end

$classes = Hash.new

class Aclass
  def initialize(nam)
    @_name = nam
    @_n = 0
    @_classified = Hash.new

    @sensitivity = false
    @ppv = false
    @specificity = false
    @npv = false
  end

  def class_name
    @_name
  end

  def found_in(nam)
    if (! @_classified.has_key?(nam))
      @_classified[nam] = 1
    else
      @_classified[nam] += 1
    end

    @_n += 1
  end

  def report(s)
    s.print "#{@_n} instances of class '#{@_name}', "

    fraction_class_correct = 0.0    # our return code

    if (@_classified.has_key?(@_name))
      n = @_classified[@_name]
      fraction_class_correct = n.to_f / @_n.to_f
      s.printf("correct #{n} %.3f\n", fraction_class_correct)
    else
      s.print "correct 0 0.0\n"
    end

    longest_class_name = 0
    $classes.each_value do |c|
      cname = c.class_name
      longest_class_name = cname.length if (cname.length > longest_class_name)
    end

    $classes.each_value do |c|
      cname = c.class_name
      if (@_classified.has_key?(cname))
        ratio = @_classified[cname].to_f / @_n.to_f
        s.printf(" %#{longest_class_name}s #{@_classified[cname]} %.2f\n", cname, ratio) if ($display_matrix)
      else
        s.printf(" %#{longest_class_name}s 0 0.0\n", cname) if ($display_matrix)
      end

    end

    return fraction_class_correct
  end

  def times_predicted_in_class(c)
    if (! @_classified.has_key?(c))
      return 0
    end

    @_classified[c]
  end

# We are the positive class in a two class problem. Compute various metrics

  def compute_ppv(negative_class)
    raise "Not a two class problem" if (@_classified.size > 2)

    true_positives = 0
    if (@_classified.has_key?(@_name))
      true_positives = @_classified[@_name]
    end

    false_positives = @_n - true_positives

    g = true_positives + negative_class.times_predicted_in_class(@_name)

    if (g > 0)
      @sensitivity = true_positives.to_f / @_n.to_f
      @ppv = true_positives.to_f / g.to_f
      $stdout.printf("Sensitivity %.3f PPV %.3f\n", @sensitivity, @ppv)
    else
      $stdout.printf("Sensitivity %.3f PPV 0\n", true_positives.to_f / @_n.to_f)
    end
  end

# We are the negative class in a two class problem. Compute various metrics

  def compute_npv(positive_class)
    raise "Not a two class problem" if (@_classified.size > 2)

    true_negatives = 0
    if (@_classified.has_key?(@_name))
      true_negatives = @_classified[@_name]
    end

    false_negatives = @_n - true_negatives

    h = positive_class.times_predicted_in_class(@_name) + true_negatives

    if (h > 0)
      @specificity = true_negatives.to_f / @_n.to_f
      @npv = true_negatives.to_f / h.to_f
      $stdout.printf("Specificity %.3f NPV %.3f\n", @specificity, @npv)
    else
      $stdout.print "Specificity undefined\n";
    end
  end

  def compute_tn_fn(positive_class)
    tn = 0
    if (@_classified.has_key?(@_name))
      tn = @_classified[@_name]
    end

    fn = @_n - tn

    return tn,fn
  end
  def compute_tp_fp(negative_class)
    tp = 0
    if (@_classified.has_key?(@_name))
      tp = @_classified[@_name]
    end

    fp = @_n - tp

    return tp,fp
  end

  def sum_measures_positive_class
    return false unless (@sensitivity && @ppv)
    return @sensitivity + @ppv
  end
  def sum_measures_negative_class
    return false unless (@specificity && @npv)
    return @specificity + @npv;
  end
end

def compute_matthews(tp, tn, fp, fn)
  m = (tp*tn - fp*fn).to_f / Math.sqrt((tp+fn)*(tn+fp)*(tp+fp)*(tn+fn))
  $stdout.printf("MatthewsCoefficient %.3f\n", m)
end

def compute_accuracy(tp, tn, fp, fn)
  acc = (tp + tn).to_f / (tp + tn + fp + fn).to_f
  $stdout.printf("Accuracy %.3f\n", acc)
end

def compute_f1(tp, tn, fp, fn)
  f1 = (2 * tp).to_f / (2 * tp + fp + fn).to_f
  $stdout.printf("F1 %.3f\n", f1)
end

number_samples = 0
number_correct = 0

ambiguous_predictions_skipped = 0

inp.each do |line|
  f = line.split

  id = f[0]

  pred = f[predicted_column]

  if (experimental_value.size > 0)
    expt = experimental_value[id]
  else
    expt = convert_to_class_if_needed1(f[experimental_column], $convert_experimental_to_two_class)
  end

# $stderr.print "id '#{id}', expt '#{expt}', pred '#{pred}'\n"

  unless (expt)
    $stderr.print "No expt value for '#{id}', line #{inp.lineno}\n"
    break
  end
  unless (pred)
    $stderr.print "No pred value for '#{id}', #{line.chomp}'\n"
    break
  end

  pred = convert_to_class_if_needed2(pred, $convert_predicted_to_two_class_lower_bound, $convert_predicted_to_two_class_upper_bound)

  if (! pred)
    ambiguous_predictions_skipped += 1
    next
  end

# $stderr.print "ID '#{id}', expt = #{expt}, pred #{pred}\n"

  if (! $classes.has_key?(expt))
    $classes[expt] = Aclass.new(expt)
    $stderr.print "First instance of class '#{expt}'\n" if (verbose)
  end

  $classes[expt].found_in(pred)

  if (expt == pred)
    number_correct += 1
  end

  number_samples += 1
end

$stdout.print "Found info on #{number_samples} items in #{$classes.size} classes\n" if (verbose)

if (0 == number_samples)
  $stderr.print "No data\n";
  exit 4
end

total_f = 0.0
min_fraction_correct = 1.0

$classes.each do |k, v|
  f = v.report($stdout)

  total_f += f
  if (f < min_fraction_correct)
    min_fraction_correct = f
  end
end

exit 1 if (0 == number_samples)

$stdout.printf("Fraction correct %.3f\n", number_correct.to_f / number_samples.to_f)
$stdout.printf("Average fraction correct %.3f\n", total_f / $classes.size.to_f)
$stdout.printf("Min fraction correct %.3f\n", min_fraction_correct)

if (ambiguous_predictions_skipped > 0)
  applicability = ((number_samples - ambiguous_predictions_skipped).to_f / number_samples.to_f)
  $stdout.printf("#{ambiguous_predictions_skipped} ambiguous predictions skipped, applicability %.3f\n", applicability)
end
exit 0 unless (cl.option_present('pos'))

positive_class = cl.value('pos')

if (! $classes.has_key?(positive_class))
  $stderr.print "Positive class is '#{positive_class}' but no instances it\n"
  $stdout.printf("Sensitivity %.3f PPV %.3f\n", 0.0, 0.0)
  $stdout.printf("Specificity %.3f NPV %.3f\n", 0.0, 0.0)
  exit 4
end

if (2 != $classes.size)
  $stderr.print "Positive class for sensitivity/specificity computation, but #{$classes.size} classes\n"
  exit 5
end

p = false
n = false
$classes.each_key do |k|
  if (positive_class == k)
    p = $classes[k]
    $stderr.print "Positive class is '#{p.class_name}'\n" if (verbose)
  else
    n = $classes[k]
    $stderr.print "Negative class is '#{n.class_name}'\n" if (verbose)
  end
end

raise "No positive class" unless (p)
raise "No negative class" unless (n)

p.compute_ppv(n)
n.compute_npv(p)

if (cl.option_present('sum'))
  s = p.sum_measures_positive_class + n.sum_measures_negative_class
  $stdout.printf("Sum all measures %.3f\n", s)
end

tp,fp = p.compute_tp_fp(n)
tn,fn = n.compute_tn_fn(p)
compute_matthews(tp, tn, fp, fn)
compute_accuracy(tp, tn, fp, fn)
compute_f1(tp, tn, fp, fn)
