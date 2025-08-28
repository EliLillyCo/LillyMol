#! /usr/bin/env ruby

require_relative "lib/iwcmdline.rb"

dname = 'D'
idstem = 'ID'

$expert = false

def usage (rc)
  $stderr.print "Generates random descriptors\n"
  $stderr.print " -c <ncol>        number of columns\n";
  $stderr.print " -r <nrow>        number of rows\n";
  $stderr.print " -minval <x>      minimum descriptor value\n";
  $stderr.print " -maxval <x>      minimum descriptor value\n";
  $stderr.print " -R <min,max> ... ranges for each column\n" if $expert
  $stderr.print " -int             produce whole numbers\n";
  $stderr.print " -char <n>        character data, as many as <n> different letters\n";
  $stderr.print " -pchar <f>       approx fraction of columns that should be character data\n";
  $stderr.print " -dname <stem>    stem name for newly created identifiers\n"
  $stderr.print " -s <n>           starting value for identifiers\n"
  $stderr.print " -o <char>        output character separator\n"
  $stderr.print " -F <fname>       take identifiers from file <fname>\n"
  $stderr.print " -fcol <col>      identifiers are in column <col> of <fname>\n"
  $stderr.print " -fapp            append the random data to the -F file\n"
  $stderr.print " -fdes            the -F file is a descriptor file\n"
  $stderr.print " -DNAME <fname>   take descriptor names from <fname>\n"
  $stderr.print " -noid            suppress identifiers in column 1\n";
  $stderr.print " -noh             skip the header record\n"

  exit(rc)
end

cl = IWCmdline.new("-v-c=ipos-r=ipos-minval=f-maxval=f-int-char=ipos-pchar=fraction-noh-s=ipos-o=s-F=sfile-fcol=ipos-fapp-fdes-dname=s-DNAME=sfile-noid-R=s-expert")

$expert = cl.option_present('expert')

if cl.unrecognised_options_encountered()
  $stderr.print "Unrecognised options encountered\n"
  usage(1)
end

verbose = cl.option_count('v')

ncols = 10;
if (cl.option_present('c'))
  ncols = cl.value('c')
end

if cl.option_present('R')
  ncols = cl.option_count('R')
end

nrows = 10
if (cl.option_present('r'))
  nrows = cl.value('r')
end

$minval = 0.0;
if (cl.option_present('minval'))
  $minval = cl.value('minval')
end

$maxval = 1.0;
if (cl.option_present('maxval'))
  $maxval = cl.value('maxval')
end

raise "Minval #{$minval} must be less than maxval #{$maxval}" unless ($minval < $maxval)

$alphabet = 'abcdefghijklmnopqrstuvwxyz'
$chars_in_string_data = false

all_columns_string_data = false
number_character_columns_to_create = 0

if (cl.option_present('char'))
  $chars_in_string_data = cl.value('char')
  $alphabet = $alphabet[0..$chars_in_string_data]
  $chars_in_string_data += 1
  all_columns_string_data = true
end

if (cl.option_present('pchar'))
  t = cl.value('pchar')
  number_character_columns_to_create = (ncols * t + 0.499).to_i
  all_columns_string_data = false
  if (verbose)
    $stderr.print "Will create #{number_character_columns_to_create} character columns\n"
  end
end

$whole_numbers = cl.option_present('int')

idstart = 0

if (cl.option_present('s'))
  idstart = cl.value('s')
end

suppress_identifiers = cl.option_present('noid')

identifiers_to_use = Array.new

existing_data = Hash.new

fapp = cl.option_present('fapp')

descriptor_header_record = ""

if (cl.option_present('F'))
  fname = cl.value('F')
  inp = File.open(fname, mode='r')
  raise "Cannot open '#{fname}'" unless inp

  c = 0
  if (cl.option_present('fcol'))
    c = cl.value('fcol')
    c = c - 1
  end

  if cl.option_present('fdes')
    descriptor_header_record = inp.gets
  end

  inp.each do |line|
    f = line.chomp.split
    id = f[c]
    identifiers_to_use.push(id)
    existing_data[id] = line.chomp if fapp
  end
  nrows = identifiers_to_use.size

  $stderr << "Read #{nrows} identifiers from '#{fname}'\n" if verbose
else
  (0...nrows).each do |i|
    identifiers_to_use.push("#{idstem}#{i+idstart}")
  end
end

# Need to increment the max if using whole numbers

if ($whole_numbers)
  if ($minval.to_i != $minval || $maxval.to_i != $maxval)
    $stderr.print "When using whole numbers, minval and maxval must be whole numbers\n"
    usage(3)
  end
  $maxval += 1
end

output_separator = ' '

if (cl.option_present('o'))
  output_separator = cl.value('o')
  if "tab" == output_separator
    output_separator = "\t"
  elsif "vbar" == output_separator
    output_separator = "|"
  end
end

if cl.option_present('dname')
  dname = cl.value('dname')
end

descriptor_names_from_elsewhere = false

if cl.option_present('DNAME')
  fname = cl.value('DNAME')
  if ! File.size?(fname)
    $stderr << "Missing or empty descriptor name file '#{fname}'\n"
    exit(1)
  end

  inp = File.open(fname, mode='r')
  line = inp.gets
  descriptor_names_from_elsewhere = line.chomp.split
  descriptor_names_from_elsewhere.shift
  ncols = descriptor_names_from_elsewhere.size
end

srand(Process.pid)

if (cl.option_present('h') || cl.option_present('help'))
  usage(1)
end

class Column
  def initialize(is_char)
    @_is_char = is_char
    @_minval = $minval
    @_maxval = $maxval
  end

  def set_minmax(x,y)
    @_minval = x
    @_maxval = y
  end

  def next_value
#   $stderr.print "next_value between #{$minval} and #{$maxval} #{@_is_char}\n"
    if (@_is_char)
      t = rand($chars_in_string_data)
      return $alphabet[t..t]
    elsif ($whole_numbers)
      return @_minval.to_i + rand(@_maxval.to_i - @_minval.to_i)
    else
      return sprintf("%.3f", @_minval + (@_maxval - @_minval)*rand)
    end
  end
end

col = Array.new

character_columns_created = 0

(0...ncols).each do |i|
  if (all_columns_string_data)
    c = Column.new(true)
    col.push(c)
  elsif (character_columns_created < number_character_columns_to_create && rand < 0.5)
    c = Column.new(true)
    col.push(c)
    character_columns_created += 1
  else
    c = Column.new(false)
    col.push(c)
  end
end

# Are there individually specified ranges

if cl.option_present('R')
  if cl.option_count('R') != ncols
    $stderr << "Specified #{ncols} columns, but only #{cl.option_count('R')} range, cannot continue\n";
    exit(1)
  end
  (0...ncols).each do |i|
    s = cl.value('R', i)
    f = s.split(',')
    x = false
    y = false
    if $whole_numbers
      x = f[0].to_i
      y = f[1].to_i
    else
      x = f[0].to_f
      y = f[1].to_f
    end
    if x > y
      $stderr << "Invalid range '#{s}'\n"
      exit 1
    end
    col[i].set_minmax(x, y)
#   $stderr << "Examining #{f}\n"
  end
end

if (! cl.option_present('noh'))
  o = output_separator.dup
  if suppress_identifiers
    o = ""
  else
    $stdout.print "Name";
  end
  (0...ncols).each do |i|
    if descriptor_names_from_elsewhere
      $stdout.print "#{o}#{descriptor_names_from_elsewhere[i]}"
    else
      $stdout.print "#{o}#{dname}#{i}"
    end
    o = output_separator.dup
  end

  $stdout.print "\n"
end

if descriptor_header_record.length > 0
  $stdout << descriptor_header_record.chomp
  (0...ncols).each do |i|
    $stdout.print "#{output_separator}#{dname}#{i}"
  end
  $stdout.print "\n"
end

identifiers_to_use.each do |i|
# $stdout.print "#{idstem}#{i+idstart}"
  o = output_separator.dup
  if fapp
    $stdout.print "#{existing_data[i]}"
  elsif suppress_identifiers
    o = ""
  else
    $stdout.print i
  end

  col.each do |c|
    $stdout << "#{o}#{c.next_value}"
    o = output_separator.dup
  end
  $stdout.print "\n"
end
