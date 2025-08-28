#!/usr/bin/env ruby

# Implementation of Lilly Medchem rules using latest LillyMol executables.

ianhome = File.dirname($0)      # location for supporting files

ianhome = ENV['LILLYMOL_HOME']
unless ianhome
  ianhome = File.dirname(File.dirname(File.dirname __FILE__))
end

require "#{ianhome}/contrib/bin/lib/iwcmdline.rb"

$expert = false

$default_lower_atom_count_cutoff = 7
$default_soft_upper_atom_count_cutoff = 25
$default_hard_upper_atom_count_cutoff = 40

def usage (rc)
  $stderr.print "Runs the Lilly medchem rules\n"
  $stderr.print " -c <n>         lower atom count cutoff (default #{$default_lower_atom_count_cutoff})\n" if $expert
  $stderr.print " -Cs <n>        soft upper atom count cuttof (default #{$default_soft_upper_atom_count_cutoff})\n" if $expert
  $stderr.print " -Ch <n>        hard upper atom count cuttof (default #{$default_hard_upper_atom_count_cutoff})\n" if $expert
  $stderr.print " -smarts <s>    optional smarts to reject\n" if $expert
  $stderr.print " -rej <q>       optional query file to reject\n" if $expert
  $stderr.print " -relaxed       relaxed rules: 7-50 heavy atoms, 160 demerit cutoff\n"
  $stderr.print " -nodemerit     hard rejections only, do not apply any demerits\n" if $expert
  $stderr.print " -S <fname>     write output to <fname> rather than stdout\n" if $expert
  $stderr.print " -B <stem>      output name stem for rejected molecules\n" if $expert
  $stderr.print " -log <stem>    name stem for log files\n" if $expert
  $stderr.print " -tp...-tp      options passed directly to mc_first_pass\n" if $expert
  $stderr.print " -iwd...-iwd    options passed directly to iwdemerit\n" if $expert
  $stderr.print " -odm <name>    omit demerit with file name <name>\n" if $expert
  $stderr.print " -edm <fname>   extra demerits to be applied - query file format only\n" if $expert
  $stderr.print " -dmrt <fname>  extra demerits to be applied\n" if $expert
  $stderr.print " -dcf  <fname>  demerit control file (the -C option to iwdemerit)\n" if $expert
  $stderr.print " -q <dir>       directory for queries\n" if $expert
  $stderr.print " -okiso         allow isotopic atoms to pass through\n";
  $stderr.print " -symm <bonds>  discard symmetric molecules where two symmetric atoms > <bonds> apart\n" if $expert
  $stderr.print " -noapdm        do not append demerit reasons\n"
  $stderr.print " -i <type>      input type\n" if $expert
  $stderr.print " -expert        more options\n" unless $expert;
  $stderr.print " -v             verbose output\n"
  exit(rc)
end

cl = IWCmdline.new("-v-noapdm-i=s-expert-b=fraction-B=s-q=dir-log=s-tp=close-iwd=close-smarts=s-rej=s-c=ipos-Cs=ipos-Ch=ipos-okiso-odm=s-edm=sfile-relaxed-nodemerit-S=s-dcf=sfile-nobadfiles-symm=ipos")

if cl.unrecognised_options_encountered()
  $stderr.print "Unrecognised options encountered\n"
  usage(1)
end

verbose = cl.option_present('v')

if cl.option_present('expert')
  $expert = true
end

if ARGV.empty?
  $stderr.print "Insufficient arguments\n"
  usage(2)
end

append_demerit_reason = ! cl.option_present('noapdm')   # convoluted inverse logic

stop_afer_completing_step = 3    # no longer optional

input_type = ""
if cl.option_present('i')
  input_type = '-i ' << cl.value('i')
end

# The ring bond ratio used by mc_first_pass

ring_bond_ratio = -1
if cl.option_present('b')
  ring_bond_ratio = cl.value('b')
end

mc_first_pass_options = ""
if cl.option_present('tp')
  mc_first_pass_options = cl.value('tp')
end

unless cl.option_present('okiso')
  mc_first_pass_options << ' -I 0'
end

mc_first_pass_options << ' -A I -A ipp'

iwdemerit      = "#{ianhome}/bin/Linux/iwdemerit"
mc_first_pass  = "#{ianhome}/bin/Linux/tp_first_pass"
tsubstructure  = "#{ianhome}/bin/Linux/tsubstructure"

query_dir = "#{ianhome}/data/LillyMedchemRules"

query_dir = cl.value('q') if cl.option_present('q')

# The default query files

query_file = Array.new

query_file[1] = 'reject1'
query_file[2] = 'reject2'
query_file[3] = 'demerits'

# The default stem of the rejected files

bad_stem = false

if cl.option_present('B')
  bad_stem = cl.value('B')
elsif cl.option_present('nobadfiles')
  true
else
  bad_stem = 'bad'
end

# The default stem for log files

logfilestem = 'ok';

if cl.option_present('log')
  logfilestem = cl.value('log')
end

optional_queries = ""

cl.values('rej').each do |q|
  optional_queries << " -q #{q}"
end

cl.values('smarts').each do |s|
  optional_queries << " -s '#{s}'"
end

extra_iwdemerit_options = ""

if cl.option_present('relaxed')
  $default_soft_upper_atom_count_cutoff = 26
  $default_hard_upper_atom_count_cutoff = 50
  extra_iwdemerit_options << " -f 160"
end

if cl.option_present('nodemerit')
  extra_iwdemerit_options << " -r"
  $default_soft_upper_atom_count_cutoff = $default_hard_upper_atom_count_cutoff - 1
end

if cl.option_present('iwd')
  extra_iwdemerit_options << " " << cl.value('iwd')
end

if cl.option_present('symm')
  extra_iwdemerit_options << " -s " << cl.value('symm').to_s
end

charge_assigner = "#{ianhome}/data/queries/charges/queries"

unless FileTest.size?(charge_assigner)
  $stderr << "Charge assigner not available, skipping\n"
else
  extra_iwdemerit_options << " -N F:#{charge_assigner}"
end

lower_atom_count_cutoff = $default_lower_atom_count_cutoff
if cl.option_present('c')
  lower_atom_count_cutoff = cl.value('c')
end

soft_upper_atom_count_cutoff = $default_soft_upper_atom_count_cutoff;
if cl.option_present('Cs')
  soft_upper_atom_count_cutoff = cl.value('Cs')
end

hard_upper_atom_count_cutoff = $default_hard_upper_atom_count_cutoff;
if cl.option_present('Ch')
  hard_upper_atom_count_cutoff = cl.value('Ch')
  unless cl.option_present('Cs')
    soft_upper_atom_count_cutoff = hard_upper_atom_count_cutoff - 1
  end
end

if hard_upper_atom_count_cutoff < soft_upper_atom_count_cutoff
  hard_upper_atom_count_cutoff = soft_upper_atom_count_cutoff + 1
end

unless FileTest.directory?(query_dir)
  $stderr.print "Cannot continue, query dir '#{query_dir}' invalid\n"
  exit(3)
end

raise "Query file '#{query_dir}/#{query_file[1]}' missing or inaccessible" unless (FileTest.size?("#{query_dir}/#{query_file[1]}") && FileTest.readable?("#{query_dir}/#{query_file[1]}"))
raise "Query file '#{query_dir}/#{query_file[2]}' missing or inaccessible" unless (FileTest.size?("#{query_dir}/#{query_file[2]}") && FileTest.readable?("#{query_dir}/#{query_file[2]}"))
raise "Query file '#{query_dir}/#{query_file[3]}' missing or inaccessible" unless (FileTest.size?("#{query_dir}/#{query_file[3]}") && FileTest.readable?("#{query_dir}/#{query_file[3]}"))

$stderr.print "Queries from '#{query_dir}'\n" if verbose

files_to_be_deleted = Array.new

query_file3 = "#{query_dir}/#{query_file[3]}"

# I want to have the odm option behave as either a regular expression
# or as exact matches.
# Note that I don't check against multiple RX= directives

class Odm
  def initialize(o)
    @rx = false
    @hash = Hash.new

    o.each do |d|
      if 'RX=' == d[0,3]
        @rx = Regexp.new(d[3, d.size - 3], Regexp::IGNORECASE)
      else
        @hash[d] = true
      end
    end
  end

  def match(d)
    if @hash.has_key?(d)
      return true
    end

    if @rx
      return @rx.match(d)
    end

    return false
  end
end

if cl.option_present('odm')

  tmpdir = "."
  if cl.option_present('tmpdir')
    tmpdir = cl.value('tmpdir')
  end

  old_demerits = File.open("#{query_dir}/#{query_file[3]}", mode='r')
  raise "Cannot open original demerit file '#{query_dir}/#{query_file[3]}'" unless (old_demerits)

  temporary_demerit_file = "#{tmpdir}/demerits" << Process.pid.to_s
  new_demerits = File.open(temporary_demerit_file, mode='w')
  raise "Cannot open temporary demerit file '#{temporary_demerit_file}'" unless (new_demerits)

  odm = Odm.new(cl.values('odm'))

  demerit_rx = Regexp.new("^(\\S+)\.qry")

  items_discarded = 0

  old_demerits.each do |line|
    m = demerit_rx.match(line)

    next unless (m)

    stem = m[1]

    if odm.match(stem)
      items_discarded += 1
    else
      new_demerits << "#{query_dir}/#{stem}.qry\n"
    end
  end

  if 0 == items_discarded
    $stderr.print "Warning, no demerits discarded\n"
  elsif verbose
    $stderr.print "Discarded #{items_discarded} demerits\n"
  end

  new_demerits.close

  query_file3 = "#{temporary_demerit_file}"

  files_to_be_deleted.push(query_file3)
end

iwdemerit_optional_control_file = false

if cl.option_present('dcf')
  iwdemerit_optional_control_file = cl.value('dcf')
end

if cl.option_present('edm')
  additional_demerits = cl.value('edm')
end

cmd = "#{mc_first_pass} ";

cmd << " -b #{ring_bond_ratio}" if ring_bond_ratio >= 0.0

cmd << " #{mc_first_pass_options}" if mc_first_pass_options.length > 0

cmd << " #{input_type} " if input_type.length > 0

cmd << " -c #{lower_atom_count_cutoff} -C #{hard_upper_atom_count_cutoff} -E autocreate -o smi -V -g all -g ltltr -i ICTE "
cmd << "-L #{bad_stem}0 -K TP1 " if bad_stem
cmd << "-a -S - #{ARGV.join(' ')} 2> #{logfilestem}0.log "

if stop_afer_completing_step >= 1
  cmd << "| #{tsubstructure} -E autocreate -b -u -i smi -o smi -A D "
  cmd << "-m #{bad_stem}1 -m QDT " if bad_stem
  cmd << "-n - -q F:#{query_dir}/#{query_file[1]} "

  cmd << optional_queries if optional_queries.length > 0

  cmd << " - 2> #{logfilestem}1.log ";

  if stop_afer_completing_step >= 2
    cmd << "| #{tsubstructure} -A D -E autocreate -b -u -i smi -o smi "
    cmd << "-m #{bad_stem}2 -m QDT " if bad_stem
    cmd << "-n - -q F:#{query_dir}/#{query_file[2]} - 2> #{logfilestem}2.log ";
    if stop_afer_completing_step >= 3
      cmd << " | #{iwdemerit} -x #{extra_iwdemerit_options} -E autocreate -A D -i smi -o smi -q F:#{query_file3} "
      cmd << "-R #{bad_stem}3 " if bad_stem
      cmd << "-G - -c smax=#{soft_upper_atom_count_cutoff} -c hmax=#{hard_upper_atom_count_cutoff} "
      cmd << "-q F:#{additional_demerits} " if additional_demerits
      cmd << "-C #{iwdemerit_optional_control_file} " if iwdemerit_optional_control_file
      cmd << "-t " if append_demerit_reason
      cmd << "- 2> #{logfilestem}3.log "
    end
  end
end

if cl.option_present('S')
  s = cl.value('S')
  cmd << " > #{s}"
end

$stderr.print "Command is '#{cmd}'\n" if verbose

system(cmd)

files_to_be_deleted.each do |f|
  File.unlink(f) if FileTest.exists?(f)
end
