#!/usr/bin/env ruby

lillymol_home = ENV['LILLYMOL_HOME']
require "#{lillymol_home}/contrib/bin/lib/iwcmdline.rb"

$expert = false

def usage (rc)
  $stderr.print "Generates peptides\n"
  $stderr.print " -cyclic        join start to end\n"
  $stderr.print " -l             specify components via single letter abbreviations\n"
  $stderr.print " -a             specify components via 3 letter abbreviations\n"
  $stderr.print " -cyclic        close the ring\n"
  $stderr.print " -I             remove isotopes from product molecules\n"
  $stderr.print " -expert        more options\n" unless ($expert)
  $stderr.print " -v             verbose output\n"
  exit(rc)
end

cl = IWCmdline.new("-v-expert-keep-cyclic-cycle-l-a-I")

$expert = cl.option_present('expert')

if cl.unrecognised_options_encountered()
  $stderr.print "Unrecognised options encountered\n"
  usage(1)
end

verbose = cl.option_present('v')

if 0 == ARGV.size
  $stderr.print "Insufficient arguments\n"
  usage(2)
end

make_these_molecules = 'make_these_molecules'
molecular_transformations = 'molecular_transformations'

datadir = "#{lillymol_home}/contrib/data/Peptides"

if cl.option_present('lib')
  datadir = cl.value('lib')
end

peptides = "#{datadir}/peptides_with_isotopic_cysteines.smi"
cysteine = File.join(datadir, 'Cysteine.smi')

if ! File.size?(peptides)
  $stderr << "Where is my peptides file '#{peptides}', cannot continue\n"
  exit(1)
end

if ! File.size?(cysteine)
  $stderr << "Where is my Cysteine file '#{cysteine}', cannot continue\n"
  exit(1)
end

join_rxn = "#{datadir}/peptide_join.rxn"

if ! File.size?(join_rxn)
  $stderr << "Where is my joining reaction '#{join_rxn}', cannot continue\n"
  exit 1
end

col = 2

if cl.option_present('l')
  col = 2
elsif cl.option_present('a')
  col = 3
end

tmpdir = '.'

files_to_be_deleted = Array.new

iwcut = 'iwcut'
tmpp = "#{tmpdir}/peptidex#{Process.pid}.smi"
cmd = "#{iwcut} -f 1,#{col} #{peptides} > #{tmpp}"
system(cmd)
if ! File.size?(tmpp)
  $stderr << "Could not create temporary smiles peptide file '#{cmd}'\n"
  exit 1
end
files_to_be_deleted.push(tmpp)
peptides = "#{tmpp}"

inp = File.open(peptides, mode='r')
raise "Cannot open '#{peptides}'" unless inp

peptide_name_to_smi = Hash.new
peptide_letter_to_smi = Hash.new
peptide_abbrev_to_smi = Hash.new

inp.each do |line|
  f = line.chomp.split

  peptide_name_to_smi[f[1]] = f[0]
  peptide_letter_to_smi[f[2]] = f[0]
  peptide_abbrev_to_smi[f[3]] = f[0]
end

inp.close

fname = ARGV.shift

if ! File.size?(fname)
  $stderr << "Missing or empty input file '#{fname}'\n"
  exit 1
end

zdata = File.readlines(fname)
if 0 == zdata.size
  $stderr << "Cannot read data from '#{fname}'\n"
  exit 1
end

ends_with_digit = Regexp.new('(\d)$')

number_components = 0

ring_openings_found = Array.new(9,0)

zdata.each do |line|
  f = line.chomp.split
  if 0 == number_components
    number_components = f.size
  elsif f.size != number_components
    $stderr << "Size mismatch, expected #{number_components}, so '#{line.chomp}' invalid\n"
    exit 1
  end

  f.each do |p|
    m = ends_with_digit.match(p)
    if m
      r = m[1].to_i
      ring_openings_found[r] += 1
    end
  end
end

rings = Array.new

ring_openings_found.each_with_index do |r, i|
  if r > 0
    $stderr << "Ring open/close #{i} #{r}\n"
    if 0 != r % 2
      $stderr << "Mismatched ring opening/closing #{i}\n"
      exit 1
    end
    rings.push(i)
  end
end

cmd = "#{make_these_molecules} -M #{fname} -W -"
if cl.option_present('I')
  cmd << ' -I'
end
(number_components-1).times do |i|
  cmd << " -R #{join_rxn}"
end

cmd << ' -v' if verbose

number_components.times do |i|
  cmd << " #{peptides}"
end

if cl.option_present('cyclic') || cl.option_present('cycle')
  ring_closing_reaction = "#{datadir}/ring_close.rxn"
  if ! File.size?(ring_closing_reaction)
    $stderr << "Where is my ring closing reaction '#{ring_closing_reaction}', cannot continue\n"
    exit 1
  end

  cmd += " | #{molecular_transformations} -z i -z w -R #{ring_closing_reaction} -i smi -"
end

if rings.size > 0
  isotopic_ring_closing = "#{ianhome}/lib/isotopic_ring_closing"

  cmd += "|#{molecular_transformations} -i smi"
  rings.each do |r|
    rxn = File.join(isotopic_ring_closing, "#{r}.rxn")
    if ! File.size?(rxn)
      $stderr << "Where is '#{rxn}', cannot continue\n"
      exit 1
    end

    cmd += " -R #{rxn}"
  end

  cmd += " -"
end

$stderr << "Executing '#{cmd}'\n" if verbose
system(cmd)

unless cl.option_present('keep')
  files_to_be_deleted.each do |fname|
    File.unlink(fname) if File.exist?(fname)
  end
end
