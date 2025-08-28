#!/usr/bin/env ruby

require_relative "lib/iwcmdline"

$expert = false

def usage (rc)
  $stderr.print "Eliminates duplicate structures in 3D\n"
  $stderr.print " -d <dist>      grid spacing - default 1A\n"
  $stderr.print " -t <dist>      radius for leader clustering\n"
  $stderr.print " -mg ... -mg    options passed directly to molecular_grid.sh\n" if $expert
  $stderr.print " -gf ... -gf    options passed directly to molecular_grid.sh\n" if $expert
  $stderr.print " -expert        more options\n" unless ($expert)
  $stderr.print " -v             verbose output\n"
  exit(rc)
end

# IWCmdlineV2  if using v2

cl = IWCmdline.new("-v-expert-keep-d=f-mg=close-gf=close-t=fraction")

$expert = cl.option_present('expert')

if cl.unrecognised_options_encountered()
  $stderr.print "Unrecognised options encountered\n"
  usage(1)
end

verbose = cl.option_present('v')

grid_sep = 1.0

if cl.option_present('d')
  grid_sep = cl.value('d')
end

mg_options = false
if cl.option_present('mg')
  mg_options = cl.value('mg')
end

gf_options = '-P UST:AY -T 1.0'
if cl.option_present('gf')
  gf_options = cl.value('gf')
end

radius = 0.05
if cl.option_present('t')
  radius = cl.value('t')
end

if 0 == ARGV.size
  $stderr.print "Insufficient arguments\n"
  usage(2)
end

fileconv = 'fileconv'
molecular_grid = 'molecular_grid'
grid_fingerprint = 'grid_fingerprint'
nplotnn = 'nplotnn'
gfp_leader = 'gfp_leader'
fetch_sdf = 'fetch_sdf'

tmpdir = '.'

common_structure_reading_options = '-E autocreate -i mdlD -i mdlT -i info -i ignore_bad_m -i mdlquiet'

tmpsdf = "#{tmpdir}/um3d#{Process.pid}"
cmd = "#{fileconv} -n 1 -S #{tmpsdf} #{common_structure_reading_options} -o sdf -o info -o mdlMEND"

ARGV.each do |fname|
  raise "Missing or empty input file '#{fname}'" unless File.size?(fname)
  cmd << " #{fname}"
end

$stderr << "Executing '#{cmd}'\n" if verbose
system(cmd)

tmpsdf << ".sdf"

if ! File.size(tmpsdf)
  $stderr << "'#{cmd}' failed, - number assignment\n"
  usage(1)
end

files_to_be_deleted = Array.new
files_to_be_deleted.push(tmpsdf)

grid_file = "#{tmpdir}/um3d#{Process.pid}.grid"
cmd = "#{molecular_grid} #{common_structure_reading_options} -S #{grid_file}"
if mg_options
  cmd << " #{mg_options}"
else
  cmd << " -d #{grid_sep} #{tmpsdf}"
end

$stderr << "Executing '#{cmd}'\n" if verbose
system(cmd)

grid_file << ".sdf"

if ! File.size?(grid_file)
  $stderr << "Grid generation '#{cmd}' failed\n"
  exit 2
end

files_to_be_deleted.push(grid_file)

cmd = "#{grid_fingerprint} #{common_structure_reading_options} #{gf_options}"

fp = "#{tmpdir}/um3d#{Process.pid}.gfp"
cmd << " #{grid_file} #{tmpsdf} > #{fp}"
$stderr << "Executing '#{cmd}'\n" if verbose
system(cmd)

if ! File.size?(fp)
  $stderr << "Did not generate fingerprints '#{cmd}'\n"
  exit 1
end

files_to_be_deleted.push(fp)

ldr = "#{tmpdir}/um3d#{Process.pid}.ldr"
cmd = "#{gfp_leader} -t #{radius} #{fp} > #{ldr}"
$stderr << "Executing '#{cmd}'\n" if verbose
system(cmd)

if ! File.size?(ldr)
  $stderr << "Did not generate leader file '#{cmd}'\n"
  exit 1
end

files_to_be_deleted.push(ldr)

nnfile = "#{tmpdir}/um3d#{Process.pid}.nn.smi"
cmd = "#{nplotnn} -n 0 #{ldr} > #{nnfile}"
$stderr << "Executing '#{cmd}'\n" if verbose
system(cmd)

if ! File.size?(nnfile)
  $stderr << "Did not generate nn file '#{cmd}'\n"
  exit 1
end

files_to_be_deleted.push(nnfile)

cmd = "#{fetch_sdf} -v -c 2 -f 1 #{nnfile} #{tmpsdf}"
$stderr << "Executing '#{cmd}'\n" if verbose
system(cmd)

if ! cl.option_present('keep')
  files_to_be_deleted.each do |fname|
    File.unlink(fname) if File.exist?(fname)
  end
end
