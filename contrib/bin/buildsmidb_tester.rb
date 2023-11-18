#!/usr/bin/env ruby

ianhome = ENV['LILLYMOL_HOME']

require 'securerandom'
require "#{ianhome}/contrib/bin/lib/iwcmdline.rb"

$expert = false

def usage (rc)
  $stderr << "Tester for buildsmidb and in_database\n"
  $stderr << " -build <fname>    executable for buildsmidb\n"
  $stderr << " -look <fname>     executable for in_database\n"
  $stderr << " -f <n>            stop after <n> test failures\n"
  $stderr << " -keep             keep temporary files\n" if $expert
  $stderr << " -tmpdir <dir>     directory for temporary files (databases)\n"
  $stderr << " -expert           more options\n" unless ($expert)
  $stderr << " -v                verbose output\n"
  exit(rc)
end

cl = IWCmdline.new("-v-expert-keep-look=xfile-build=xfile-tmpdir=dir-f=ipos")

$expert = cl.option_present('expert')

if cl.unrecognised_options_encountered()
  $stderr << "Unrecognised options encountered\n"
  usage(1)
end

verbose = cl.option_present('v')

$keep_temporary_files = cl.option_present('keep')

if ! cl.option_present('build') || ! cl.option_present('look')
  $stderr << "Must specify path to both buildsmidb_bdb (-build) and in_database (-look)\n";
  usage(1)
end

buildsmidb = cl.value('build')
in_database = cl.value('look')

tmpdir = '.'

if cl.option_present('tmpdir')
  tmpdir = cl.value('tmpdir')
end

if ARGV.empty?
  $stderr << "Insufficient arguments\n"
  usage(2)
end

smiles = ARGV[0]
if ! File.size?(smiles)
  $stderr << "Missing or empty smiles file '#{smiles}'\n"
  exit(1)
end

def run_test(buildsmidb, in_database, tmpdir, load_options, lookup_options, smiles, verbose)

  dbname = File.join(tmpdir, "#{SecureRandom.hex}.bdb")

  cmd = "#{buildsmidb} #{load_options} -d #{dbname} #{smiles}"
  $stderr << "Executing '#{cmd}'\n" if verbose
  system(cmd)

  if ! File.size?(dbname)
    $stderr << "'#{cmd}' failed to create '#{dbname}'\n"
    return false
  end

  u = SecureRandom.hex
  ufile = File.join(tmpdir, u)

  cmd = "#{in_database} -U #{ufile} -d #{dbname} #{lookup_options} #{smiles}"
  $stderr << "Executing '#{cmd}'\n" if verbose
  system(cmd)

  ufile << ".smi"

  if ! File.exist?(ufile)
    $stderr << "'#{cmd}' failed to create '#{ufile}'\n"
    return false
  end

  if 0 == File.size(ufile)
    File.unlink(dbname)
    File.unlink(ufile)
    return true
  end

  $stderr << "Build/Lookup failed. Build options '#{load_options}' lookup options '#{lookup_options}' db '#{dbname}' -U file #{ufile}\n"
  return false
end

max_failures_allowed = 4

if cl.option_present('f')
  max_failures_allowed = cl.value('f')
end

failures = 0
if ! run_test(buildsmidb, in_database, tmpdir, '-g all -E autocreate', '-g all -E autocreate', smiles, verbose)
  failures += 1 
  exit 1 if failures >= max_failures_allowed
end

if ! run_test(buildsmidb, in_database, tmpdir, '-g all -b -E autocreate', '-g all -b -E autocreate', smiles, verbose)
  failures += 1 
  exit 1 if failures >= max_failures_allowed
end

if ! run_test(buildsmidb, in_database, tmpdir, '-g all -H allg -E autocreate', '-g all -H allg -E autocreate', smiles, verbose)
  failures += 1 
  exit 1 if failures >= max_failures_allowed
end

if ! run_test(buildsmidb, in_database, tmpdir, '-g all -H allg -E autocreate', '-g all -j -E autocreate', smiles, verbose)
  failures += 1 
  exit 1 if failures >= max_failures_allowed
end

if 0 == failures
  $stderr << "All tests successful\n"
  exit(0)
end

$stderr << "Encountered #{failures} test failures\n"
exit(1)
