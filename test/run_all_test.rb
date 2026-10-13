#!/usr/bin/env ruby

require 'digest'
require 'etc'
require 'fileutils'
require 'pathname'
require 'set'
require 'shellwords'
require 'socket'
require 'tmpdir'

require 'google/protobuf'

require_relative 'lillymol_tests_pb'

# Run unit tests in this directory.

unless ENV.key?('LILLYMOL_HOME')
  $stderr << "Shell variable \${LILLYMOL_HOME} should be defined, taking default\n"
  pn = Pathname.new(File.expand_path(File.path(__FILE__)))
  ENV['LILLYMOL_HOME'] = pn.parent.parent.to_s
end

unless ENV.key?('BUILD_DIR')
  # $stderr << "Shell variable \${BUILD_DIR} should be defined, taking default\n"
  ENV['BUILD_DIR'] = Etc.uname[:sysname]
end

# I think build_dir would work just as well here.
unless ENV.key?('OSTYPE')
# $stderr << "Shell variable \${OSTYPE} should be defined, taking default\n"
  ENV['OSTYPE'] = 'Linux'
end

# Some tests only run inside Lilly.
$inside_lilly = Socket.gethostname.match(/\.lilly.com$/)

require "#{ENV['LILLYMOL_HOME']}/contrib/bin/lib/iwcmdline"

def usage
  $stderr << "Runs regression tests for LillyMol executables\n"
  $stderr << "Each executable is in its own directory.\n"
  $stderr << "In each directory there will be one or more unit tests for that executable\n"
  $stderr << "For those tests that have been converted to use this script, there will be a file 'tests.json'\n"
  $stderr << "in the directory for that executable. That file controls what gets run.\n"

  $stderr << " -rx <rx>                 Only run tests that match regular expression <rx>\n"
  $stderr << " -bindir <dir>            Directory containing executables to test\n"
  $stderr << "                          Default: ${LILLYMOL_HOME}/bin/${BUILD_DIR}\n"
  $stderr << " -copy_fail <dir>         Upon failure, copy the test outcomes to <dir> - dir must already exist\n"

  exit(0)
end

# https://stackoverflow.com/questions/32178998/using-environment-variables-in-a-file-path
# Expand environment variables in `str`, return expanded form
def expand_env(str)
  str.gsub(/\$([a-zA-Z_][a-zA-Z0-9_]*)|\${\g<1>}|%\g<1>%/) do 
    ENV.fetch(Regexp.last_match(1), nil) 
    if ENV.fetch(Regexp.last_match(1), nil) == nil
      "${" + Regexp.last_match(1) + "}"
    else
       ENV.fetch(Regexp.last_match(1), nil)
    end

  end
end

def difftool_options(proto)
  if proto.has_difftool_options?
    return proto.difftool_options
  else
    return ""
  end
end

# Return true if `fname1` and `fname2` have identical contents.
# If their sizes are the same, compare the md5 sums. Otherwise
# shell out to diff.
def files_the_same(proto, fname1, fname2)
  unless File.exist?(fname1)
    $stderr << "Missing file #{fname1}\n"
    return false
  end
  unless File.exist?(fname2)
    $stderr << "Missing file #{fname2}\n"
    return false
  end


  if File.size?(fname1) == File.size?(fname2)
    d1 = Digest::MD5.file(fname1).hexdigest
    d2 = Digest::MD5.file(fname2).hexdigest
    return true if d1 == d2
  end

  rc = if proto.has_difftool?
         $stderr << "Executing #{proto.difftool} #{proto.difftool_options} #{fname1} #{fname2}\n"
         system("#{proto.difftool} #{difftool_options(proto)} #{fname2} #{fname1}")
       else
         system("diff -w #{difftool_options(proto)} #{fname1} #{fname2}")
       end

  if rc 
    return true
  else
    system("wc #{fname1}")
    system("wc #{fname2}")
    return false
  end
end

def to_single_line(lines)
  result = []
  lines.each do |line|
    result << expand_env(line)
  end
  return result.join(' ')
end

class Options
  attr_accessor :lillymol_home, :build_dir, :binary_dir, :default_binary_dir, :tmpdir, :ostype, :verbose, :copy_failures_dir, :only_process_rx

  def initialize(cl)
    @lillymol_home = ENV['LILLYMOL_HOME']
    @build_dir = ENV['BUILD_DIR']
    @default_binary_dir = File.join(@lillymol_home, 'bin', @build_dir)
    @binary_dir = @default_binary_dir

    if cl.option_present('bindir')
      @binary_dir = File.expand_path(cl.value('bindir'))
    end

    @verbose = cl.option_present('v')

    # We do all our work in a tmpdir, which gets removed when we are done.
    # User can specify a directory into which we copy files with diffs.
    if cl.option_present('copy_fail')
      @copy_failures_dir = cl.value('copy_fail')
      FileUtils.mkdir_p(@copy_failures_dir) unless File.directory?(@copy_failures_dir)
    else
      @copy_failures_dir = ""
    end
    $stderr << "copy_failures_dir starts at #{@copy_failures_dir}\n" if verbose
    $stderr << "Executables read from #{@binary_dir}\n" if verbose

    @tmpdir = ""

    @ostype = ENV['OSTYPE'] if ENV.key?('OSTYPE')

    # We can select which tests to run via any number of regex's on the command line.
    @only_process_rx = []

    if cl.option_present('rx')
      set_only_process_regex(cl.values('rx'))
    end
  end

  def set_only_process_regex(rx)
    rx.each do |pattern|
      @only_process_rx << Regexp.new(pattern)
    end
  end

  def perform_test?(test_name)
    return true if @only_process_rx.empty?

    @only_process_rx.each do |rx|
      return true if rx.match(test_name)
    end
    return false
  end
end

def get_test_proto(options, dirname, test_json)
  config = File.join(dirname, test_json)
  unless File.size?(config)
    $stderr << "No config #{config}\n"
    return nil
  end
  contents = File.read(config)
  $stderr << "Contents #{contents}\n" if options.verbose
  begin
    proto = LillymolTests::TestCase.decode_json(contents)
  rescue => e
    $stderr << "get_test_proto: #{dirname} bad json #{e.message}\n"
    $stderr << contents
    return nil
  end

  $stderr << "proto #{proto}\n" if options.verbose
  return proto
end

def get_directory_contents(dirname)
  result = Set.new
  unless File.directory?(dirname)
    $stderr << "#{dirname} not found\n"
    return result
  end

  Dir.open(dirname).each do |fname|
    next if fname[0] == '.'
    result.add(File.join(dirname, fname))
  end

  return result
end

def make_my_tmpdir(tmpdir, dir)
  # $stderr << "make_my_tmpdir #{tmpdir} dir '#{dir}'\n"
  path = File.join(tmpdir, dir)
  if File.directory?(path)
    $stderr << "Temp directory #{path} already exists, is there a duplicate test name?\n"
    return path
  end

  unless Dir.mkdir(path) == 0
    $stderr << "Could not create temporary directory #{path}\n"
    return tmpdir
  end

  return path
end


# If `exe` is the name of a LillyMol executable, return the full path.
# This is designed for things like same_structures, jfilecompare...
def maybe_lillymol_executable(options, exe)
  path = find_lillymol_executable(options, exe)
  return path if path

  # Not a LillyMol executable, return unchanged
  return exe
end

# If `exe` is the name of a LillyMol executable, return the full path.
# If -bindir was specified and lookup there fails, fall back to the historical
# ${LILLYMOL_HOME}/bin/${BUILD_DIR} location.
def find_lillymol_executable(options, exe)
  path = File.join(options.binary_dir, exe)
  return path if File.executable?(path)

  return nil if options.binary_dir == options.default_binary_dir

  path = File.join(options.default_binary_dir, exe)
  return path if File.executable?(path)

  nil
end

# Getting output files is complicated by the fact that some tests
# may have separate directories based on UNAME.
# If there are files to be ignored in `proto` enforce that.
def get_output_files(options, proto, dirname)
  outdir = File.join(dirname, 'out')

  platform_specific = File.join(outdir, options.ostype)

  dc = if File.directory?(platform_specific)
         get_directory_contents(platform_specific)
       else
         get_directory_contents(outdir)
       end

  return dc if proto.ignore_file.empty?

  return dc.select {|fname| ! proto.ignore_file.include?(File.basename(fname))}
end

def output_directory_for_test(options, dirname)
  outdir = File.join(dirname, 'out')
  platform_specific = File.join(outdir, options.ostype)

  return platform_specific if File.directory?(platform_specific)

  outdir
end


# Test JSON is trusted. This intentionally supports existing Ruby interpolation
# forms like #{indir}, #{datadir}, #{ostype}, #{exe}, and #{outdir}.
def interpolate_command_text(text, indir, datadir, ostype, exe, outdir)
  eval("\"" + expand_env(text) + "\"")
end

# Build the command fragment for TestCase.shell_script. The script is found in the
# specific test case directory. The caller is responsible for changing to the
# temporary directory and capturing stdout/stderr.
def command_from_shell_script(options, shell_script, test_case_dir, exe, indir, outdir, datadir, ostype)
  if shell_script.script_name.empty?
    $stderr << "Empty shell_script.script_name in #{test_case_dir}\n"
    return nil
  end

  script = File.join(test_case_dir, shell_script.script_name)
  unless File.file?(script)
    $stderr << "Shell script #{script} missing\n"
    return nil
  end
  unless File.executable?(script)
    $stderr << "Shell script #{script} is not executable\n"
    return nil
  end

  cmd = Shellwords.escape(script)
  cmd << " #{Shellwords.escape(exe)}"
  cmd << " #{Shellwords.escape(indir)}"
  cmd << " #{Shellwords.escape(outdir)}"

  shell_script.args.each do |arg|
    value = interpolate_command_text(arg, indir, datadir, ostype, exe, outdir)
    cmd << " #{value}"
  end

  cmd
end

def command_for_test_case(options, proto, exe, args, mytmp, indir, outdir, datadir, ostype, test_case_dir)
  cmd = ""

  # Preamble entries are shell fragments inserted before the test command.
  # They run in the same shell as the test command so shell variables set here
  # are visible to args/default_command_components.
  if proto.preamble.size > 0
    cmd << "cd #{Shellwords.escape(mytmp)} && "
    cmd << interpolate_command_text(proto.preamble.join("\n") << "\n", indir, datadir, ostype, exe, outdir)
  end

  # Both normal and shell-script tests execute from the per-test temp directory,
  # and both have stdout/stderr captured by the runner.
  cmd << "cd #{Shellwords.escape(mytmp)} && "
  if proto.has_shell_script?
    shell_cmd = command_from_shell_script(options, proto.shell_script, test_case_dir, exe, indir, outdir, datadir, ostype)
    return nil unless shell_cmd
    cmd << shell_cmd
  else
    cmd << Shellwords.escape(exe)
    if proto.default_command_components.size > 0
      cmd << " " << interpolate_command_text(proto.default_command_components.join(' '), indir, datadir, ostype, exe, outdir)
    end

    cmd << " #{args}"
  end

  cmd << " > stdout 2> stderr"
  expand_env(cmd)
end

def command_successful?(cmd, non_zero_rc_expected, verbose)
  $stderr << "Executing #{cmd}\n" if verbose

  return true if system(cmd)
  return true if non_zero_rc_expected

  $stderr << "Warning: #{cmd} failed rc #{$?}\n"
  false
end

def output_files_match?(proto, output_file, mytmp, verbose)
  all_files_same = true

  system("/bin/ls -l #{mytmp}") if verbose

  output_file.each do |correct|
    in_tmp = File.join(mytmp, File.basename(correct))
    $stderr << "Checking diffs #{correct} #{in_tmp}\n" if verbose
    next if files_the_same(proto, correct, in_tmp)

    $stderr << "Diffs btw #{correct}\n      and #{in_tmp}\n"
    system("ls -l #{correct} #{in_tmp}")
    all_files_same = false
  end

  all_files_same
end

def copy_failure_artifacts(options, proto, mytmp, cmd)
  $stderr << "Got failed test where #{options.copy_failures_dir}\n"
  return if options.copy_failures_dir.empty?

  destination = File.join(options.copy_failures_dir, proto.executable, proto.name)
  FileUtils.mkdir_p(destination)

  get_directory_contents(mytmp).each do |fname|
    f = File.basename(fname)
    FileUtils.copy_file(fname, File.join(destination, f))
  end

  File.write(File.join(destination, 'cmd'), "#{cmd}\n")
end

def run_case_proto(options, proto, test_dir, test_name, parent_tmpdir)
  return true if proto.broken_do_not_evaluate

  proto.name = test_name unless proto.has_name?

  mytmp = make_my_tmpdir(parent_tmpdir, proto.name)

  # output files are not expanded.
  # $stderr << "run_case_proto test dir #{test_dir} name #{proto.name}\n"
  output_file = get_output_files(options, proto, File.join(test_dir, proto.name))
  if output_file.empty?
    $stderr << "Nothing to test in #{test_dir} #{proto.name}\n"
    system("/bin/ls -l #{File.join(test_dir, proto.name)}")
    return false
  end
  # $stderr << "Outfiles #{output_file}\n"

  args = to_single_line(proto.args)
  # $stderr << "Args are #{args}\n"

  exe = find_lillymol_executable(options, proto.executable)
  unless exe
    $stderr << "Executable #{File.join(options.binary_dir, proto.executable)} missing"
    if options.binary_dir != options.default_binary_dir
      $stderr << "; fallback #{File.join(options.default_binary_dir, proto.executable)} also missing"
    end
    $stderr << "\n"
    return false
  end

  indir = File.join(test_dir, test_name, 'in')
  if File.directory?(File.join(test_dir, 'data'))
    datadir = File.join(test_dir, 'data')
  else
    datadir = ""
  end
  ostype = options.ostype

  test_case_dir = File.join(test_dir, proto.name)
  outdir = output_directory_for_test(options, test_case_dir)

  args = interpolate_command_text(args, indir, datadir, ostype, exe, outdir)

  cmd = command_for_test_case(options, proto, exe, args, mytmp, indir, outdir, datadir, ostype, test_case_dir)
  return false unless cmd

  successful_execution = command_successful?(cmd, proto.non_zero_rc_expected, options.verbose)
  all_files_same = output_files_match?(proto, output_file, mytmp, options.verbose)

  return true if successful_execution && all_files_same

  copy_failure_artifacts(options, proto, mytmp, cmd)
  return false
end

def run_case(options, dirname, test_name, parent_tmpdir)
  proto = get_test_proto(options, File.join(dirname, test_name), 'test.json')

  if proto
    proto.executable = File.basename(dirname)
    return run_case_proto(options, proto, dirname, test_name, parent_tmpdir)
  end

  run_case = File.join(dirname, test_name, 'run_case.sh')
  return false unless File.executable?(run_case)

  test_dir = File.join(dirname, test_name)
  return system("cd #{test_dir} && ./run_case.sh")
end

def run_tests_in_dir(options, dirname, parent_tmpdir)
  json_fname = File.join(dirname, 'tests.json')
  contents = File.read(json_fname)
  $stderr << "Contents #{contents}\n" if options.verbose
  begin
    proto = LillymolTests::TestCases.decode_json(contents)
  rescue => e
    $stderr << "run_tests_in_dir: #{dirname} bad json #{e.message}\n"
    $stderr << contents
    return 0, 0
  end

  if proto.test.empty?
    $stderr << "No tests in #{proto}\n"
    return 0, 0
  end

  return 0, 0 if proto.only_inside_lilly && ! $inside_lilly

  return 0, 0 if proto.broken_do_not_evaluate

  if proto.has_difftool?
    proto.difftool = maybe_lillymol_executable(options, proto.difftool)
  else
    difftool = 'diff -w'
#   proto.difftool_options = '-w'
  end

  passed = 0
  failed = 0
  proto.test.each do |test_case|
    test_case.executable = proto.executable
    test_case.default_command_components.concat(proto.default_command_components.to_a)
    unless test_case.has_name?
      $stderr << "Skipping test with no name\n"
      $stderr << test_case << "\n"
      next
    end
    test_dir = File.join(dirname, test_case.name)
    unless File.directory?(test_dir)
      $stderr << "NO directory for test #{test_case.name} in #{test_dir}\n"
      failed += 1
      next
    end
    next unless options.perform_test?(test_case.name)

    # Inherit difftool settings from the parent TestCases message.
    if proto.has_difftool? && ! test_case.has_difftool?
      test_case.difftool = proto.difftool
    end
    if proto.has_difftool_options? && ! test_case.has_difftool_options?
      test_case.difftool_options = proto.difftool_options
    end
    # $stderr << test_case << "\n"

    if run_case_proto(options, test_case, dirname, test_case.name, parent_tmpdir)
      passed += 1
    else
      failed += 1
    end
  end

  passfail = if failed == 0
               'PASS'
             else
               'FAIL'
             end
  $stdout << File.basename(dirname) << " #{passed} passed #{failed} failed #{passfail}\n"

  return passed, failed
end

# Return a list of the sub-directories in `dir`
def dirs_in_dir(dir) 
  result = []
  Dir.open(dir).each do |fname|
    next if fname[0] == '.'
    next unless File.directory?(File.join(dir, fname))
    result << fname
  end

  return result
end

# `dir` will be the full path name of the directory containing one or
# more tests, which will be in directories like 'case_*'
# basename(dir) will be the name of the tool being tested.
# return the number of passed and failed tests.
def run_tests(options, dir)
  $stderr << "Begin processing #{dir}'\n" if options.verbose
  tool = File.basename(dir)

  # A temporary directory for all tests of this tool.
  tmpdir = File.join(options.tmpdir, tool)
  Dir.mkdir(tmpdir)

  # If there is a tests.json in the directory, that is how this directory is structured.
  return run_tests_in_dir(options, dir, tmpdir) if File.size?(File.join(dir, 'tests.json'))

  # Test config files in each directory.

  failures = 0
  test_passes = 0
  dirs_in_dir(dir).each do |subdir|
    if run_case(options, dir, subdir, tmpdir)
      $stderr << "#{tool} #{subdir} TEST #{subdir} PASS\n";
      test_passes += 1
    else
      $stderr << "#{tool} #{subdir} TEST #{subdir} FAIL\n";
      failures += 1
    end
  end

  return test_passes, failures
end

# A path name `pathname` has come out of a directory scan.
# There may be zero or more file names in `argv.
# If empty, return true, all tests are to be done
# Check to see if the file name component of `pathname` is in `argv`
def selected_for_processing(pathname, argv)
  return true if argv.empty?
  fname = File.basename(pathname)

  return argv.include?(fname)
end

# Return a sorted list of the subdirectories in `dirname`.
def get_sorted_subdirectories(dirname)
  result = []

  unless File.directory?(dirname)
    $stderr << "#{dirname} not found\n"
    return []
  end

  Dir.open(dirname).each do |fname|
    next if fname[0] == '.'
    pathname = File.join(dirname, fname)
    next unless File.directory?(pathname)
    result << pathname
  end

  return result.sort
end

def main
  cl = IWCmdline.new('-v-copy_fail=s-rx=s-bindir=s-help')

  if cl.option_present('help')
    usage
  end

  if cl.unrecognised_options_encountered
    $stderr << "unrecognised_options_encountered\n"
    usage
  end

  options = Options.new(cl)

  # Strip trailing directory separators
  ARGV.map! { |fname| 
    if fname.end_with?('/')
      fname.chomp('/')
    else
      fname
    end
  }

  ntests = 0
  passing_count = 0
  tools_failing = {}
  failure_count = 0

  pn = Pathname.new(File.expand_path(__FILE__))
  test_dir = pn.parent
  Dir.mktmpdir { |dir| 
    options.tmpdir = dir
    $stderr << "Temp files in #{dir}\n" if options.verbose
    get_sorted_subdirectories(pn.parent).each do |fname|
      next unless selected_for_processing(fname, ARGV)
      next if File.basename(fname) == options.copy_failures_dir
      ntests += 1
      p, f = run_tests(options, fname)
      # $stderr << "p #{p} f #{f}\n"
      passing_count += p
      failure_count += f
      if f > 0
        tools_failing[fname] = 1
      end
    end
  }

# if options.verbose || tools_failing > 0
  total_tests = passing_count + failure_count
  $stderr << "Tested #{ntests} tools with #{total_tests} total tests\n"
  $stderr << tools_failing.size << " tools had failing tests\n"
  if tools_failing.size > 0
    tools_failing.each do |k|
      $stderr << " #{k}\n"
    end
  end
  $stdout << "#{passing_count} tests passed, #{failure_count} tests failed\n"
# end
end

main
