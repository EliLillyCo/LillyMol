!#/usr/bin/env ruby

require 'digest'
require 'etc'
require 'fileutils'
require 'pathname'
require 'set'
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
         # $stderr << "Executing #{proto.difftool} #{fname1} #{fname2}\n"
         system("#{proto.difftool} #{fname2} #{fname1}")
       else
         system("diff -w #{fname1} #{fname2}")
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
  attr_accessor :lillymol_home, :build_dir, :tmpdir, :ostype, :verbose, :copy_failures_dir, :only_process_rx, :same_structures

  def initialize(cl)
    @lillymol_home = ENV['LILLYMOL_HOME']
    @build_dir = ENV['BUILD_DIR']

    @verbose = cl.option_present('v')

    # We do all our work in a tmpdir, which gets removed when we are done.
    # User can specify a directory into which we copy files with diffs.
    if cl.option_present('copy_fail')
      @copy_failures_dir = cl.value('copy_fail')
    else
      @copy_failures_dir = ""
    end
    $stderr << "copy_failures_dir starts at #{@copy_failures_dir}\n" if verbose

    @tmpdir = ""

    @ostype = ENV['OSTYPE'] if ENV.key?('OSTYPE')

    # We can select which tests to run via any number of regex's on the command line.
    @only_process_rx = []

    if cl.option_present('rx')
      set_only_process_regex(cl.values('rx'))
    end

    # Tests using same_structures can use this as a short-cut for getting the path.
    @same_structures = File.join(@lillymol_home, 'bin', @build_dir, 'same_structures')
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
  proto = LillymolTests::TestCase.decode_json(contents)
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

# When scanning a directory holding output file(s), two file names
# are special, stdout and stderr.
# This class scans a directory and records the list of files present,
# while recording whether or not the two special files are present.
class OutDirContents
  attr_reader :has_stdout, :has_stderr, :files
  def initialize(dirname)
    @dir = dirname
    @has_stdout = false
    @has_stderr = false
    @files = []

    Dir.open(dirname).each do |fname|
      next if fname[0] == '.'
      if fname == 'stdout'
        @has_stdout = true
        next
      end
      if fname == 'stderr'
        @has_stderr = true
        next
      end
      @files << File.join(dirname, fname)
    end
  end
  def get_stdout
    return File.join(@dir, 'stdout')
  end
  def get_stderr
    return File.join(@dir, 'stderr')
  end
end

def get_out_directory_contents(dirname)
  result = Set.new
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
    $stderr << "Cound not create temporary directory #{path}\n"
    return tmpdir
  end

  return path
end

# We have read the name of an input or output file from the proto.
# We want to know the full path name.
# Try various things.
# `in_out` will be either `in` or `out` and if we can find the
# file in that directory, we return it.
def path_for_file(dirname, in_out, fname)
  return fname if File.exist?(fname)
  return fname if fname == "stdout"
  return fname if fname == "stderr"

  # First try top level directory.
  maybe_in_dir = File.join(dirname, fname)
  return maybe_in_dir if File.exist?(maybe_in_dir)

  maybe_in_dir = File.join(dirname, in_out, fname)
  return maybe_in_dir if File.exist?(maybe_in_dir)

  $stderr << "Possibly missing input file: dir #{dirname} fname #{fname}\n"
  return fname
end

# We have formed a command to execute. 
# append stdout and stderr captures and return the new variant
def maybe_append_stdout_stderr(options, cmd, output_file)

  return "#{cmd} > stdout 2> stderr"
  # remove the code below...

  cmd = "#{cmd} > stdout" if output_file.include?('stdout')

  return "#{cmd} 2> stderr" if output_file.include?('stderr')

  # If verbose, allow stderr to pass through.
  return cmd if options.verbose

  # Not captured, send it to somewhere - could use /dev/null
  return "#{cmd} 2> .stderr"
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

def discard_ignored_files(proto, files)
  return files if proto.ignore_file.empty?
end

# We have read an array of strings that are either input or output files
# from the proto. Return an array with every member replaced by shell
# expansion, or presence in either `in` or `out` directories.
def maybe_expanded(lines, dirname, in_out) 
  lines.map{ |line| path_for_file(dirname, in_out, expand_env(line)) }
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

  exe = File.join(options.lillymol_home, 'bin', options.build_dir, proto.executable)
  unless File.executable?(exe)
    $stderr << "Executable #{exe} missing\n"
    return false
  end

  indir = File.join(test_dir, test_name, 'in')
  if File.directory?(File.join(test_dir, 'data'))
    datadir = File.join(test_dir, 'data')
  else
    datadir = ""
  end

  args = eval("\"" + args + "\"")

  if (proto.preamble.size > 0)
    cmd = proto.preamble.join("\n") << "\n"
  else
    cmd = ""
  end

  cmd = "cd #{mytmp} && #{exe} "
  if proto.default_command_components.size > 0
    cmd << eval("\"" + proto.default_command_components.join(' ') + "\"")
  end

  cmd << " #{args} > stdout 2> stderr"

  cmd = expand_env(cmd)

  if proto.preamble.size > 0
    cmd = proto.preamble.join("\n") << "\n" << cmd
  end

  $stderr << "Executing #{cmd}\n" if options.verbose

  successful_execution = true
  if system(cmd)
  elsif proto.non_zero_rc_expected
  else
    successful_execution = false
    $stderr << "Warning: #{cmd} failed rc #{$?}\n"
  end

  all_files_same = true

  system("/bin/ls -l #{mytmp}") if options.verbose

  if proto.has_difftool? && proto.difftool == 'same_structures'
    proto.difftool = options.same_structures
  end

  # Check that the files match.
  output_file.each do |correct|
    in_tmp = File.join(mytmp, File.basename(correct))
    $stderr << "Checking diffs #{correct} #{in_tmp}\n" if options.verbose
    next if files_the_same(proto, correct, in_tmp)
    $stderr << "Diffs btw #{correct}\n      and #{in_tmp}\n"
    system("ls -l #{correct} #{in_tmp}")
    all_files_same = false
  end

  return true if successful_execution && all_files_same

  # Failed test, do we need to copy the temporary directory

  $stderr << "Got failed test where #{options.copy_failures_dir}\n"
  return false unless options.copy_failures_dir.length > 0

  # Copy the contents of the result directory to the saved location.
  destination = File.join(options.copy_failures_dir, proto.executable, proto.name)
  # $stderr << "Making #{destination}\n"
  FileUtils.mkdir_p(destination)

  get_directory_contents(mytmp).each do |fname|
    f = File.basename(fname)
    # $stderr << "Copying #{fname} to #{destination}/#{f}\n"
    FileUtils.copy_file(fname, File.join(destination, f))
  end

  File.write(File.join(destination, 'cmd'), "#{cmd}\n")

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
  proto = LillymolTests::TestCases.decode_json(contents)
  if proto.test.empty?
    $stderr << "No tests in #{proto}\n"
    return 0, 0
  end

  return 0, 0 if proto.only_inside_lilly && ! $inside_lilly

  return 0, 0 if proto.broken_do_not_evaluate

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

    # Inherit difftool from parent - I think this is a good idea...
    if (proto.has_difftool? && ! test_case.has_difftool?)
      test_case.difftool = proto.difftool
      if proto.has_difftool_options?
        test_case.difftool = "#{test_case.difftool} #{proto.difftool_options}"
      end
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

# Return a sorted list of the subsirectories in `dirname`.
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
  cl = IWCmdline.new('-v-copy_fail=dir-rx=s-help')

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
