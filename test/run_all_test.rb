!#/usr/bin/env ruby

require 'digest'
require 'etc'
require 'pathname'
require 'set'
require 'tmpdir'
require 'google/protobuf'

require_relative 'test_case_pb'

# Run unit tests in this directory.

unless ENV.key?('LILLYMOL_HOME')
  $stderr << "Shell variable \${LILLYMOL_HOME} should be defined, taking default\n"
  pn = Pathname.new(File.expand_path(File.path(__FILE__)))
  ENV['LILLYMOL_HOME'] = pn.parent.parent.to_s
end

unless ENV.key?('BUILD_DIR')
  $stderr << "Shell variable \${BUILD_DIR} should be defined, taking default\n"
  ENV['BUILD_DIR'] = Etc.uname[:sysname]
end

# I think build_dir would work just as well here.
unless ENV.key?('OSTYPE')
  $stderr << "Shell variable \${OSTYPE} should be defined, taking default\n"
  ENV['OSTYPE'] = 'linux'
end

require "#{ENV['LILLYMOL_HOME']}/contrib/script/ruby/lib/iwcmdline"

# https://stackoverflow.com/questions/32178998/using-environment-variables-in-a-file-path
# Expand environment variables in `str`, return expanded form
def expand_env(str)
  str.gsub(/\$([a-zA-Z_][a-zA-Z0-9_]*)|\${\g<1>}|%\g<1>%/) do 
    ENV.fetch(Regexp.last_match(1), nil) 
  end
end

# Return true if `fname1` and `fname2` have identical contents.
# If their sizes are the same, compare the md5 sums. Otherwise
# shell out to diff.
def files_the_same(fname1, fname2)
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

  rc = system("diff -w #{fname1} #{fname2}")
  if rc == 0
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
  attr_accessor :lillymol_home, :build_dir, :tmpdir, :ostype, :verbose

  def Initialize
    @lillymol_home = ""
    @build_dir = ""
    @tmpdir = ""
    # Does not work, not sure why
    @ostype = ENV['OSTYPE'] if ENV.key?('OSTYPE')
    # We do all our work in a tmpdir, which gets removed when we are done.
    # User can specify a directory into which we copy files with diffs
    @copy_failures_dir = ""
  end

  def set_copy_failures_dir(s)
    @copy_failures_dir = s;
  end
end

def get_test_proto(options, dirname, test_json)
  config = File.join(dirname, test_json)
  unless File.size?(config)
    $stderr << "No config #{config}\n"
    return nil
  end
  contents = IO.read(config)
  $stderr << "Contents #{contents}\n" if options.verbose
  proto = LillyMolTest::TestCase.decode_json(contents)
  $stderr << "proto #{proto}\n" if options.verbose
  return proto
end

def get_directory_contents(dirname)
  result = Set.new
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
  path = File.join(tmpdir, dir)
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

# We have formed a command to execute. If we must capture stdout or stderr
# add redirections to `cmd` and return the new variant.
def maybe_append_stdout_stderr(options, cmd, output_file)

  cmd = "#{cmd} > stdout" if output_file.include?('stdout')

  return "#{cmd} 2> stderr" if output_file.include?('stderr')

  # If verbose, allow stderr to pass through.
  return cmd if options.verbose

  # Not captured, send it to somewhere - could use /dev/null
  return "#{cmd} 2> .stderr"
end

# We have read an array of strings that are either input or output files
# from the proto. Return an array with every member replaced by shell
# expansion, or presence in either `in` or `out` directories.
def maybe_expanded(lines, dirname, in_out) 
  lines.map{ |line| path_for_file(dirname, in_out, expand_env(line)) }
end

# We are testing if output file `fname` is correct. We need to retrieve
# the full path name of it. If `test_dir/out/fname` exists, that is
# the result. Otherwise look OS specific places.
def get_correct_file(options, test_dir, fname)
  result = File.join(test_dir, 'out', fname)
  return result if File.exist?(result)

  return result unless options.ostype

  result = File.join(test_dir, 'out', options.ostype, fname)
  return result if File.exist?(result)

  $stderr << "Missing output file #{test_dir} #{fname}\n"
  return File.join(test_dir, 'out', fname)
end

def run_case(options, dirname, test_name, parent_tmpdir)
  $stderr << "Dirname #{dirname} #{test_name}\n" if options.verbose
  proto = get_test_proto(options, File.join(dirname, test_name), 'test.json')
  return false unless proto

  test_dir = File.join(dirname, test_name)

  mytmp = make_my_tmpdir(parent_tmpdir, test_name)

  input_file = maybe_expanded(proto.input_file, test_dir, 'in')
  # output files are not expanded.
  output_file = proto.output_file
  opts = to_single_line(proto.options)
  # $stderr << "input #{input_file} output #{output_file}\n"

  exe = File.join(options.lillymol_home, 'bin', options.build_dir, proto.executable)
  unless File.executable?(exe)
    $stderr << "Executable #{exe} missing\n"
    return false
  end

  opts = eval("\"" + opts + "\"")
  cmd = "cd #{mytmp} && #{exe} #{opts}"

  cmd = maybe_append_stdout_stderr(options, cmd, output_file)

  cmd = expand_env(cmd)
  $stderr << "Execting #{cmd}\n" if options.verbose

  unless system(cmd)
    $stderr << "Warning: #{cmd} failed\n"
  end

  rc = true

  system("/bin/ls -l #{mytmp}") if options.verbose

  # Check that the files match.
  output_file.each do |fname|
    correct = get_correct_file(options, test_dir, fname)
    in_tmp = File.join(mytmp, File.basename(fname))
    $stderr << "Checking diffs #{correct} #{in_tmp}\n" if options.verbose
    next if files_the_same(correct, in_tmp)
    $stderr << "Diffs btw #{fname}\n      and #{in_tmp}\n"
    rc = false
  end

  rc
end

# `dir` will be the full path name of the directory containing one or
# more tests, which will be in directories like 'case_*'
# basename(dir) will be the name of the tool being tested.
# return the number of failed tests.
def run_tests(options, dir)
  $stderr << "Begin processing #{dir}'\n" if options.verbose
  tool = File.basename(dir)

  # A temporary directory for all tests of this tool.
  tmpdir = File.join(options.tmpdir, tool)
  Dir.mkdir(tmpdir)

  failures = 0
  Dir.open(dir).each do |fname|
    next if fname[0] == '.'
    next unless fname =~ /^case_/

    if run_case(options, dir, fname, tmpdir)
      $stderr << "#{tool} #{fname} TEST #{fname} PASS\n";
    else
      $stderr << "#{tool} #{fname} TEST #{fname} FAIL\n";
      failures += 1
    end
  end

  return failures
end

def main
  cl = IWCmdline.new('-v-copy_fail=dir')

  options = Options.new
  options.lillymol_home = ENV['LILLYMOL_HOME']
  options.build_dir = ENV['BUILD_DIR']
  options.verbose = cl.option_present('v')
  if cl.option_present('copy_fail')
    options.set_copy_failures_dir(cl.value('copy_fail'))
  end
  # Kludge, not sure why things do not work in the Options constructor.
  options.ostype = 'linux' unless options.ostype

  ntests = 0
  tools_failing = 0
  failure_count = 0

  pn = Pathname.new(File.expand_path(__FILE__))
  test_dir = pn.parent
  Dir.mktmpdir { |dir| 
    options.tmpdir = dir
    $stderr << "Temp files in #{dir}\n" if options.verbose
    pn.parent.children.each do |fname|
      next unless File.directory?(fname)
      ntests += 1
      failed = run_tests(options, fname)
      if failed > 0
        failure_count += failed
        tools_failing += 1
      end
    end
  }

  if options.verbose || tools_failing > 0
    $stderr << "Across #{ntests} tools #{tools_failing} failed, total #{failure_count} individual test failures\n"
  end
end

main
