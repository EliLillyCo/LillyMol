#!/usr/bin/env ruby

# Descritpor generator
# Works by examining all the descriptors specified on the
# command line and creating a Makefile that generates the
# needed descriptors.
# Some descriptors not suppored.

require 'fileutils'
require 'set'

c3tk_home = ENV['C3TK_HOME']
raise 'C3TK_HOME not defined' unless c3tk_home

require "#{c3tk_home}/bin/ruby/lib/iwcmdline"

def usage(descriptors)
msg = <<-END
Descriptor generator.
Descriptors are specified on the command line and a Makefile describing the dependencies
is constructed in a temporary directory. The 'make' command then makes all the descriptors
and uses concat_files to generate a concatenated set to stdout.
Descriptors can be specified either separately
-abr -ap -cmi -jurs     or via the -c option    -c abr,ap,cmi,jurs
Uppercase variants can be used to pass optional arguments to an underlying tool
-abr -W -O none -O alogp -W -jurs...
The following options are recognised (* means 3D)
END
  $stderr << msg
  descriptors.each do |k, v|
    $stderr << " -#{k}"
    $stderr << "*" if v.threed
    $stderr << "\n"
  end

  $stderr << " -o <fname>       output file containing the descriptors\n"
  $stderr << "                  or omit this and output will be to stdout.\n"
  $stderr << "                  But note that if output is to stdout, `make` will run silently\n"
  $stderr << "                  That is because `make` echo's commands to stdout\n"
  $stderr << " -j <n>           degree of parellelism 'make -j <n>\n"
  $stderr << " -speed           do nothing, but report the relative speed of each descriptor computation\n"
  $stderr << " -keep            do NOT delete the temporary directory when done, keep it\n"
  $stderr << " -tmpdir <dir>    specify your own temporary directory - will not be deleted\n"
  $stderr << " -v               verbose output\n"

  exit(1)
end

class Descriptor
  attr_accessor :name
  attr_accessor :threed
  attr_accessor :programme
  attr_accessor :extra
  attr_accessor :speed
  def initialize(programme, threed, speed)
    # The name of this descriptor - what the user enters as a command line option.
    @name = ""

    # The tool that implements this desctiptor
    @programme = programme

    # True if this is a 3D descriptor
    @threed = threed

    # Any extra arguments that might need to be passed to the executable.
    @extra = ""

    # A relative measure of speed, with smaller numbers being better.
    # In reality it is the time taken to compute 50k molecules.
    @speed = speed
  end
end

# Should we also remove isotopes? 
# Upper atom count cutoff?
def fileconv_options
  return '-f lod -E autocreate -e -V -g all -O def -c 2'
end

# The flush option is not needed if there are >1 3D descriptors being computed.
def corina_options
  return '-r 5 -x -Y flush'
end

def special_handling_cmi(fname)
  cmd = "fileconv.sh -f lod -E autocreate -e -V -g all -O def -C 99 -S - #{fname} | cmi.sh -"
  $stderr << "executing #{cmd}\n"
  system(cmd)
  
  return true
end

def special_handling_marvin(fname)
  cmd = "fileconv.sh -f lod -E autocreate -e -V -g all -O def -C 99 -S - #{fname} | marvin.sh -"
  $stderr << "executing #{cmd}\n"
  system(cmd)
  
  return true
end

# Request for just a single descriptor. We do not need to create a Makefile
# just issue a pipelined command.
def do_single_descriptor(d, fname, cl)
  return special_handling_cmi(fname) if d.name == 'cmi'
  return special_handling_marvin(fname) if d.name == 'marvin'
  if cl.option_present('nostd')
    cmd = ""
    input_for_next_stage = fname
  else
    cmd = "fileconv.sh #{fileconv_options} -S - #{fname}|"
    input_for_next_stage = '-i smi -'
  end

  cmd << "rcorina.sh #{corina_options} -u #{input_for_next_stage}|" if d.threed

  cmd << "#{d.programme} -i smi -"
  cmd << ' > ' << cl.value('o') if cl.option_present('o')

  $stderr << cmd << "\n" if cl.option_present('v')
  return(system(cmd))
end

# Creates a Makefile in `tmpdir` which forces descriptors in `to_compute` to be
# computed and sets up concat_files to generate a single output.

def create_makefile(to_compute, tmpdir, cl)
  # Spit into 2d and 3d descriptors
  d2 = to_compute.filter { |d|! d.threed }
  d3 = to_compute.filter { |d|  d.threed }

  $stderr << "Have #{d2.size} 2D and #{d3.size} 3D descriptors\n" if cl.option_present('v')

  makefile = File.open(File.join(tmpdir, 'Makefile'), "w")
  makefile << "all:"
  to_compute.each do |d|
    makefile << " std.#{d.name}"
  end
  makefile << "\n"

  if cl.option_present('v')
    verbose = '-v'
  else
    verbose = ""
  end

  # If standardisation not requested make a link from start.smi to the standardised file.
  if cl.option_present('nostd')
    system("cd #{tmpdir} && ln -s start.smi std.smi")
    #FileUtils.ln(File.join(tmpdir, 'start.smi'), File.join(tmpdir, 'std.smi'))
  else
    makefile << "std.smi:start.smi\n"
    makefile << "	fileconv.sh #{verbose} #{fileconv_options} -S std $<\n"
  end

  if ! d3.empty?
    makefile << "std.sdf: std.smi\n"
    makefile << "	rcorina.sh #{verbose} #{corina_options} $< > $@\n"
    d3.each do |d|
      makefile << "std.#{d.name}: std.sdf\n"
      makefile << "	#{d.programme} #{d.extra} $< > $@\n"
    end
  end

  d2.each do |d|
    makefile << "std.#{d.name}: std.smi\n"
    makefile << "	#{d.programme} #{d.extra} $< > $@\n"
  end

  makefile << "concat: all\n"
  makefile << "	concat_files.sh -I #{verbose}"
  to_compute.each do |d|
    makefile << " std.#{d.name}"
  end

  if cl.option_present('o')
    o = cl.value('o')
    if o.match(/^\//)
       output_file = o
    else
      output_file = File.join(Dir.pwd, o)
    end
    makefile << " > #{output_file}\n"
  end

  makefile.close
end

# Read descriptor classes from `fname` and add to `to_compute`.
def descriptors_from_file(fname, descriptors, to_compute)
  seen_here = Set.new()
  File.open(fname).each_line do |line|
    line.chomp.split.each do |d|
      d.gsub!(/_.*/, "")
      unless descriptors.key?(d)
        raise "Unrecognised descriptor class #{d}"
      end
      next if seen_here.include?(d)
      to_compute << descriptors[d]
      seen_here.add(d)
    end
  end
end

# Read descriptors from descriptor file `fname` and add to `to_compute`.
def descriptors_from_desctiptor_file(fname, descriptors, to_compute)
  seen_here = Set.new()
  File.open(fname).each_line do |line|
    line.chomp.split[1..-1].each do |d|
      d.gsub!(/_.*/, "")
      unless descriptors.key?(d)
        raise "Unrecognised descriptor class #{d}"
      end
      next if seen_here.include?(d)
      to_compute << descriptors[d]
      seen_here.add(d)
    end
    return
  end
end

def show_relative_speeds(descriptors)
  $stderr << "Seconds to compute 50k molecules. Lower numbers are better.\n";
  $stderr << "It takes corina 287 seconds to compute 3D structures\n"
  $stderr << "so all 3D descriptors should have 287 seconds added to their score.\n";
  descriptors.sort_by { |k, v| v.speed }.each do |k, v|
    $stderr << k << "\t";
    if v.threed
      $stderr << '*'
    else
      $stderr << ' '
    end
    $stderr << "\t" << v.speed << "\n"
  end
  $stderr << "Note that marvin ran 8 way parallel!\n"
  exit 0
end

def main
  # Corina takes 287 seconds to process 50k molecules
  descriptors = {}
  descriptors['abr'] = Descriptor.new('abraham.sh', false, 17)
  descriptors['ap'] = Descriptor.new('ap.sh', false, 10.6)
  descriptors['chgfp'] = Descriptor.new('chgfp.sh', false, 8.2)
  descriptors['cmi'] = Descriptor.new('cmi.sh', false, 24.4)
  descriptors['cnk'] = Descriptor.new('cnk.sh', false, 23.1)
  descriptors['comma'] = Descriptor.new('comma.sh', true, 31.1)
  descriptors['dbf'] = Descriptor.new('dbf.sh', true, 12.2)
  descriptors['estate'] = Descriptor.new('jwestate.sh', false, 4.4)
  descriptors['ghose'] = Descriptor.new('ghose_crippen.sh -a -a', false, 16.6)
  descriptors['ha'] = Descriptor.new('ha.sh', false, 14.8)
  descriptors['hb'] = Descriptor.new('hb.sh', false, 14.7)
  descriptors['hpo'] = Descriptor.new('hydrophobic_sections.sh', false, 4.0)
  descriptors['jurs'] = Descriptor.new('jurs.sh', true, 169.0)
  descriptors['jwc'] = Descriptor.new('jwcats.sh', false, 5.9)
  descriptors['jwdist'] = Descriptor.new('jwdist.sh', true, 9.7)
  descriptors['jwdp'] = Descriptor.new('jwdip.sh', true, 3.9)
  descriptors['jwmc'] = Descriptor.new('jwmolconn.sh', false, 62.1)
  descriptors['marvin'] = Descriptor.new('marvin.sh', false, 461)
  descriptors['medv'] = Descriptor.new('jwmedv.sh', false, 7.1)
  descriptors['mk'] = Descriptor.new('maccskeys.sh', false, 5.4)
  descriptors['morse'] = Descriptor.new('jwmorse.sh -f', true, 95.1)
  # descriptors['ms'] = Descriptor.new('unknown.sh', true)
  descriptors['pd'] = Descriptor.new('iwpathd.sh', false, 252)
  descriptors['sh'] = Descriptor.new('tshadow.sh', true, 11.9)
  descriptors['tt'] = Descriptor.new('topotorsion.sh -H 1600', false, 2.25)
  descriptors['w'] = Descriptor.new('w.sh', false, 15.2)

  # Set the name field of each descriptor
  descriptors.each do |k, v|
    v.name = k
  end

  # Build the options for the command line
  opts = ""
  descriptors.each do |key, v|
    opts << "-#{key}-#{key.upcase}=close"
  end
  opts <<"-o=s-c=s-all-v-nostd-keep-tmpdir=s-j=ipos-speed-dfile=sfile-dlist=sfile"

  # $stderr << opts << "\n"

  cl = IWCmdline.new(opts)

  if cl.unrecognised_options_encountered
    $stderr << "unrecognised_options_encountered\n";
    usage(descriptors)
  end

  if cl.option_present('speed')
    show_relative_speeds(descriptors)
  end

  verbose = cl.option_present('v')

  already_present = Set.new()

  to_compute = []

  cl.values('c').each do |c|
    c.split(',').each do |f|
      unless descriptors.key?(f)
        $stderr << "Unrecognised feature '#{f}'\n"
        return 1
      end
      raise "Duplicate #{f}" if already_present.include?(f)
      already_present.add(f)
      to_compute << descriptors[f]
    end
  end

  if cl.option_present('all')
    descriptors.each do |k, v|
      to_compute << v
    end
  end

  if cl.option_present('dlist')
    descriptors_from_file(cl.value('dfile'), descriptors, to_compute)
  end
  if cl.option_present('dfile')
    descriptors_from_desctiptor_file(cl.value('dfile'), descriptors, to_compute)
  end

  descriptors.each do |k, v|
    if cl.option_present(k)
      raise "Duplicate #{k}" if already_present.include?(k)
      already_present.add(k)
      to_compute << v
    elsif cl.option_present(k.upcase)
      # don't bother checking duplicates, too rare.
      v.extra = cl.value(k.upcase)
    end
  end

  if to_compute.empty?
    $stderr << "Nothing to compute\n"
    usage(descriptors)
  end

  $stderr << "Will compute #{to_compute.size} descriptor sets\n" if verbose

  if ARGV.empty?
    $stderr << "Insufficient arguments\n"
    usage(descriptors)
  end

  # If just one descriptor being computed, no need for a Makefile, just
  # a single command.
  if to_compute.size == 1
    return do_single_descriptor(to_compute[0], ARGV[0], cl)
  end

  # It would probably be doable to support multiple input files, but beware nostd...
  if ARGV.size > 1
    $stderr << "Only supports one input file\n"
    usage(descriptors)
  end

  if cl.option_present('tmpdir')
    tmpdir = cl.value('tmpdir')
  else
    tmpdir = "/tmp/make_descriptors_#{ENV['USER']}_#{Process.pid}"
  end

  if File.directory?(tmpdir)
    ;
  elsif FileUtils.mkdir_p(tmpdir)
    ;
  else
    raise "Cannot create #{tmpdir}"
  end

  FileUtils.cp(ARGV[0], File.join(tmpdir, 'start.smi'))

  create_makefile(to_compute, tmpdir, cl)

  j = 1
  if cl.option_present('j')
    j = cl.value('j')
  end

  if cl.option_present('o')
    silent = ""
  else
    silent = "--silent"
  end

  system("cd #{tmpdir} && make #{silent} -j #{j} concat")

  # If the user specified a tmpdir, do not delete it.
  exit 0 if cl.option_present('tmpdir')
  exit 0 if cl.option_present('keep')

  FileUtils.rm_rf(tmpdir)

  return 0
end

if RUBY_VERSION.match(/^2/)
  $stderr << "Ruby version too old #{RUBY_VERSION} needs at least version 3\n"
  exit 1
end

main
