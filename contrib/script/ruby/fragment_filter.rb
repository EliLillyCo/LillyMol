#! /usr/bin/env ruby

ianhome = ENV['LILLYMOL_HOME']

require_relative "lib/iwcmdline.rb"

keep_tmp_files = false

lower_atom_count_cutoff = 7

pass_dash_f_option_to_fragment_filter = false

$expert = false

def usage (rc)
  $stderr.print "Implementation of fragment filtering\n"
  $stderr.print " -tp1 ... -tp1  tp1_pipe.sh options (default '#{$tp1_pipe_options}')\n"
  $stderr.print " -c <natoms>    lower atom count cutoff\n" if ($expert)
  $stderr.print " -C <natoms>    upper atom count cutoff\n" if ($expert)
  $stderr.print " -dmrt <fname>  run tp1_summarise and write demerit data to <fname>\n"
  $stderr.print " -ff .... -ff   replace fragment_filter options\n" if ($expert)
  $stderr.print " -ffa ... -ffa  append to fragment_filter options\n" if ($expert)
  $stderr.print " -h <frac>      min heteroatom fraction\n"
  $stderr.print " -H <frac>      max heteroatom fraction\n"
  $stderr.print " -X <frac/num>  max fraction or number of any particular heteroatom\n"
  $stderr.print " -cn ... -cn    common_names options (default '#{$common_names_options}')\n" if ($expert)
  $stderr.print " -ZOF <script>  ZOF script\n" if ($expert)
  $stderr.print " -XC <fname>    smiles file of excluded counterions\n";
  $stderr.print " -relax <n>     allow <n> extra atoms to all ring sizes\n"
  $stderr.print " -f             pass the -f option to fragment_filter and perform subsequent processing\n" if ($expert)
  $stderr.print " -keep          keep temporary files\n" if ($expert)
  $stderr.print " -expert        more options\n" unless ($expert)
  $stderr.print " -v             verbose output\n"
  exit(rc)
end

cl = IWCmdline.new("-v-tmpdir=s-tp1=close-ff=close-ffa=close-keep-cn=close-dmrt=s-expert-ZOF=xfile-f-c=ipos-C=ipos-XC=sfile-relax=ipos-h=fraction-H=fraction-X=f")

$expert = cl.option_present('expert')

if cl.unrecognised_options_encountered()
  $stderr.print "Unrecognised options encountered\n"
  usage(1)
end

$verbose = cl.option_count('v')

if (0 == ARGV.size)
  $stderr.print "Insufficient arguments\n"
  usage(2)
end

tmpdir = "."

if (cl.option_present('tmpdir'))
  tmpdir = cl.value('tmpdir')

  if (! FileTest.directory?(tmpdir))
    Dir.mkdir(tmpdir)
  end
end

if (cl.option_present('c'))
  lower_atom_count_cutoff = cl.value('c')
end

atom_counts = "-M 0=14 -M 1=12 -M 2=16 -M 3=20"

if (cl.option_present('relax'))
  a = cl.value('relax').to_i

  atom_counts = "-M 0=#{14+a} -M 1=#{12+a} -M 2=#{16+a} -M 3=#{20+a}"
end

$tp1_options = "-c #{lower_atom_count_cutoff}"

min_heteroatom_fraction = 0.15
max_heteroatom_fraction = 0.70

min_heteroatom_fraction = cl.value('h') if cl.option_present('h')
max_heteroatom_fraction = cl.value('H') if cl.option_present('H')

if (min_heteroatom_fraction >= max_heteroatom_fraction)
  $stderr << "Inconsistent min #(min_heteroatom_fraction) and max #(max_heteroatom_fraction) heteroatom fractions\n";
  usage(2)
end

fragment_filter_options = "-c #{lower_atom_count_cutoff} -G 2 -E autocreate -b 3 -A D -A I -g all -l -g ltltr -h #{min_heteroatom_fraction} -H #{max_heteroatom_fraction} #{atom_counts} -Y 2 -r 1 -R 3 -D 3 -a -s 4 -p 1 -B failfrag"

if (cl.option_present('X'))
  fragment_filter_options << ' -X ' << cl.value('X').to_s
end

$common_names_options = "-T I=Cl -T Br=Cl -v -c -l"

if (cl.option_present('fc'))
  $fileconv_options = cl.value('fc')
end

if (cl.option_present('tp1'))
  $tp1_pipe_options = cl.value('tp1')
end

if (cl.option_present('ff'))
  fragment_filter_options = cl.value('ff')
end

if (cl.option_present('ffa'))
  f = cl.value('ffa')
  fragment_filter_options = "#{fragment_filter_options} #{f}"
end

if (cl.option_present('XC'))
  x = cl.value('XC')

  fragment_filter_options << " -J #{x}"
end

if (cl.option_present('cn'))
  $common_names_options = cl.value('cn')
end

if (cl.option_present('f'))
  pass_dash_f_option_to_fragment_filter = false
  fragment_filter_options << ' -f'
end

keep_tmp_files = cl.option_present('keep')

files_to_be_removed = Array.new

input_files = ""

ARGV.each do |fname|
  raise "MIssing or empty input file '#{fname}'" unless (FileTest.size?(fname))
  input_files << " #{fname}"
end

fragment_filter = 'fragment_filter.sh'
common_names = 'common_names.sh'
tp1_summarise = 'tp1_summarise.sh'
fileconv = 'fileconv.sh'

good_stem = 'okf'

if (cl.option_present('g'))
  good_stem = cl.value('g');
end

bad_stem = 'badf'

if (cl.option_present('b'))
  bad_stem = cl.value('b');
end

def execute_system(cmd, must_be_present = false)
  if ($verbose)
    $stderr.print "Executing '#{cmd}'\n"
  end

  system(cmd)

  return true unless must_be_present

  return true if (File.size?(must_be_present))

  $stderr << "'#{cmd}' failed or produced no output, cannot continue\n"
  exit 0
end

# Kind of tricky to get the order of things right. First run fileconv
# to get rid of some undesirable fragment things

ok1 = "#{tmpdir}/#{good_stem}_S1"
bad1 = "#{tmpdir}/#{bad_stem}_B1"

upper_atom_count_cutoff = 18

upper_atom_count_cutoff = cl.value('C') if (cl.option_present('C'))

cmd = "#{fileconv} -c #{lower_atom_count_cutoff} -Y aclf -g all -C #{upper_atom_count_cutoff} -F 2 -L #{bad1} -f rmxt=7 -O def -V -E autocreate -I 0 -S #{ok1} #{input_files}"

ok1 << ".smi"

execute_system(cmd, ok1)

files_to_be_removed.push(ok1)

# In order to get bad counterions recognised, we need to run the medchem
# rules before anything else

okmedchem = "#{tmpdir}/#{good_stem}_1.smi"

files_to_be_removed.push(okmedchem)

cmd = "#{ianhome}/bin/tp1_pipe.sh "

if (! cl.option_present('dmrt'))
  cmd << "-noapdm "
end

cmd << "#{$tp1_options} -B #{bad_stem} #{ok1} > #{okmedchem}"

execute_system(cmd, okmedchem)

if (cl.option_present('dmrt'))
  dmrt_fname = cl.value('dmrt')

  cmd = "#{tp1_summarise} -B #{bad_stem} -D #{okmedchem} > #{dmrt_fname}"

  system(cmd)

  unless (FileTest.size?(dmrt_fname))
    $stderr.print "Warning, '#{cmd}' failed to produce output\n"
  end 

  tmpx = "#{good_stem}_nd.smi"

  cmd = "sed 's/: D([0-9][0-9]*).*//' #{okmedchem} > #{tmpx}"
  system(cmd)

  if (! FileTest.size?(tmpx))
    $stderr.print "'#{cmd}' failed\n" 
    exit 1
  end

  File.rename(tmpx, okmedchem)
end

if ($verbose)
  system("wc #{okmedchem} >&2")
end

# Now run the fragment filter

tmp2 = "#{tmpdir}/#{good_stem}_2"
bad2 = "#{tmpdir}/#{bad_stem}_f2"

files_to_be_removed.push(tmp2)

cmd = "#{fragment_filter} #{fragment_filter_options} -B #{bad2} -S #{tmp2} #{okmedchem}"

tmp2 << ".smi"

execute_system(cmd, tmp2)

system("wc #{tmp2} >&2") if ($verbose)

# Now substructure filters

cmd = "#{ianhome}/bin/ZOF.v2.sh"
if (cl.option_present('ZOF'))
  tmp = cl.value('ZOF')
  cmd = "#{tmp}"
end

cmd << " #{tmp2}"

execute_system(cmd)

if (! FileTest.size?('okZOF.smi'))
  $stderr.print "'#{cmd} produced no output\n" 
  exit 1
end

if (pass_dash_f_option_to_fragment_filter)
  tmp3 = "#{tmpdir}/#{good_stem}_3"

  cmd = "#{common_names} #{$common_names_options} -S #{tmp3} okZOF.smi"

  execute_system(cmd)

  tmp3 << ".smi"

  files_to_be_removed.push(tmp3)

  if (! FileTest.size?(tmp3))
    $stderr.print "'#{cmd} produced no output\n"
    exit 2
  end

  cmd = "#{ianhome}/bin/large_counterion_first.rb #{tmp3}"

  system(cmd)
else
  system("/bin/cat okZOF.smi")
end

unless (keep_tmp_files)
  files_to_be_removed.each do |f|
    File.unlink(f) if (FileTest.exists?(f))
  end
end
