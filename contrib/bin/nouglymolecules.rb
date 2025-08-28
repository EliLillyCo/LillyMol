#!/usr/bin/env ruby

lillymol_home = ENV['LILLYMOL_HOME']
datadir = "#{lillymol_home}/data/nouglymolecules"

require_relative('lib/iwcmdline')

$expert = false

def usage(rc)
  $stderr << "Gets rid of \"ugly\" molecules. These are defined as:\n";
  $stderr << "  molecules having too many rotatable bonds\n";
  $stderr << "  molecules having spiro fusions\n";
  $stderr << "  molecules with too many rings\n";
  $stderr << "  molecules with ring systems containing 4 or more fused rings\n";
  $stderr << "  molecules with complex fused ring structures\n";
  $stderr << "\n";

  $stderr << "The following options are recognised\n";
  $stderr << " -min_rotb <int>      minimum number of rotatable bonds needed\n";
  $stderr << " -max_rotb <int>      maximum number of rotatable bonds needed\n";
  $stderr << " -max_flxc <int>      maximum number of rotatable bonds along a flexible chain\n";
  $stderr << " -max_atom <int>      maximum number of atoms\n";
  $stderr << " -breaku              flexible chains don't pass through unsaturated atoms\n";
  $stderr << " -maxc <nn>           maximum connectivity along a flexible chain\n";
  $stderr << "                      use '-maxc 2' to to follow flexible chains following branches\n";
  $stderr << " -min_rings <int>     minimum number of rings\n";
  $stderr << " -max_rings <int>     maximum number of rings\n";
  $stderr << " -min_ring_size <int> minimum allowable size of a ring\n";
  $stderr << " -max_ring_size <int> maximum allowable size of a ring\n";
  $stderr << " -max_sys_size <int>  maximum number of rings in a ring system (default $max_sys_size)\n";
  $stderr << " -okspiro             allow spiro fusions\n";
  $stderr << " -nospirosys          suppress spiro fusions if they are part of a ring system\n";
  $stderr << " -symm <n>:<f>        discard if a molecule has at least <n> atoms and <f> of atoms are symmetric\n"
  $stderr << "                        e.g. -symm 20:0.95\n"
  $stderr << " -smarts <smt>        extra smarts to discard\n"
  $stderr << " -long <n>:<d>        discard if a molecule has at least <n> atoms and is a long chain\n";
  $stderr << " -maxchiral <n>       discard if more than <n> chiral centres\n"
  $stderr << " -lib <dir>           directory with queries\n" if $expert
  $stderr << " -expert              more options\n"
  $stderr << " -v                   verbose output\n";

  exit (rc);
end

# IWCmdlineV2  if using v2

cl = IWCmdline.new("-v-expert-keep-min_rotb=i-max_rotb=ipos-max_flxc=i-max_atom=ipos-min_atom=ipos-breaku-maxc=ipos-min_rings=ipos-max_rings=ipos-max_ring_size=ipos-max_sys_size=ipos-okspiro-nospirosys-lib=dir-symm=s-long=s-smarts=s-maxchiral=i")

$expert = cl.option_present('expert')

if cl.unrecognised_options_encountered()
  $stderr.print "Unrecognised options encountered\n"
  usage(1)
end

verbose = cl.option_present('v')

# set this to '.sh' if using shell wrappers
sh = ""
tsubstructure = "tsubstructure#{sh}"
rotatable_bonds = "rotatable_bonds#{sh}"
fileconv = "fileconv#{sh}"
tsymmetry = "tsymmetry#{sh}"
long_molecules = "long_molecules#{sh}"

max_sys_size = 5

if cl.option_present('max_sys_size')
  max_sys_size = cl.value('max_sys_size')
end

if cl.option_present('lib')
  datadir = cl.value('lib')
end

if 0 == ARGV.size
  $stderr.print "Insufficient arguments\n"
  usage(2)
end


spiro = "#{datadir}/spiro.qry";
raise "Where is '$spiro'" unless File.size?(spiro);

spiro_sys = "#{datadir}/spiro_sys.qry";
raise "Where is '$spiro_sys'" unless File.size?(spiro_sys);

large_ring_system = "#{datadir}/large_ring_system_#{max_sys_size}.qry";
raise "Where is '#{large_ring_system}'" unless File.size?(large_ring_system);

strongly_fused = "#{datadir}/strongly_fused.qry";
raise "Where is '#{strongly_fused}'" unless File.size?(strongly_fused);

polyfluorinated = "#{datadir}/polyfluorinated.qry";
raise "Where is '#{polyfluorinated}'" unless File.size?(polyfluorinated);

def append_option_value_if_present(cl, flag1, flag2, cmd)
  if ! cl.option_present(flag1)
    return cmd
  end

  t = cl.value(flag1)

  return cmd << ' ' << flag2 << ' ' << t.to_s
end

max_rings = 7
if cl.option_present('max_rings')
  max_rings = cl.value('max_rings')
end

cmd = "#{fileconv} -V -L ugly_ring_count -i smi -S - ";
cmd << " -R #{max_rings}"


cmd = append_option_value_if_present(cl, 'max_atom', '-C', cmd)
cmd = append_option_value_if_present(cl, 'min_atom', '-c', cmd)

cmd = append_option_value_if_present(cl, 'min_rings', '-r', cmd)

if cl.option_present('maxchiral')
  m = cl.value('maxchiral')
  cmd << " -s maxc=#{m}"
end

ARGV.each do |fname|
  if ! File.size?(fname)
    $stderr << "Missing or empty file '#{fname}' ignored\n"
    next
  end
  cmd << ' ' << fname
end

min_rot_bonds = cl.value('min_rotb')
max_rot_bonds = cl.value('max_rotb')
max_flxc = cl.value('max_flxc')
maxc = cl.value('maxc')
breaku = cl.option_present('breaku')

#$stderr << "Values #{min_rot_bonds} #{max_rot_bonds} #{breaku} #{max_flxc}\n"
if min_rot_bonds || max_rot_bonds || max_flxc || maxc
  cmd << "|#{rotatable_bonds} -e";
  cmd << " -m #{min_rot_bonds-1}" if min_rot_bonds && min_rot_bonds > 0
  cmd << " -M #{max_rot_bonds+1}" if max_rot_bonds && max_rot_bonds > 0
  cmd << " -f #{max_flxc}" if max_flxc
  cmd << " -f maxc=#{maxc}" if maxc
  cmd << " -f breaku" if breaku
  cmd << " -S - -B ugly_rotbond -i smi -";
end

min_ring_size = cl.value('min_ring_size')
max_ring_size = cl.value('max_ring_size')

okspiro = cl.option_present('okspiro')

cmd  << "| #{tsubstructure} -m QDT ";
cmd  << "-q #{spiro} "  unless (okspiro);
cmd  << "-q #{spiro_sys} " if cl.option_present('nospirosys')
cmd  << "-q #{large_ring_system} -q #{strongly_fused} -m ugly_ring ";
cmd  << "-q #{polyfluorinated} " unless cl.option_present("okfluorine")

cmd  << "-s '[r3]' " if (3 == min_ring_size);
cmd  << "-s '[r<5]' " if (4 == min_ring_size);

cmd  << "-s '[r>#{max_ring_size}]' " if max_ring_size

if cl.option_present('smarts')
  cl.values('smarts').each do |s|
    cmd << "-s '#{s}' "
  end
end

cmd  << "-i smi -n - -";

if cl.option_present('symm')
  s = cl.value('symm')
  m = /^(\d+):(\S+)/.match(s)
  if m
    natoms = m[1]
    symm = m[2]
  else
    natoms = '20'
    symm = s
  end

  cmd << "| #{tsymmetry} -W 'minatm=#{natoms}' -W 'fsymmat<=#{symm}' -n ugly_symmetry -m - -i smi -"
end

if cl.option_present('long')
  s = cl.value('long')
  m = /^(\d+):(\S+)/.match(s)
  if m
    natoms = m[1]
    l = m[2]
  else
    matoms = 20
    l = s
  end
  cmd << "| #{long_molecules} -i smi -m #{natoms} -D #{l} -F ugly_long -P - -"
end

$stderr << "Executing '#{cmd}'\n" if (verbose);

system (cmd);
