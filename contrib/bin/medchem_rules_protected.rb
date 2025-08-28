#!/usr/bin/env ruby

# A common task is to run the medchem rules on a set of molecules that
# come with an undesirable functional group. The typical trick is to
# add some benign group to that, often a phenyl, and run the rules on
# those 'protected' molecules

require_relative('lib/iwcmdline_v2')
require_relative('reaction_pb')

def usage
  exit(0)
end

def main
  cl = IWCmdlineV2.new('-v-s=s-rm=int-j=ipos-protect=sfile-z=s-trxn=close-tp1_pipe=close')

  verbose = cl.option_present('v')

# if cl.unrecognised_options_encountered
#   $stderr << "unrecognised_options_encountered\n"
#   usage
# end

  if ARGV.empty?
    $stderr << "Must specify input file\n"
    usage
  end

  unless cl.option_present('s')
    $stderr << "Must specify smarts of the group to be protected/removed via the -s option\n"
    usage
  end

  if ARGV.empty?
    $stderr << "Must specify input file as argument\n"
    usage
  end
  input = ARGV[0]

  smarts = cl.value('s')

  if cl.option_present('protect')
    fname = cl.value('protect')
    protect = File.read(fname)
    protect
  else
    protect = '[1cH]1ccccc1'
  end

  atoms_to_remove = if cl.option_present('rm')
                      cl.values('rm')
                    else
                      []
                    end

  join = if cl.option_present('j')
              cl.value('j')
            else
              0
            end

  if atoms_to_remove.include?(join)
    $stderr << "Atoms to remove #{atoms_to_remove} includes join #{join}, impossible\n"
    exit(1)
  end

  files_to_remove = []

  tmprxn = "/tmp/mcrp#{Process.uid}.#{Process.pid}.dat"  # must end in .dat so trxn knows it is binary
  files_to_remove << tmprxn
  
  rxn = ReactionProto::Reaction.new
  rxn.scaffold = ReactionProto::ScaffoldReactionSite.new
  rxn.scaffold.id = 0
  rxn.scaffold.smarts << smarts
  atoms_to_remove.each do |a|
    rxn.scaffold.remove_atom << a
  end

  rxn.sidechain << ReactionProto::SidechainReactionSite.new
  rxn.sidechain[0].id = 1
  rxn.sidechain[0].reagent << protect
  rxn.sidechain[0].smarts << "[1]"

  rxn.sidechain[0].join << ReactionProto::InterParticleBond.new()
  rxn.sidechain[0].join[0].a1 = join
  rxn.sidechain[0].join[0].a2 = 0

# $stdout << rxn << "\n"

  File.open(tmprxn, 'wb') do |file|
    file.write(ReactionProto::Reaction.encode(rxn))
  end

  okmedchem = "/tmp/#{Process.uid}#{Process.pid}.smi"
  files_to_remove << okmedchem

  cmd = "trxn -W NONE -S - -P #{tmprxn}"
  cl.values('z').each do |z|
    cmd << " -z #{z}"
  end
  cmd << ' ' << cl.value('trxn') if cl.option_present('trxn')
  cmd << " #{input}"
  cmd << "| tp1_pipe.sh -okiso"
  cmd << ' ' << cl.value('tp1_pipe') if cl.option_present('tp1_pipe')
  cmd << " - > #{okmedchem}"
  $stderr << "Executing #{cmd}\n" if verbose
  unless system(cmd)
    $stderr << "Warning, non zero rc from #{cmd}\n"
  end

  unless File.size?(okmedchem)
    $stderr << "#{cmd} failed #{okmedchem} missing or empty\n";
    return 1
  end

  cmd = "fetch_smiles_quick.sh -c 2 #{okmedchem} #{input}"
  system(cmd)

  files_to_remove.each do |fname|
    File.unlink(fname)
  end
end

main
