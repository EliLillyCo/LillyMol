#!/usr/bin/env ruby

# Replace a sidechain with sidechains extracted from Chembl

require_relative 'lib/iwcmdline'

def usage
  help = <<-HELP
Adds Chembl aromatic substituents to molecules.
LillyMol includes a file of aromatic substituents extracted from Chembl 35

data/chembl_sidechains.smi

That file includes isotopically labelled sidechains as well as the number of times the
sidechain is found in Chembl.

If you have specific requirements on what kinds of substituents to use, use
molecule_filter on that file to impose restrictions, and then use the -substituents option to this script.

This tool can either replace an existing substituent, or add substituents to a site that currently has a Hydrogen substituent.
Use -smarts1 'smarts'   to specify the smarts of an atom to which substituents are added. '[$([aD2H](:[aD2]):[aD2])]' for example.
Use -smarts2 'smarts'   to specify two atoms. The bond between them is broken, the second atom removed,
                        and the sidechain attached to the first atom. '[$([aD3](:[aD2]):[aD2])]-!@{a<10}[D2]' for example.
In either case '-smarts1 def' or '-smarts2 def' will apply the default smarts listed above.
Options
 -sidechains <fname>    file of sidechains, by default LILLYMOL_HOME/data/chembl_sidechains.textproto.
 -smarts1 <smarts>      smarts for an unsubstitued atom to which the Chembl sidechains are attached.
 -smarts2 <smarts>      smarts for a bond that will be broken. The new sidechain will be attached to the first atom
                        and the fragment defined by the second atom will be removed.
 -smarts{1,2} def       use default smarts listed above.
 -support <n>           only use sidechains that occur <n> or more times in the sidechains file.
 -S <fname>             output file to produce, stdout by default.
 -trxn ... -trxn        options passed directly to trxn.
 -v                     verbose output
HELP
  $stderr << help << "\n"

  exit(0)
end

# If the -sidechains option is specified, return that.
# Otherwise look in LILLYMOL_HOME for chembl_sidechains.textproto
def get_sidechains(cl)
  return cl.value('sidechains') if cl.option_present('sidechains')

  lillymol_home = if ENV['LILLYMOL_HOME']
                    "#{ENV['LILLYMOL_HOME']}"
                  else
                    File.dirname(File.dirname(File.dirname(__FILE__)))
                  end

  sidechains = "#{lillymol_home}/data/chembl_sidechains.textproto"

  return sidechains if File.size?(sidechains)

  unless ENV['LILLYMOL_HOME']
    $stderr << "Must specify a value for LILLYMOL_HOME in order to locate Chembl sidechains file\n"
    $stderr << "export LILLYMOL_HOME=/path/to/LillyMol\n"
    exit 1
  end

  $stderr << "Did not find sidechains file #{sidechains}, is LILLYMOL_HOME set properly?\n"
  exit 1
end

def main
  cl = IWCmdline.new("-v-sidechains=sfile-support=ipos-smarts1=s-smarts2=s-trxn=close-S=s")
  if cl.unrecognised_options_encountered
    $stderr << "unrecognised_options_encountered\n"
    usage
  end

  verbose = cl.option_present('v')

  if ARGV.empty?
    $stderr << "Must specify starting smiles file(s) as argument(s)\n"
    usage
  end

  if cl.option_present('smarts1') && cl.option_present('smarts2')
    $stderr << "Sorry, only one of -smarts1 or -smarts2 can be specified.\n"
    usage
  end

  if cl.option_present('smarts1')
    s = cl.value('smarts1')
    if s == 'def'
      smarts1 = '[$([cD2H](:[aD2]):[aD2])]'
    else
      smarts1 = cl.value('smarts1')
    end
  end

  if cl.option_present('smarts2')
    s = cl.value('smarts2')
    if s == 'def'
      smarts2 = '[$([cD3](:[aD2]):[aD2])]-!@{a<10}[D2]'
    else
      smarts2 = cl.value('smarts2')
    end
  end

  if ! cl.option_present('smarts1') && ! cl.option_present('smarts2')
    smarts1 = '[$([cD2H](:[aD2]):[aD2])]'
  end

  reaction1 = <<-END_REACTION1
name: "add_sidechain"
scaffold {
  id: 0
  smarts: "#{smarts1}"
}
sidechain {
  id: 1
  smarts: "[1]"
  join {
    a1: 0
    a2: 0
    btype: SS_SINGLE_BOND
  }
}
END_REACTION1

  reaction2 = <<-END_REACTION2
name: "replace_sidechain"
scaffold {
  id: 0
  smarts: "#{smarts2}"
  break_bond {
    a1: 0
    a2: 1
  }
  remove_fragment: 1
}
sidechain {
  id: 1
  smarts: "[1]"
  join {
    a1: 0
    a2: 0
    btype: SS_SINGLE_BOND
  }
}
END_REACTION2

  rxnfile = File.join('/tmp', "replace_sidechain#{Process.pid}.textproto")

  files_to_delete = []

# $stderr << "Smarts #{smarts1} #{smarts2}\n"

  if smarts1
    if File.write(rxnfile, reaction1) != reaction1.length
      $stderr << "Could not write reaction #{rxnfile}\n"
      return 1
    end
  else
    if File.write(rxnfile, reaction2) != reaction2.length
      $stderr << "Could not write reaction #{rxnfile}\n"
      return 1
    end
  end

  files_to_delete << rxnfile

  sidechains = get_sidechains(cl)

  if cl.option_present('support')
    support = cl.value('support')
    tmpsidechains = File.join('/tmp', "sidechains#{Process.pid}.textproto")
    cmd = "gawk '$11 >= #{support}' #{sidechains} > #{tmpsidechains}"
    $stderr << "Executing #{cmd}\n" if verbose
    system(cmd)
    sidechains = "#{tmpsidechains}"
    files_to_delete << sidechains
    if verbose
      lines = `wc -l #{sidechains}`
      $stderr << "With support #{support} have #{lines} substituents\n"
    end
  end

  # cmd = "trxn.sh -z i -m each -P #{rxnfile} -J rpt=100000 -J nomshmsg"
  cmd = "/lrlhps/users/rx87690/LillyMolPrivate/src/bazel-bin/Molecule_Tools/trxn -z i -m each -P #{rxnfile} -J rpt=100000 -J nomshmsg"

  cmd << ' ' << cl.value('trxn') if cl.option_present('trxn')

  if cl.option_present('S')
    s = cl.value('S')
    cmd << " -S #{s}"
  else
    cmd << ' -S -'
  end
  cmd << ' -v' if verbose
  cmd << " #{ARGV[0]} #{sidechains}"
  $stderr << "Executing #{cmd}\n" if verbose
  system(cmd)

  files_to_delete.each do |fname|
#   File.unlink(fname)
  end

  return 0
end

main
