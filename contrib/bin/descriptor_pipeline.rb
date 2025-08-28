#!/usr/bin/env ruby

# Computes molecular descriptors that have been enabled for pipelined processing

lillymol_home = ENV['LILLYMOL_HOME']

require_relative('lib/iwcmdline')

def usage
  $stderr << "Generates molecular descriptors in a pipelined manner\n"
  $stderr << "A subset of make_descriptors descriptor generators are supported - more may be added\n"
  $stderr << "While this can be used as an efficient descriptor generator, it was designed\n"
  $stderr << "to be used in conjunction with xgboost_model_evaluate which scores an existing\n"
  $stderr << "xgboost model via the C API, and is therefore fast.\n";
  $stderr << "Pipe the output of this tool to that and you can score molecules indefinitely\n"
  $stderr << "descriptor_pipeline.sh -w -abr file.smi | xgboost_model_evaluate ... -b 1000 -\n"
  $stderr << "Options are of two forms\n"
  $stderr << "  lowercase      takes no qualifiers\n";
  $stderr << "  uppercase      takes    qualifiers - open and close syntax\n";
  $stderr << "To compute default -w descriptors\n"
  $stderr << "   descriptor_pipeline -w file.smi\n"
  $stderr << "To pass some options to the underlying w.sh script\n"
  $stderr << "   descriptor_pipeline -W -O none -O alogp -W file.smi\n"
  $stderr << "The options within the opening and closing -W are passed to the underlying executable\n";
  $stderr << " -w -W ... -W            iwdescr interpretable descriptors\n"
  $stderr << " -sh -SH ... -SH         shadow descriptors (3D, calls rcorina)\n"
  $stderr << " -abr -ABR ... -ABR      abraham\n"
  $stderr << " -hpo -HPO ... -HPO      hydrophobic sections\n"
  $stderr << " -jwmc -JWMC ... -JWMC   molecular connectivity descriptors\n"
  $stderr << " -v               verbose output\n"

  exit(0)
end

def main
  cl = IWCmdline.new('-v-sh-SH=close-hpo-HPO=close-w-W=close-jwmc-JWMC=close-abr-ABR=close')
  if cl.unrecognised_options_encountered
    $stderr << "unrecognised_options_encountered\n"
    usage()
  end

  if ARGV.empty?
    $stderr << "Insufficient arguments\n"
    usage
  end

  # We need to know which descriptor set is the last. This depends
  # on the order in which we examine the options
  if cl.option_present('jwmc') || cl.option_present('JWMC')
    last = 'jwmc'
  elsif cl.option_present('w') || cl.option_present('W')
    last = 'w'
  elsif cl.option_present('hpo') || cl.option_present('HPO')
    last = 'hpo'
  elsif cl.option_present('abr') || cl.option_present('ABR')
    last = 'abr'
  elsif cl.option_present('sh') || cl.option_present('SH')
    last = 'sh'
  else
    $stderr << "NO descriptors requested\n";
    return 1
  end

  gen_3d = cl.option_present('sh') || cl.option_present('SH')

  verbose = cl.option_present('v')

  is_first = true

  if gen_3d
    cmd = "rcorina.sh -u #{ARGV[0]}"
    is_first = false
  else
    cmd = ""
  end

  prefix='new_'

  # $stderr << "Last is #{last}\n"

  if cl.option_present('sh') || cl.option_present('SH')
    cmd << '|' unless cmd.empty?
    cmd << "#{prefix}tshadow.sh -i smi"
    cmd << ' ' << cl.value('SH') if cl.option_present('SH')
    cmd << ' -B wpipe' unless last == 'sh'
    if is_first
      cmd << ' ' << ARGV[0]
      is_first = false
    else
      cmd << ' -'
    end
  end

  if cl.option_present('abr') || cl.option_present('ABR')
    cmd << '|' unless cmd.empty?
    cmd << "#{prefix}abraham.sh -Y flush"
    cmd << ' ' << cl.value('ABR') if cl.option_present('ABR')
    cmd << ' -Y wpipe' unless last == 'abr'
    if is_first
      cmd << ' ' << ARGV[0]
      is_first = false
    else
      cmd << ' -Y rpipe -'
    end
  end

  if cl.option_present('hpo') || cl.option_present('HPO')
    cmd << '|' unless cmd.empty?
    cmd << "#{prefix}hydrophobic_sections.sh -Y flush"
    cmd << ' ' << cl.value('HPO') if cl.option_present('HPO')
    cmd << ' -Y wpipe' unless last == 'hpo'
    if is_first
      cmd << ' ' << ARGV[0]
      is_first = false
    else
      cmd << ' -Y rpipe -'
    end
  end

  if cl.option_present('w') || cl.option_present('W')
    cmd << '|' unless cmd.empty?
    cmd << "#{prefix}iwdescr.sh -B flush"
    cmd << ' ' << cl.value('W') if cl.option_present('W')
    cmd << ' -B wpipe' unless last == 'w'
    if is_first
      cmd << ' ' << ARGV[0]
      is_first = false
    else
      cmd << ' -B rpipe -'
    end
  end

  if cl.option_present('jwmc')
    cmd << '|' unless cmd.empty?
    cmd << "#{prefix}jwmolconn.sh"
    cmd << ' -Y wpipe' unless last == 'jwmc'
    if is_first
      cmd << ' ' << ARGV[0]
      is_first = false
    else
      cmd << ' -Y rpipe -'
    end
  end

  $stderr << "Executing\n#{cmd}\n" if verbose
  $stderr << "#{cmd} failed" unless system(cmd)
end


main
