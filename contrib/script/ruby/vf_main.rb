# view chemical structures, update to vf

require 'set'

require_relative 'lib/iwcmdline'

def usage(expert)
  $stderr << "Views smiles in a file or on the command line\n"
  $stderr << " -cactvs         use cactvs csps to display structures\n"
  $stderr << " -mview          use Chemaxon mview to display structures\n"
  $stderr << " -rdkit          use RDKt to display structures\n"
  $stderr << " -l <number>     view Lilly number <number>\n"
  $stderr << " -L <file>       view file of Lilly numbers\n"
  $stderr << " -S <smiles>     specify smiles on command line\n"
  $stderr << " -c <cols>       number of columns of structures\n"
  $stderr << " -r <rows>       number of rows    of structures\n"
  $stderr << " -n <number>     number of structures to show (default 30)\n"
  $stderr << " -a              all smiles in the file, rather than first 30\n"
  $stderr << " -p <prob>       sample records with probability <prob> (0,1)\n"
  $stderr << " -f <number>     skip first <number> smiles in file\n"
  $stderr << " -Z              fast random selection throughout file\n"
  $stderr << " -e <number>     process every <number> smiles (<number> = 2 means every 2nd smiles)\n"
  $stderr << " -j <number>     for each record selected, also show next <number> smiles\n" if expert
  $stderr << " -N <...>        pass to numbered_smiles with -n <...> (-N 1 for atom numbering)\n" if expert
  $stderr << " -grep <pat>     only select records that regular expression        match <pat>\n"
  $stderr << " -grepv <pat>    only select records that regular expression DO NOT match <pat>\n"
  $stderr << " -title <...>    title for each page (CACTVS)\n"
  $stderr << " -P <fname>      file name for the PostScript file (CACTVS)\n" if expert
  $stderr << " -pdf            convert the CACTVS PostScript file to pdf\n" if expert
  $stderr << " -X <pixels, -Y <pixels>   geometry for RDKit and Marvin\n"
  $stderr << " -WIDE           landscape mode (CACTVS)\n"
  $stderr << " -BOW -WOB -COW -COB colours. WOB = White on Black, COB = Colour on Black (CACTVS)\n"
  $stderr << " -CACTVS ... -CACTVS options passed directly to csps (CACTVS)\n" if expert
  $stderr << " -N <number>     call numbered_smiles to label the atoms\n" if expert
  $stderr << " -tdt            input is a TDT file (.gfp for example)\n"
  $stderr << " -tmp <dir>      use <dir> for temporary files\n" if expert
  $stderr << " -expert         expert mode, more options\n"
  $stderr << " -v              verbose output\n"
  exit 1
end

class FileDeleter
  def initialize
    @to_delete = []
  end
  def delete(fname)
    @to_delete << fname
  end
end

def files_that_exist(files)
  result = []
  files.each do |fname|
    if fname == '-'
      result << fname
      next
    end
    if File.size?(fname)
      result << fname
      next
    end
    with_ext = "#{fname}.smi"
    if File.size?(with_ext)
      result << with_ext
      next
    end
    raise "Missing or empty input file #{fname}"
  end
  result
end

# Execute system(cmd), maybe logging and returning success/fail
def execute_command(cl, cmd)
  $stderr << "Executing #{cmd}\n" if cl.option_present('v')
  return true if system(cmd)
  $stderr << "#{cmd} failed\n"
  false
end

def vf_rdkit(cl, smiles_file)
  cmd = []
  cmd << 'python'
  cmd << File.join(File.dirname(__FILE__), '../py/smiles2png.py')
  cmd << '-align' << cl.value('s') if cl.option_present('s') 
  cmd << '-x' << cl.value('X') if cl.option_present('X') 
  cmd << '-y' << cl.value('Y') if cl.option_present('Y') 
  cmd << '-mols_per_row' << cl.value('c') if cl.option_present('c')
  cmd << '-keep' if cl.option_present('k')
  cmd << '-verbose' if cl.option_present('v')
  cmd << '--input' << smiles_file
  execute_command(cl, cmd.join(' '))
end

def vf_mview(cl, smiles_file)
  mview = 'mview'
  cmd = []
  cmd << 'mview'
  if cl.option_present('r') && cl.option_present('c')
    cmd << '-r' << cl.value('r')
    cmd << '-c' << cl.value('c')
  elsif cl.option_present('spreadsheet')
    cmd << '--spreadsheet'
  end
  if cl.option_present('X') && cl.option_present('Y')
    cmd << '--geometry'
    cmd << "#{cl.value('X')}x#{cl.value('Y')}"
  end
  cmd << smiles_file
  execute_command(cl, cmd.join(' '))
end

def lines_in_file(fname)
  File.foreach(fname).count
end

# Return a hash containing :nrows and :ncols given that there are
# nlines of smiles to be shown
def rows_and_cols(cl, nlines)
  if cl.option_present('r') && cl.option_present('c')
    return {:nrows => cl.value('r'), :ncols => cl.value('c')}
  end

  return {:nrows => 1, :ncols => 1} if nlines == 1

  if nlines <= 3
    if cl.option_present('r')
      nrows = cl.value('r')
      nrows = nlines if nrows > nlines
      ncols = nlines / nrows
      ncols = 1 if ncols == 0
    elif cl.option_present('c')
      ncols = cl.value('c')
      ncols = nlines if ncols > nlines
      nrows = nlines / ncols
      nrows = 1 if nrows == 0
    else
      nrows = nlines
      ncols = 1
    end
    return {:nrows => nrows, :ncols => 1} if cl.option_present('r')
  end

  if cl.option_present('r')
    nrows = cl.value('r')
    ncols = nlines / nrows
    if ncols == 0
      ncols = 1
    elsif ncols > 3
      ncols = 3
    end
    return {:nrows => nrows, :ncols => ncols}
  end

  if cl.option_present('c')
    ncols = cl.value('c')
    $stderr << "Read #{nlines} lines, ncols #{ncols}\n"
    nrows = nlines / ncols
    if nrows == 0
      nrows = 1
    elsif nrows > 3
      nrows = 5
    end
    return {:nrows => nrows, :ncols => ncols}
  end

  # neither specified
  return {:nrows => 5, :ncols => 3} if nlines >= 15
  return {:nrows => 3, :ncols => 3} if nlines >= 9
  return {:nrows => 2, :ncols => 2} if nlines >= 6
  return {:nrows => 1, :ncols => 1} 

end

def get_colormode(cl)
  return 'bow' if cl.option_present('BOW')
  return 'cow' if cl.option_present('COW')
  return 'wob' if cl.option_present('WOB')
  return 'cob' if cl.option_present('COB')
  return 'cow'
end

def vf_cactvs(cl, tmpdir, smiles_file)
  nlines = lines_in_file(smiles_file)
  rc = rows_and_cols(cl, nlines)

  cmds = []
  cmds << 'csps'
  cmds << '-space 1'
  cmds << '-postprocess 0'
  cmds << '-ncols ' << rc[:ncols]
  cmds << '-nrows ' << rc[:nrows]
  cmds << '-colormode ' << get_colormode(cl)
  if cl.option_present('s')
    cmds << '-matchmode distinct -matchsmiles ' << cl.value('s')
  end
  cmds << "-orientation landscape" if cl.option_present('WIDE')
  cmds << '-title ' << cl.value('title') if cl.option_present('title')
  cmds << '--symbolfontsize ' << cl.value('sfs') if cl.option_present('sfs')
  cmds << cl.value('CACTVS') if cl.option_present('CACTVS')
  cmds << smiles_file
  if cl.option_present('P')
    psfile = cl.value('P')
  else
    psfile = File.join(tmpdir, "vftmp#{Process.pid}.ps")
  end
  cmds << " > #{psfile}"
  cmd = cmds.join(' ')
  $stderr << cmd << "\n" if cl.option_present('v')
  system(cmd)
  if ! File.size?(psfile)
    $stderr << "#{cmd} did not produce psfile\n"
    return false
  end

  if cl.option_present('pdf')
    system("ps2pdf #{psfile}")
  end

  system("evince #{psfile}") unless cl.option_present('K')

  File.unlink(psfile) unless cl.option_present('P')
  true
end

def vf(cl, tmpdir, smiles)
  if cl.option_present('rdkit') || cl.option_present('RDKIT')
    return vf_rdkit(cl, smiles)
  end
  if cl.option_present('mview') || cl.option_present('MVIEW') || cl.option_present('marvin')
    return vf_mview(cl, smiles)
  end
  vf_cactvs(cl, tmpdir, smiles)
end

# Call numbered_smiles on `input` with output to `destination`.
def make_numbered_smiles(n_option, verbose, input, destination)
  cmd = "numbered_smiles.sh -S - -n #{n_option} #{input} > #{destination}"
  $stderr << "Executing #{cmd}\n" if verbose
  system(cmd)
end

def smiles_to_tmp(smiles, fname)
  output = File.open(fname, 'w') do |output|
    smiles.each do |smi|
      output << "#{smi}\n"
    end
  end
end

def lsn_to_tmp(lsns, destination)
  opts = lsns.map{ |lsn| "-K #{lsn}" }
  cmd = 'selimsteg.sh ' << opts.join(' ') << " > #{destination}"
  raise "#{cmd} failed" unless system(cmd)
  $stderr << "Just did #{cmd}\n"
  true
end

def files_of_lsns_to_tmp(fnames, desination)
  cmd = 'selimsteg.sh ' << fnames.join(' ') << " >#{desination}"
  raise "#{cmd} failes" unless system(cmd)
  true
end

def smiles_from_cmdline(cl, files, destination)
  files = files_that_exist(ARGV)
  selector = File.join(File.dirname(__FILE__), 'vf_record_selector.rb')
  raise "Where is @{selector}" unless File.size?(selector)
  if cl.option_present('a')
    nsel = nil
  elsif cl.option_present('n')
    nsel = cl.value('n')
  else
    nsel = 30
  end

  cmd = []
  cmd << 'ruby'
  cmd << selector
  cmd << '-v' if cl.option_present('v')
  cmd << '-n' << nsel if nsel
  cmd << '-a' if cl.option_present('a')
  cmd << '-p' << cl.value('p') if cl.option_present('p')
  cmd << '-z' << cl.value('z') if cl.option_present('z')
  cmd << '-Z' if cl.option_present('Z')
  cmd << '-m' << cl.value('m') if cl.option_present('m')
  cmd << '-j' << cl.value('j') if cl.option_present('j')
  cmd << '-col' << cl.value('col') if cl.option_present('col')
  cmd << '-e' << cl.value('e') if cl.option_present('e')
  cmd << '-slen' << cl.value('slen') if cl.option_present('slen')
  cmd << '-SLEN' << cl.value('SLEN') if cl.option_present('SLEN')
  cmd << '-f' << cl.value('f') if cl.option_present('f')
  cmd << '-grep' << cl.value('grep').to_s if cl.option_present('grep')
  cmd << '-grepv' << cl.value('grepv') if cl.option_present('grepv')

  command = cmd.join(' ') << ' ' << files.join(' ')  << " > #{destination}"
  execute_command(cl, command)
end

def main
  # copy the first two lines from vf_record_selector
  cl = IWCmdline.new('-v-a-n=ipos-p=float-z=ipos-Z-m=ipos-j=ipos-col=s' \
                     '-e=ipos-slen=ipos-SLEN=ipos-f=ipos-grep=s-grepv=s' \
                     '-title=s' \
                     '-c=ipos-r=ipos' \
                     '-X=ipos-Y=ipos'  \
                     '-cactvs-CACTVS=close-sfs=ipos-WIDE-BOW-WOB-COB-COW' \
                     '-rdkit-RDKIT=close' \
                     '-mview-marvin-MVIEW=close-spreadsheet' \
                     '-V=xfile-G' \
                     '-k-keep-K' \
                     '-M=s-P=s-pdf' \
                     '-l=list-L=list-S=list' \
                     '-s=s' \
                     '-N=s' \
                     '-tdt-tdtsmi=s' \
                     '-tmp=sfile-expert')
  if cl.unrecognised_options_encountered
    $stderr << "unrecognised_options_encountered\n"
    usage(1)
  end

  expert = cl.option_present('expert')

  verbose = cl.option_present('v')

  if cl.option_present('S') && ARGV.size > 0
    $stderr << "Cannot specify both smiles -S and files\n"
    usage(expert)
  end

  tmpdir = if cl.option_present('tmp')
             cl.value('tmp')
           else
             "/tmp"
           end

  files_to_remove = Set.new

  if cl.option_present('M')
    tmpsmi = cl.value('M')
  else
    tmpsmi = File.join(tmpdir, "vftmp#{Process.pid}.smi")
    files_to_remove << tmpsmi
  end

  
  if cl.option_present('k') || cl.option_present('K')
    files_to_remove.delete(tmpsmi)
  end

  if cl.option_present('S')
    smiles_to_tmp(cl.values('S'), tmpsmi)
  elsif cl.option_present('l')
    lsn_to_tmp(cl.values('l'), tmpsmi)
  elsif cl.option_present('L')
    files_of_lsns_to_tmp(cl.values('L'), tmpsmi)
  elsif ARGV.empty?
    usage(expert)
  else
    smiles_from_cmdline(cl, ARGV, tmpsmi)
  end

  if cl.option_present('N')
    another_tmp_smiles_file = File.join(tmpdir, "vftmp#{Process.pid}n.smi")
    files_to_remove << another_tmp_smiles_file
    make_numbered_smiles(cl.value('N'), verbose, tmpsmi, another_tmp_smiles_file)
    tmpsmi = "#{another_tmp_smiles_file}"
  end

  vf(cl, tmpdir, tmpsmi)

  files_to_remove.each do |fname|
    File.unlink(fname)
  end

  sleep(1)
end

main
