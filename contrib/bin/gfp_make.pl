#!/usr/bin/env perl

use strict;

# General script for making composite fingerprints.
# This script began life in the mid 1990's, and has evolved over time.
# It has become successively larger, and less elegant, but remains in
# widespread use today. It is the starting point for most similarity computations.
# Note that parts of this are specific to the Lilly environment, and will take
# effort to make work outside there. But most functionality is available.
# LILLYMOL_HOME=/path/to//LillyMolPrivate perl ../contrib/script/perl/gfp_make.pl -v -STD file.smi
# works.

# The entire purpose of this script is to construct a pipelined command
# where each pipeline component adds one, or more, fingerprints to the output.

die "Environment Variable LILLYMOL_HOME is not defined!!!" unless defined $ENV{LILLYMOL_HOME};
my $ianhome = $ENV{LILLYMOL_HOME};

die "Where is the LILLYMOL_HOME directory '$ianhome'" unless (-d $ianhome);

my $hostname = `hostname`;
chomp $hostname;

my $uname = `uname`;
chomp $uname;

my @bindir;

unshift (@bindir, "$ianhome/bin/Linux");
unshift (@bindir, "$ianhome/contrib/bin");

# Flags for each kind of fingerprint possible

my $hbd = 0;
my $mk = 0;
my $mk2 = 0;
my $ncmk = 0;

my @ttust = ();

my $mpr = 0;
my $chg = 0;

my $iwfp = 0;
my $nciwfp = 0;
my @nciwfp = ();
my @iwfpust = ();

my $gc = 0;
my $pubchem = 0;

my @ecust = ();

my @mapust = ();

my $es = 0;
my $es_options = "";

my $dbf = 0;
my $shadow = 0;
my $need3d = 0;
my $jurs = 0;

my $tsubstructure_fingerprints_non_colliding = 0;
my $tsubstructure_bitrep = 0;

my $molecular_abstraction_options = "";
my $molecular_abstraction_bfile_contents = "";
my $fingerprint_substructure_options = "";
my $tsubstructure_fingerprint_options = "";
my $rand = 0;
my $random_fp_options = "-d 0.3";
my @EZ = ();
my $ezf = 0;
my $ezf_write = "";
my $ezf_read = "";
my $ct = 0;
my $descriptors = "";
my $rs = 0;
my $rze = 0;
my $pathd = 0;
my $pathdt = 0;
my $pathdc = 0;
my $pathdtc = 0;

my $dicer = 0;

my $tsub = 0;
my $mapO = "";
my $erg = 0;
my $fperg = 0;
my $ncerg = 0;

my $cats = 0;
my $cats_min_path = 0;
my $cats_max_path = 0;
my $catsp = 0;
my $catsp_min_path = 0;
my $catsp_max_path = 0;

my $clogp = 0;
my $clogd = 0;
my $clogpd = 0;
my $clogp_bit_replicates = 0;

my $marvin = 0;
my $marvin_string = "";

my $psa = 0;
my $psa_bit_replicates = 0;
my $abr = 0;
my $abr_bit_replicates = 0;
my $abp = 0;
my $abp_bit_replicates = 0;

my $temperature = 0;
my $natoms = 0;
my $nrings = 0;
my $aroma = 0;

my $d2f = 0;
my $d2f_string = "";

my $i2f = 0;
my $i2f_string = "";

my $w = 0;
my $w_specification = "";

my $atp = 0;
my $atp_options = "";
my $atr = 0;
my $atr_options = "";
my $ata = 0;
my $ata_options = "";

my $ring_fingerprint = 0;
# my $ring_fingerprint_options = "";

my $similarity_to = 0;
my @similarity_to_targets = ();
my $similarity_to_options = "";

my @model;

my @yoyo;

my $flatten_counted = 0;

# We can insert arbitrary text into the pipeline

my $user_insert_string = "";
my $user_insert = 0;

my $convert_to_unique_smiles = 0;

my $tmpdir = ".";
if (! -w $tmpdir)
{
  $tmpdir = "/tmp";
}

# Reactions for molecular transformations

my $reactions = "";

my $max_ec_shell_length = "";
my $min_ec_shell_length = "";

my $multiplicative_ec = 1;

my $mapmin = "";
my $mapmax = "";

my $input_qualifiers = "";

my $work_as_filter = 0;
my $work_as_tdt_filter = 0;

my $sp2fb = "";

my $translate_elements = "";

my $join_existing_file = "";

my $chirality_fingerprint_opts = "";

# Flush all fingerprint generators after each molecule.
my $flush = 0;

my $verbose = 0;

my $expert = 0;

sub usage
{
  print STDERR "Generates composite fingerprints, smiles > stdout\n";
  print STDERR " -HBD           FP Hydrogen bonding fingerprints DNU\n" if ($verbose);
  print STDERR " -STD           shortcut for -MK -MK2 -IW -MPR\n";
  print STDERR " -MK            FP MACCS keys\n";
  print STDERR " -MK2           FP MACCS keys level 2 - multiple feature determination\n";
  print STDERR " -NCMK          FP non-colliding counted MACCS keys\n" if ($expert);
  print STDERR " -MPR           FP molecular properties (natoms, nrings, ...)\n";
  print STDERR " -CHG           apply formal charges then fingerprint\n" if ($expert);
  print STDERR " -IW            FP IW fingerprints\n";
  print STDERR " -w <option> -w option(s) passed to iwfp (two -w's are needed)\n" if ($expert);
  print STDERR " -GHOSE         FP Ghose Crippen atom types\n" if ($expert);
  print STDERR " -EC            FP Jibo's Extended Connectivity fingerprints\n" if ($expert);
  print STDERR " -EC...         FP variants on EC fingerprints, C TT PP ...\n" if ($expert);
  print STDERR " -EC<r>:<T>     FP Extended Connectivity fingerprints, radius <r> type <T> (ZYHAPCR)\n";
  print STDERR " -EX...         FP Extended Connectivity fingerprints (V2), custom atom type DNU\n" if ($expert);
  print STDERR " -ecmax <num>   max shell length for -EC fingerprints\n" if ($expert);
  print STDERR " -ecmin <num>   min shell length for -EC fingerprints\n" if ($expert);
  print STDERR " -ecbig         generate larger, but more precise -EC fingerprints\n" if ($expert);
  print STDERR " -extec ... -extec extra options passed to all EC fingerprints\n" if ($expert);
# print STDERR " -extra ... -extra extra options passed to -EX\n" if ($expert);
# print STDERR " -ESA           Jibo's Estate fingerprints (Atom Estate)\n";
# print STDERR " -ESH           Jibo's Estate fingerprints (Hydrogen Estate)\n";
  print STDERR " -TT            FP Topological Torsion fingerprints DNU\n" if ($expert);
  print STDERR " -TT...         FP Topological Torsion fingerprints (variants) DNU\n" if ($expert);
  print STDERR " -MAPZ          FP Merck Atom Pairs (atomic number atom types)\n" if ($expert);
  print STDERR " -MAP           FP Merck Atom Pairs (custom atom types)\n";
  print STDERR " -MAP<r>:<T><n> FP Merck Atom Pairs radius <r> type <T> (ZYHAPCR), iter <n>\n";
  print STDERR " -mapmin <n>    min path length for MAP type fingerprints\n" if ($expert);
  print STDERR " -mapmax <n>    max path length for MAP type fingerprints\n" if ($expert);
  print STDERR " -extmap ... -extmap extra options passed to map\n" if ($expert);
  print STDERR " -GRF           FP fingerprint molecular graph - all bonds single\n" if ($expert);
  print STDERR " -SKE           FP fingerprint molecular skeleton - all bonds single, all atoms Carbon\n" if ($expert);
  print STDERR " -RING          FP largest ring/ring system fingerprints\n" if ($expert);
  print STDERR " -PD -PDT -PDC -PDTC  FP various iwpathd fingerprints\n" if ($expert);
# print STDERR " -dicer         dicer fragment fingerprints\n" if ($expert);
  print STDERR " -RZE           FP ring size fingerprints DNU\n" if ($expert);
  print STDERR " -RS            FP ring substitution fingerprints DNU\n" if ($expert);
  print STDERR " -rings         FP fingerprints of just the rings (all)\n" if ($expert);
  print STDERR " -ringsI        FP fingerprints of just the rings (all) labelled\n" if ($expert);
  print STDERR " -ARING         FP fingerprints of aromatic rings\n" if ($expert);
  print STDERR " -SIM <fname>   FP fingerprint  of similarity to molecules in <fname>\n" if ($expert);
  print STDERR " -simopt ... -simopt passed directly to similarity_to_fingerprint (def -T 0.5)\n" if ($expert);
# print STDERR " -nspch         FP fingerprints of spinach trimmed form\n" if ($expert);
  print STDERR " -SCF0          FP fingerprint of the molecular scaffold\n";
  print STDERR " -SPINACH       spinach fingerprint\n" if ($expert);
  print STDERR " -MABS ... -MABS FP options for molecular abstractions\n" if ($expert);
  print STDERR " -TSUB ...      FP query/queries for tsubstructure, S:file for smarts\n" if ($expert);
  print STDERR " -tsubnc        generate tsubstructure fingerprints as counted\n" if ($expert);
  print STDERR " -FPS ...       FP query/queries for fingerprint_substructure, S:file for smarts DNU\n" if ($expert);
  print STDERR " -fpsopt ... -fpsopt extra options for -FPS (suggest -x 3)\n" if ($expert);
  print STDERR " -ts ... -ts    options for tsubstructure\n" if ($expert);
  print STDERR " -W <dname>     include fingerprint of descriptor <dname>\n";
  print STDERR " -R <rxn>       isostere reaction DNU\n" if ($expert);
  print STDERR " -R FRED        transform to isosteric form based on Fred's rules DNU\n" if ($expert);
  print STDERR " -T def         translate elements default is I=Cl, Br=Cl. Beware, smiles in GFP changed!!\n" if ($expert);
  print STDERR " -EZ <smarts>   FP query for EZ fingerprints: *-*=*-* DNU\n" if ($expert);
  print STDERR " -CT            FP fingerprint all cis-trans bonds\n" if ($expert);
  print STDERR " -DSC <prog>    FP descriptor programme to insert descriptors DNU\n" if ($expert);
  print STDERR " -CLOGP         FP Biobyte clogp and clogd as fingerprints\n" if ($expert);
  print STDERR " -MCLP, -MPKA, -MPKB, Marvin clogp, pKa, pKb - use -MALL for all\n" if ($expert);
  print STDERR " -MVD           FP Marvin logd 7.4\n" if ($expert);
  print STDERR " -PSA           FP Novartis Polar Surface Area\n" if ($expert);
  print STDERR " -ABR           FP Abraham fingerprint\n" if ($expert);
  print STDERR " -ABP           FP Abraham and Platts fingerprint\n" if ($expert);
  print STDERR " -INS ... -INS  FP insert arbitrary pipelined commands, F:... to read from file\n" if ($expert);
  print STDERR " -D2F ... -D2F  FP use descriptors_to_fingerprint to insert bucketised descriptors (very flexible)\n" if ($expert);
  print STDERR " -d2f <fname,r,b> FP use descriptors to insert bucketised descriptors (replicates,buckets) - less flexible\n" if ($expert);
  print STDERR " -I2F ... -I2F  FP use descriptor_file_to_01_fingerprints to insert integer descriptors\n" if ($expert);
  print STDERR " -SHD           FP fingerprints based on shadow descriptors (3D)\n" if ($expert);
  print STDERR " -DBF           FP fingerprints based on distance between features descriptors (3D)\n" if ($expert);
  print STDERR " -JURS          FP fingerprints based on Jurs descriptors (3d)\n" if ($expert);
# print STDERR " -nostd         do NOT use '-g all' to standardise structures\n" if ($expert);
  print STDERR " -SP2FB ... -SP2FB convert sparse fingerprints to fixed binary (gfp_sparse_to_fixed.sh)\n" if ($expert);
  print STDERR " -fltc          flatten counted fingerprints to count = 1\n" if ($expert);
  print STDERR " -JOIN <fname>  use tdt_join to merge in the contents of an existing .gfp file\n" if ($expert);
  print STDERR " -CHIRAL ...    chirality fingerprint\n" if $expert;
  print STDERR " -bindir <dir>  search in <dir> for binaries\n" if ($expert);
  print STDERR " -flush         flush after each molecule\n" if ($expert);
  print STDERR " -i <qualifier> input qualifiers to the first pipe stage\n" if ($expert);
  print STDERR " -f             work as a filter. Input assumed to be smiles\n" if ($expert);
  print STDERR " -tdt           work as a filter. Input must       be a TDT\n" if ($expert);
  print STDERR " -expert        more options\n";
  print STDERR " -v             verbose output\n";

  exit ($_[0]);
}

sub display_subset_specification_stuff
{
  print STDERR "There are five levels of abstraction supported\n";
  print STDERR " 0      the skeleton of the subset, all atoms single, all bonds carbon\n";
  print STDERR " 1      the skeleton of the subset, all atoms single, all heteroatoms sulphur\n";
  print STDERR " 2      all aromatic atoms reduced to 'a' types\n";
  print STDERR " 3      all aromatic atoms reduced to 'c' if possible. Some may be 'o' or 'n'\n";
  print STDERR " 4      complete form - subset unchanged\n";
  print STDERR "\n";
  print STDERR "These numbers can be followed by the letter I which means generate fingerprints that are\n";
  print STDERR "sensitive to the substitution patterns - this is done with isotopic labels\n";

  exit ($_[0]);
}

# The atom typing specification is kind of strange. There are a number of atom types
# that are recognised natively. But others need to be specified as UST:type.
# This is inherently unstable

sub ust_needed_for_atype
{
  use strict 'vars';

  my $atype = $_[0];

  if ($atype =~ /^TT\d*$/ ||
      $atype =~ /^SB\d*$/ ||
      $atype =~ /^C\d*$/ ||
      $atype =~ /^SFX\d*$/ ||
      $atype =~ /^PP\d*$/ || $atype =~ /^PP\d*=\S+$/)
  {
    return "";
  }
  else
  {
    return "UST:";
  }
}

my $fingerprints_specified = 0;

# Feb 2000. Too many people making mistakes, make dash_g a default

my $dash_g = "-g all -l";     

my $aromatic_smiles = "-A I";   # make default, IL: not -A D?

my $extra_iwfp_options = "";
my $extra_ec_options = "";
my $extra_exc_options = "";
my $extra_map_options = "";

# Temperature will barf on multi-fragment molecules

my $temperature_dash_m = "";

my $cmd_fname = "";

my @files_to_be_deleted = ();
 
my $argptr = 0;
OPTION: while ($argptr < @ARGV)
{
  my $opt = $ARGV[$argptr++];

  if ($opt eq "-v")
  {
    $verbose++;
  }
  elsif ($opt eq "-expert")
  {
    $expert = 1;
  }
  elsif ($opt eq "-CHG")
  {
    $chg = 1;
    $fingerprints_specified++;
    $dash_g = "";
  }
  elsif ($opt eq "-HBD" || $opt eq "-HB" || $opt eq "-FPHB")
  {
    die "Cannot specify multiple HBD fingerprints" if ($hbd);
    $hbd = 1;
    $fingerprints_specified++;
  }
  elsif ($opt eq "-STD")
  {
    die "Cannot multiply specify standard fingerprints" if ($mpr || $mk || $mk2 || $iwfp);
    $mpr = 1;
    $mk = 1;
    $mk2 = 1;
    $iwfp = 1;
    $fingerprints_specified += 4;
  }
  elsif ($opt eq "-MK" || $opt eq "-FPMK")
  {
    die "Cannot specify multiple MK fingerprints" if ($mk);
    $mk = 1;
    $fingerprints_specified++;
  }
  elsif ($opt eq "-MK2")
  {
    $mk2 = 1;
    $fingerprints_specified++;
  }
  elsif ($opt eq "-NCMK")
  {
    $ncmk = 1;
    $fingerprints_specified++;
  }
  elsif ($opt eq "-TT")
  {
    push(@ttust, "-J NCTT");
    $fingerprints_specified++;
  }
  elsif ($opt =~ /^-TT:(\S+)/)
  {
    my $atype = $1;
    my $ust = ust_needed_for_atype($atype);
    push (@ttust, "-P ${ust}${atype} -J NCTT${atype}");
    $fingerprints_specified++;
  }
  elsif ($opt =~ /^-TT(\S+)/)
  {
    my $atype = $1;
    push(@ttust, "-P ${atype}");
    $fingerprints_specified++;
  }
  elsif ($opt eq "-MAPO")
  {
    my $fname = $ARGV[$argptr++];
    if (! -s $fname)
    {
      print STDERR "Missing or empty MAP subset file '${fname}'\n";
    }
    if (length($mapO))
    {
      print STDERR "Sorry, only one MAPO option allowed, see Ian\n";
      exit(3);
    }
    $fingerprints_specified++;
    $mapO = $fname;
  }
  elsif ($opt =~ /^-MAPC([0-9]+)$/)
  {
    my $rad = $1;
    push(@mapust, "-C ${rad} -P C -J NCMAP${rad}C");
    $fingerprints_specified++;
  }
  elsif ($opt eq "-MAP" || $opt eq "-MAPZ")
  {
    push(@mapust, "-P Z -J NCMAP");
    $fingerprints_specified++;
  }
  elsif ($opt eq "-FAP" || $opt eq "-FAPZ")
  {
    push(@mapust, "-P Z -J FPMAP -m 1024");
    $fingerprints_specified++;
  }
  elsif ($opt =~ /^-MAP(\d+)$/)
  {
    my $maxrad = $1;
    push(@mapust, "-C ${maxrad} -P UST:Z -J NCMAP${maxrad}Z");
    $fingerprints_specified++;
  }
  elsif ($opt =~ /^-FAP(\d+)$/)
  {
    my $maxrad = $1;
    push(@mapust, "-C ${maxrad} -P UST:Z -J FPMAP${maxrad}Z -m 1024");
    $fingerprints_specified++;
  }
  elsif ($opt eq "-MAPC")
  {
    push(@mapust, "-P C -J NCMAPC");
    $fingerprints_specified++;
  }
  elsif ($opt eq "-FAPC")
  {
    push(@mapust, "-P C -J FPMAPC -m 1024");
    $fingerprints_specified++;
  }
  elsif ($opt =~ /^-MAP(\d+)-(\d+)$/)
  {
    my $minrad = $1;
    my $maxrad = $2;
    push(@mapust, "-c ${minrad} -C ${maxrad} -P UST:Z -J NCMAP${minrad}${maxrad}Z");
    $fingerprints_specified++;
  }
  elsif ($opt =~ /^-MAP(\d+)-(\d+):(\S+)$/)
  {
    my $minrad = $1;
    my $maxrad = $2;
    my $atype = $3;
    my $ust = ust_needed_for_atype($atype);
    push(@mapust, "-c ${minrad} -C ${maxrad} -P ${ust}${atype} -J NCMAP${minrad}${maxrad}${atype}");
    $fingerprints_specified++;
  }
  elsif ($opt =~ /^-FAP(\d+)-(\d+):(\S+)$/)
  {
    my $minrad = $1;
    my $maxrad = $2;
    my $atype = $3;
    my $ust = ust_needed_for_atype($atype);
    push(@mapust, "-c ${minrad} -C ${maxrad} -P ${ust}${atype} -J FPMAP${minrad}${maxrad}${atype} -m 1024");
    $fingerprints_specified++;
  }
  elsif ($opt =~ /^-MAP(\d+)-(\d+):(\S+)(\d+)$/)
  {
    my $minrad = $1;
    my $maxrad = $2;
    my $atype = $3;
    my $ust = ust_needed_for_atype($atype);
    my $iter = $4;
    push(@mapust, "-c ${minrad} -C ${maxrad} -P ${ust}${atype}${iter} -J NCMAP${minrad}${maxrad}${atype}${iter}");
    $fingerprints_specified++;
  }
  elsif ($opt =~ /^-MAP(\d+):(\S+)(\d+)$/)
  {
    my $rad = $1;
    my $atype = $2;
    my $ust = ust_needed_for_atype($atype);
    my $iter = $3;
    push(@mapust, "-C ${rad} -P ${ust}${atype}${iter} -J NCMAP${rad}${atype}${iter}");
    $fingerprints_specified++;
  }
  elsif ($opt =~ /^-FAP(\d+):(\S+)(\d+)$/)
  {
    my $rad = $1;
    my $atype = $2;
    my $ust = ust_needed_for_atype($atype);
    my $iter = $3;
    push(@mapust, "-C ${rad} -P ${ust}${atype}${iter} -J FPMAP${rad}${atype}${iter} -m 1024");
    $fingerprints_specified++;
  }
  elsif ($opt =~ /^-MAP(\d+):(\S+)/)
  {
    my $rad = $1;
    my $atype = $2;
    my $ust = ust_needed_for_atype($atype);
    push(@mapust, "-C ${rad} -P ${ust}${atype} -J NCMAP${rad}${atype}");
    $fingerprints_specified++;
  }
  elsif ($opt =~ /^-FAP(\d+):(\S+)/)
  {
    my $rad = $1;
    my $atype = $2;
    my $ust = ust_needed_for_atype($atype);
    push(@mapust, "-C ${rad} -P ${ust}${atype} -J FPMAP${rad}${atype} -m 1024");
    $fingerprints_specified++;
  }
  elsif ($opt =~ /^-MAP:(\S+)/)
  {
    my $atype = $1;
    my $ust = ust_needed_for_atype($atype);
    push(@mapust, "-P ${ust}${atype} -J NCMAP${atype}");
    $fingerprints_specified++;
  }
  elsif ($opt =~ /^-FAP:(\S+)/)
  {
    my $atype = $1;
    my $ust = ust_needed_for_atype($atype);
    push(@mapust, "-P ${ust}${atype} -J FPMAP${atype} -m 1024");
    $fingerprints_specified++;
  }
  elsif ($opt =~ /^-MAP(\d+)(\S+)/)
  {
    my $rad = $1;
    my $atype = $2;
    my $ust = ust_needed_for_atype($atype);
    push(@mapust, "-C ${rad} -P ${ust}${atype} -J NCMAP${rad}${atype}");
    $fingerprints_specified++;
  }
  elsif ($opt =~ /^-FAP(\d+)(\S+)/)
  {
    my $rad = $1;
    my $atype = $2;
    my $ust = ust_needed_for_atype($atype);
    push(@mapust, "-C ${rad} -P ${ust}${atype} -J FPMAP${rad}${atype} -m 1024");
    $fingerprints_specified++;
  }
  elsif ($opt =~ /^-MAP(\S+)/)
  {
    my $atype = $1;
    push(@mapust, "-P ${atype} -J NCMAP${atype}");
    $fingerprints_specified++;
    print STDERR " Atype is ${atype}\n";
  }
  elsif ($opt =~ /^-FAP(\S+)/)
  {
    my $atype = $1;
    push(@mapust, "-P ${atype} -J FPMAP${atype} -m 1024");
    $fingerprints_specified++;
    print STDERR " Atype is ${atype}\n";
  }
  elsif ($opt eq "-mapmin")
  {
    $mapmin = $ARGV[$argptr++];
  }
  elsif ($opt eq "-mapmax")
  {
    $mapmax = $ARGV[$argptr++];
  }
  elsif ($opt =~ /-mapmax([0-9][0-9]*)$/)
  {
    $mapmax = $1;
  }
  elsif ($opt eq "-IW" || $opt eq "-FPIW" )
  {
    die "Cannot specify multiple IWFP fingerprints" if ($iwfp);
    $iwfp = 1;
    $fingerprints_specified++;
  }
  elsif ($opt =~ /^-IW:(\S+)/)
  {
    my $atype = $1;
    my $ust = ust_needed_for_atype($atype);
    push(@iwfpust, "-P ${ust}${atype} -J FPIW${atype}");
    $fingerprints_specified++;
  }
  elsif ($opt =~ /^-IW(\d+):(\S+)/)
  {
    my $max_path_len = $1;
    my $atype = $2;
    my $ust = ust_needed_for_atype($atype);
    push(@iwfpust, "-R ${max_path_len} -P ${ust}${atype} -J FPIW${max_path_len}${atype}");
    $fingerprints_specified++;
  }
  elsif ($opt eq "-NCIW")
  {
    die "Cannot specify multiple IWFP fingerprints" if ($iwfp || $nciwfp);
    $nciwfp = 1;
    $fingerprints_specified += 1;
  }
  elsif ($opt =~ /^-NCIW:(\S+)/)
  {
    my $atype = $1;
    my $ust = ust_needed_for_atype($atype);
    push(@nciwfp, " -P ${ust}${atype} -J NCIW${atype}");
    $fingerprints_specified += 1;
  }
  elsif ($opt =~ /^-NCIW(\d+)$/)
  {
    my $max_path_len = $1;
    push(@nciwfp, " -R ${max_path_len} -J NCIW${max_path_len}");
    $fingerprints_specified += 1;
  }
  elsif ($opt =~ /^-NCIW(\d+):(\S+)/)
  {
    my $max_path_len = $1;
    my $atype = $2;
    my $ust = ust_needed_for_atype($atype);
    push(@nciwfp, " -P ${ust}${atype} -R ${max_path_len} -J NCIW${max_path_len}${atype}");
    $fingerprints_specified += 1;
  }
  elsif ($opt eq "-RAND")
  {
    $rand = 1;
    $fingerprints_specified++;
  }
  elsif ($opt eq "-GC" || $opt eq "-GHOSE")
  {
    die "Cannot specify multiple GC fingerprints" if ($gc);
    $gc = 1;
    $fingerprints_specified++;
  }
  elsif ($opt eq "-PUBCHEM")
  {
     $pubchem = 1;
     $fingerprints_specified++;
  }
  elsif ($opt eq "-EC" || $opt eq "-ECZ")    # strange, no radius
  {
    my $tmp = "-P Z -J NCECZ";
    push(@ecust, $tmp);
    $fingerprints_specified++;
  }
  elsif ($opt eq "-FEC")
  {
    push(@ecust, "-P Z -J FPECZ");
    $fingerprints_specified++;
  }
  elsif ($opt =~ /^-EC(\d)$/)
  {
    my $rad = $1;
    my $tmp = "-R ${rad} -P Z -J NCEC${rad}Z";
    push(@ecust, $tmp);
    $fingerprints_specified++;
  }
  elsif ($opt =~ /^-FEC(\d)$/)
  {
    my $rad = $1;
    my $tmp = "-R ${rad} -P Z -J FPEC${rad}Z";
    push(@ecust, $tmp);
    $fingerprints_specified++;
  }
  elsif ($opt =~ /^-ECC([0-9]+)$/)
  {
    my $rad = $1;
    my $tmp = "-R ${rad} -P C -J NCEC${rad}C";
    push(@ecust, $tmp);
    $fingerprints_specified++;
  }
  elsif ($opt =~ /^-EC(\d+):(\S+)/)
  {
    my $rad = $1;
    my $atype = $2;
    my $ust = ust_needed_for_atype($atype);
    push (@ecust, "-R ${rad} -P ${ust}${atype} -J NCEC${rad}${atype}");
    $fingerprints_specified++;
  }
  elsif ($opt =~ /^-FEC(\d+):(\S+)/)
  {
    my $rad = $1;
    my $atype = $2;
    my $ust = ust_needed_for_atype($atype);
    push (@ecust, "-R ${rad} -P ${ust}${atype} -J FPEC${rad}${atype}");
    $fingerprints_specified++;
  }
  elsif ($opt =~ /^-EC:(\S+)/)
  {
    my $atype = $1;
    my $ust = ust_needed_for_atype($atype);
    my $tmp = "-P ${ust}${atype} -J NCEC${atype}";
    push(@ecust, $tmp);
    $fingerprints_specified++;
  }
  elsif ($opt =~ /^-FEC:(\S+)/)
  {
    my $atype = $1;
    my $ust = ust_needed_for_atype($atype);
    my $tmp = "-P ${ust}${atype} -J FPEC${atype}";
    push(@ecust, $tmp);
    $fingerprints_specified++;
  }
  elsif ($opt =~ /^-EC(\d)(\S+)/)
  {
    my $rad = $1;
    my $atype = $2;
    my $tmp = "-P ${atype} -R ${rad} -J NCEC${rad}${atype}";
    push(@ecust, $tmp);
    $fingerprints_specified++;
  }
  elsif ($opt =~ /^-EC(\S+)(\d)$/)     # historically we interpreted ECTT3 to be EC3TT
  {
    my $atype = $1;
    my $rad = $2;
    if ('1' ne $rad)     # then it looks like a radius
    {
      $rad = $2;
      $atype = $1;
    }
    push(@ecust, "-P ${atype} -R ${rad} -J NCEC${rad}${atype}");
    $fingerprints_specified++;
  }
  elsif ($opt =~ /^-EC(\S+)/)
  {
    my $atype = $1;
    my $tmp = "-P ${atype} -J NCEC${atype}";
    push(@ecust, $tmp);
    $fingerprints_specified++;
  }
  elsif ($opt eq "-ESA")
  {
    $es++;
    $es_options = "-J estate -J NCESTA";
    $fingerprints_specified++;
  }
  elsif ($opt eq "-ESH")
  {
    $es++;
    $es_options = "-J hestate -J NCESTH";
    $fingerprints_specified++;
  }
  elsif ($opt eq "-MPR")
  {
    die "Cannot specify multiple MPR fingerprints" if ($mpr == 1);
    $mpr = 1;
    $fingerprints_specified++;       # well, not really a fingerprint
  }
  elsif ($opt eq "-GRF")
  {
    $molecular_abstraction_options .= " -c -a 'allbonds(- FP)'";
    $fingerprints_specified++;
  }
  elsif ($opt eq "-SKE")
  {
    $molecular_abstraction_options .= " -c -a 'charge.allbonds(-).allatoms(C FP)'";
    $fingerprints_specified++;
  }
  elsif ($opt eq "-ts")
  {
    my $gotclose;    # do we find a closing -ts

    $gotclose = 0;

    GET_TS: while (1)
    {
      if ($ARGV[$argptr] eq "-ts")
      {
        $gotclose = 1;
        $argptr++;
        last GET_TS;
      }

      $tsubstructure_fingerprint_options .= " $ARGV[$argptr]";

      $argptr++;

      last GET_TS if ($argptr >= @ARGV);
    }

    die "The -TS option grouping must be closed" unless ($gotclose);
#   print STDERR "Options '${tsubstructure_fingerprint_options}'\n";
  }
  elsif ($opt eq "-FPS")
  {
    $fingerprint_substructure_options = $ARGV[$argptr++];
    $fingerprints_specified++;
  }
  elsif ($opt eq "-fpsopt")
  {
    my $gotclose;
    $gotclose = 0;

    GET_FPS: while (1)
    {
      if ($ARGV[$argptr] eq "-fpsopt")
      {
        $gotclose = 1;
        $argptr++;
        last GET_FPS;
      }

      $fingerprint_substructure_options .= " $ARGV[$argptr]";

      $argptr++;

      last GET_FPS if ($argptr >= @ARGV);
    }

    die "The -fpsopt option grouping must be closed" unless ($gotclose);
#   print STDERR "Options '${tsubstructure_fingerprint_options}'\n";
  }
  elsif ($opt =~ /^-NCTSUB(\d*)$/ || $opt =~ /^-TSUBNC(\d*)$/)
  {
    $tsubstructure_bitrep = $1 if (length($1) > 0);
    if (0 == length($tsubstructure_fingerprint_options))
    {
      $tsubstructure_fingerprint_options = '-u -q ' . $ARGV[$argptr++];
    }
    else
    {
      $tsubstructure_fingerprint_options .= ' -q ' . $ARGV[$argptr++];
    }
    if (0 == $tsub)
    {
      $fingerprints_specified++;
      $tsub = 1;
    }
    $tsubstructure_fingerprints_non_colliding = 1;
  }
  elsif ($opt =~ /^-TSUB(\d*)$/)
  {
    $tsubstructure_bitrep = $1 if (length($1) > 0);
    if (0 == length($tsubstructure_fingerprint_options))
    {
      $tsubstructure_fingerprint_options = '-u -q ' . $ARGV[$argptr++];
    }
    else
    {
      $tsubstructure_fingerprint_options .= ' -q ' . $ARGV[$argptr++];
    }
    if (0 == $tsub)
    {
      $fingerprints_specified++;
      $tsub = 1;
    }
  }
  elsif ($opt eq "-tsubnc")
  {
    $tsubstructure_fingerprints_non_colliding = 1;
  }
  elsif ($opt eq "-bitrep")
  {
    $tsubstructure_bitrep = $ARGV[$argptr++];
  }
  elsif ($opt eq "-MABS")
  {
    my $gotclose;    # do we find a closing -ABS

    $gotclose = 0;

    GET_MABS: while (1)
    {
      if ($ARGV[$argptr] eq "-MABS")
      {
        $gotclose = 1;
        $argptr++;
        last GET_MABS;
      }

      $molecular_abstraction_options .= " $ARGV[$argptr]";

      $argptr++;

      last GET_MABS if ($argptr >= @ARGV);
    }

    die "The -MABS option grouping must be closed" unless ($gotclose);

    $fingerprints_specified++;
  }
  elsif ($opt eq "-SCFCP")   # without cyclopropyl
  {
    $molecular_abstraction_options .= " -Y nbits=2048 -a 'rmatoms(\[/IWxR\]!\@C1\[CD2\]\[CD2\]1).scaffold(FP)'";
    $fingerprints_specified++;
  }
  elsif ($opt eq "-SPINACH")
  {
    $molecular_abstraction_options .= "-Y nbits=2048 -c -a 'spinach(AROM=Ar ALIPH=C CHAIN=C FP)'";
    $fingerprints_specified++;
  }
  elsif ($opt eq "-RING")
  {
    $molecular_abstraction_options .= " -a 'bigring(FP)'";
    $fingerprints_specified++;
  }
  elsif ($opt eq "-RINGI")
  {
    $molecular_abstraction_options .= " -a 'bigring(ISO FP)'";
    $fingerprints_specified++;
  }
  elsif ($opt eq "-rings" || $opt eq "-RINGS")
  {
    $molecular_abstraction_options .= " -a 'rings(FP)'";
    $fingerprints_specified++;
  }
  elsif ($opt eq "-ringsI" || $opt eq "-RINGSI")
  {
    $molecular_abstraction_options .= " -a 'rings(ISO FP)'";
    $fingerprints_specified++;
  }
  elsif ($opt eq "-ARING")
  {
    $ring_fingerprint = 1;
    $fingerprints_specified++;
  }
  elsif ($opt eq "-SIM")
  {
    $similarity_to = 1;
    push(@similarity_to_targets, $ARGV[$argptr++]);
    $fingerprints_specified++;
  }
  elsif ($opt eq "-simopt")
  {
    my $tmp = $ARGV[$argptr++];
    while ($tmp ne "-simopt" && $argptr < @ARGV)
    {
      $similarity_to_options .= " $tmp";
      $tmp = $ARGV[$argptr++];
    }
  }
  elsif ($opt eq "-SCF0" || $opt eq "-SCF")
  {
    $molecular_abstraction_options .= " -Y nbits=2048 -a 'scaffold(FP)'";
    $fingerprints_specified++;
  }
  elsif ($opt eq "-SCFI")
  {
    $molecular_abstraction_options .= " -Y ip=5 -Y nbits=2048 -a 'scaffold(ISO FP)'";
    $fingerprints_specified++;
  }
  elsif ($opt =~ /-SCF:(\S+)/)
  {
    if ($1 =~ /UST/)
    {
      $molecular_abstraction_bfile_contents .= "scaffold(AT=$1 FP.2048)\n";
    }
    else
    {
      my $atype = $1;
      my $ust = ust_needed_for_atype($atype);
      $molecular_abstraction_bfile_contents .= "scaffold(AT=${ust}$1 FP.2048)\n";
    }

    $fingerprints_specified++;
  }
  elsif ($opt eq "-RS")
  {
    $rs = 1;
    $fingerprints_specified++;
  }
  elsif ($opt eq "-RZE")
  {
    $rze = 1;
    $fingerprints_specified++;
  }
  elsif ($opt eq "-PD")
  {
    $pathd = 1;
    $fingerprints_specified++;
  }
  elsif ($opt eq "-PDT")
  {
    $pathdt = 1;
    $fingerprints_specified++;
  }
  elsif ($opt eq "-PDC")
  {
    $pathdc = 1;
    $fingerprints_specified++;
  }
  elsif ($opt eq "-PDTC")
  {
    $pathdtc = 1;
    $fingerprints_specified++;
  }
  elsif ($opt eq "-DICER" || $opt eq "-dicer")
  {
    $dicer = 1;
    $fingerprints_specified++;
  }
  elsif ($opt =~ /^-ERG$/i)
  {
    $erg = 1;
    $fingerprints_specified++;
  }
  elsif ($opt eq "-FPERG")
  {
    $fperg = 1;
    $fingerprints_specified++;
  }
  elsif ($opt eq "-NCERG")
  {
    $ncerg = 1;
    $fingerprints_specified++;
  }
  elsif ($opt =~ /^-CATS$/i)
  {
    $cats = 1;
    $fingerprints_specified++;
  }
  elsif ($opt =~ /^-CATSP$/i)
  {
    $catsp = 1;
    $fingerprints_specified++;
  }
  elsif ($opt =~ /-CATS(\d+)$/)
  {
    $cats = 1;
    $cats_max_path = $1;
    $fingerprints_specified++;
  }
  elsif ($opt =~ /-CATSP(\d+)$/)
  {
    $catsp = 1;
    $catsp_max_path = $1;
    $fingerprints_specified++;
  }
  elsif ($opt =~ /-CATS(\d+)-(\d+)$/)
  {
    $cats = 1;
    $cats_min_path = $1;
    $cats_max_path = $2;
    $fingerprints_specified++;
  }
  elsif ($opt =~ /-CATSP(\d+)-(\d+)$/)
  {
    $catsp = 1;
    $catsp_min_path = $1;
    $catsp_max_path = $2;
    $fingerprints_specified++;
  }
  elsif ($opt eq "-g")
  {
    $dash_g = "-g all -l";
  }
  elsif ($opt eq "-nostd")
  {
    $dash_g = "";
  }
  elsif ($opt eq "-usmi")
  {
    $convert_to_unique_smiles = 1;
  }
  elsif ($opt eq "-A")
  {
    $aromatic_smiles = "-A I";
  }
  elsif ($opt eq "-w")
  {
    my $tmp = $ARGV[$argptr++];
    while ($tmp ne "-w" && $argptr < @ARGV)
    {
      $extra_iwfp_options .= " $tmp";
#     print STDERR "TMP '$tmp', extra_iwfp_options expanded to '$extra_iwfp_options'\n";
      $tmp = $ARGV[$argptr++];
    }
  }
  elsif ($opt eq "-extra")
  {
    my $tmp = $ARGV[$argptr++];
    while ($tmp ne "-extra" && $argptr < @ARGV)
    {
      $extra_exc_options .= " $tmp";
      $tmp = $ARGV[$argptr++];
    }
  }
  elsif ($opt eq "-extec")
  {
    my $tmp = $ARGV[$argptr++];
    while ($tmp ne "-extec" && $argptr < @ARGV)
    {
      $extra_ec_options .= " $tmp";
      $tmp = $ARGV[$argptr++];
    }
  }
  elsif ($opt eq "-extmap")
  {
    my $tmp = $ARGV[$argptr++];
    while ($tmp ne "-extmap" && $argptr < @ARGV)
    {
      $extra_map_options .= " $tmp";
      $tmp = $ARGV[$argptr++];
    }
  }
  elsif ($opt =~ /-ATP(\d)*:(\S+)/)
  {
    my $d = $1;
    my $atype = $2;
    my $ust = ust_needed_for_atype($atype);
    $atp_options = "-J NCATP";
    $atp_options .= "${d} -C $d" if $d && length($d);
    $atp_options .= " -P ${ust}$atype" if length($atype);    # should build into fingerprint
    $atp = 1;
    $fingerprints_specified++;
  }
  elsif ($opt =~ /^-ATP(\d+)/)
  {
    my $d = $1;
    $atp = 1;
    $atp_options = "-J NCATP -C $1";
    $fingerprints_specified++;
  }
  elsif ($opt eq "-ATP")
  {
    $atp_options = "-J NCATP";
    $atp = 1;
    $fingerprints_specified++;
  }
  elsif ($opt =~ /-ATR(\d)*:(\S+)/)
  {
    my $d = $1;
    my $atype = $2;
    my $ust = ust_needed_for_atype($atype);
    $atr_options = "-a -J NCATR";
    $atr_options .= "${d} -C $d" if $d && length($d);
    $atr_options .= " -P ${ust}$atype" if length($atype);    # should build into fingerprint
    $atr = 1;
    $fingerprints_specified++;
  }
  elsif ($opt eq "-ATR")
  {
    $atr_options = "-J NCATR -a";
    $atr = 1;
    $fingerprints_specified++;
  }
  elsif ($opt =~ /^-ATR(\d+)/)
  {
    my $d = $1;
    $atr = 1;
    $atr_options = "-J NCATR${d} -a -C ${d}";
    $fingerprints_specified++;
  }
  elsif ($opt =~ /-ATA(\d+)*:(\S+)/)
  {
    my $d = $1;
    my $atype = $2;
    my $ust = ust_needed_for_atype($atype);
    $ata_options = "-h -J NCATA";
    $ata_options .= "${d} -C $d" if $d && length($d);
    $ata_options .= " -P ${ust}$atype" if length($atype);    # should build into fingerprint
    $ata = 1;
    $fingerprints_specified++;
  }
  elsif ($opt eq "-ATA")
  {
    $atr_options = "-J NCATA -h";
    $ata = 1;
    $fingerprints_specified++;
  }
  elsif ($opt =~ /^-ATA(\d+)/)
  {
    my $d = $1;
    $ata = 1;
    $ata_options = "-J NCATA${d} -C $1 -h";
    $fingerprints_specified++;
  }
  elsif ($opt eq "-i")
  {
    $input_qualifiers .= " -i " . $ARGV[$argptr++];
  }
  elsif ($opt eq "-f" || $opt eq "-pipe")
  {
    $work_as_filter = 1;
  }
  elsif ($opt eq "-fltc")
  {
    $flatten_counted = 1;
    $fingerprints_specified++;
  }
  elsif ($opt eq "-tdt")
  {
    $work_as_tdt_filter = 1;
  }
  elsif ($opt eq "-flush")
  {
    $flush = 1;
  }
  elsif ($opt eq "-multi")
  {
    $temperature_dash_m = "-m";
  }
  elsif ($opt eq "-ecmax")
  {
    $max_ec_shell_length = $ARGV[$argptr++];
  }
  elsif ($opt eq "-ecmin")
  {
    $min_ec_shell_length = $ARGV[$argptr++];
  }
  elsif ($opt =~ /^-ecmax([0-9][0-9]*)$/)
  {
    $max_ec_shell_length = $1;
  }
  elsif ($opt eq "-ecbig")
  {
    $multiplicative_ec = 0;
  }
  elsif ($opt eq "-R")
  {
    my $tmp = $ARGV[$argptr++];
    if ($tmp eq "FRED")
    {
      $reactions = gather_isostere_reactions ("$ianhome/data/isosteres")
    }
    else
    {
      die "Missing or empty reaction '$tmp'" unless (-s $tmp);
      $reactions .= " -R $tmp";
    }
  }
  elsif ($opt eq "-T")
  {
    my $tmp = $ARGV[$argptr++];
    if ("def" eq $tmp)
    {
      $translate_elements .= " -t I=Cl -t Br=Cl";
    }
    else
    {
      $translate_elements .= " -t ${tmp}";
    }
  }
  elsif ($opt eq "-EZ")
  {
    push(@EZ, $ARGV[$argptr++]);
  }
  elsif ($opt eq "-ezf")
  {
    $ezf = 1;
    $fingerprints_specified++;
  }
  elsif ($opt eq "-ezf_write")
  {
    $ezf = 1;
    $ezf_write = $ARGV[$argptr++];
  }
  elsif ($opt eq "-ezf_read")
  {
    $ezf = 1;
    $ezf_read = $ARGV[$argptr++];
    die "Missing or empty ezf_read file '$ezf_read'" unless (-s $ezf_read);
  }
  elsif ($opt eq "-CT")
  {
    $ct = 1;
    $fingerprints_specified++;
  }
  elsif ($opt eq "-DSC")
  {
    $descriptors = $ARGV[$argptr++];
  }
  elsif ($opt eq "-CLOGP")
  {
    $fingerprints_specified++ if (0 == $clogp);
    $clogp = 1;
  }
  elsif ($opt eq "-CLOGD")
  {
    $fingerprints_specified++ if (0 == $clogp);
    $clogd = 1;
    $clogp = 1;
  }
  elsif ($opt eq "-CLOGPD")
  {
    $fingerprints_specified++ if (0 == $clogp);
    $clogpd = 1;
    $clogp = 1;
  }
  elsif ($opt =~ /-CLOGPD(\d+)$/)
  {
    $clogp_bit_replicates = $1;
    $clogpd = 1;
    $clogp = 1;
    $fingerprints_specified++ if (0 == $clogp);
  }
  elsif ($opt =~ /-CLOGP(\d+)$/)
  {
    $clogp = 1;
    $clogp_bit_replicates = $1;
    $fingerprints_specified++;
  }
  elsif ($opt eq "-MCLP")
  {
    $marvin_string .= ' -O clogp';
    $fingerprints_specified++ if (0 == $marvin);
    $marvin = 1;
  }
  elsif ($opt eq "-MPKA")
  {
    $marvin_string .= ' -O apK';
    $fingerprints_specified++ if (0 == $marvin);
    $marvin = 1;
  }
  elsif ($opt eq "-MPKB")
  {
    $marvin_string .= ' -O bpK';
    $fingerprints_specified++ if (0 == $marvin);
    $marvin = 1;
  }
  elsif ($opt eq "-MVD")
  {
    $marvin_string .= ' -O clogd';
    $fingerprints_specified++ if (0 == $marvin);
    $marvin = 1;
  }
  elsif ($opt eq "-MVALL")
  {
    $marvin_string .= ' -O clogp -O apK -O bpK';
    $fingerprints_specified++ if (0 == $marvin);
    $marvin = 1;
  }
  elsif ($opt eq "-NATOMS")
  {
    $natoms = 1;
    $fingerprints_specified++ if (0 == $temperature);
    $temperature++;
  }
  elsif ($opt eq "-NRINGS")
  {
    $nrings = 1;
    $fingerprints_specified++ if (0 == $temperature);
    $temperature++;
  }
  elsif ($opt eq "-AROMA")
  {
    $aroma = 1;
    $fingerprints_specified++ if (0 == $temperature);
    $temperature++;
  }
  elsif ($opt eq "-PSA")
  {
    $psa = 1;
    $fingerprints_specified++;
  }
  elsif ($opt =~ /^-PSA(\d+)/)
  {
    $psa = 1;
    $psa_bit_replicates = $1;
    $fingerprints_specified++;
  }
  elsif ($opt eq "-W")
  {
    $w++;
    $fingerprints_specified++ if (1 == $w);
    $w_specification .= " -G " . $ARGV[$argptr++];
  }
  elsif ($opt eq "-MODEL" || $opt eq "-MDIR" || $opt eq "-mdir")
  {
    push(@model, $ARGV[$argptr++]);
    $fingerprints_specified++; # if (1 == @model);
  }
  elsif ($opt eq "-YOYO")
  {
    push(@yoyo, $ARGV[$argptr++]);
    $fingerprints_specified++;
  }
  elsif ($opt eq "-ABR")
  {
    $abr = 1;
    $fingerprints_specified++;
  }
  elsif ($opt =~ /^-ABR(\d+)/)
  {
    $abr = 1;
    $fingerprints_specified++;
    $abr_bit_replicates = $1;
  }
  elsif ($opt eq "-ABP")
  {
    $abp = 1;
    $fingerprints_specified++;
  }
  elsif ($opt =~ /^-ABP(\d+)/)
  {
    $abp = 1;
    $fingerprints_specified++;
    $abp_bit_replicates = $1;
  }
  elsif ($opt eq "-INS")
  {
    $user_insert_string .= "|";
    my $tmp = $ARGV[$argptr++];
    while ($tmp ne "-INS" && $argptr < @ARGV)
    {
      $user_insert_string .= " $tmp";
      $tmp = $ARGV[$argptr++];
    }

    if (0 == $user_insert)
    {
      $user_insert = 1;
      $fingerprints_specified++;
    }
  }
  elsif ($opt eq "-ins")
  {
    my $tmp = $ARGV[$argptr++];
    if ( ! -x $tmp)
    {
      $tmp = find_executable($tmp);
    }
    $user_insert_string .= "| ${tmp} -";

    if (0 == $user_insert)
    {
      $user_insert = 1;
      $fingerprints_specified++;
    }
  }
  elsif ($opt eq "-D2F")
  {
    my $descriptors_to_fingerprint = find_executable('descriptors_to_fingerprint');

    $d2f_string .= "|${descriptors_to_fingerprint} -f ";
    my $tmp = $ARGV[$argptr++];
    while ($tmp ne "-D2F" && $argptr < @ARGV)
    {
      $d2f_string .= " $tmp";
      $tmp = $ARGV[$argptr++];
    }

    if (0 == $d2f)
    {
      $d2f = 1;
      $fingerprints_specified++;
    }

    $d2f_string .= " -";
  }
  elsif ($opt eq "-d2f")
  {
    my $descriptors_to_fingerprint = find_executable('descriptors_to_fingerprint');

    $d2f_string .= "|${descriptors_to_fingerprint} -f";

    my $tmp = $ARGV[$argptr++];
    if ($tmp =~ /(\S+),(\d+),(\d+)$/)
    {
      $tmp = $1;
      $d2f_string .= " -r $2 -b $3";
    }
    elsif ($tmp =~ /(\S+),(\d+)$/)
    {
      $tmp = $1;
      $d2f_string .= " -r $2";
    }

    $d2f_string .= " -P ${tmp} -";

    if (0 == $d2f)
    {
      $d2f = 1;
      $fingerprints_specified++;
    }
  }
  elsif ($opt eq "-I2F")
  {
    my $descriptor_file_to_01_fingerprints = find_executable('descriptor_file_to_01_fingerprints');

    $i2f_string .= "|${descriptor_file_to_01_fingerprints} -f ";
    my $tmp = $ARGV[$argptr++];
    while ($tmp ne "-I2F" && $argptr < @ARGV)
    {
      $i2f_string .= " $tmp";
      $tmp = $ARGV[$argptr++];
    }

    if (0 == $i2f)
    {
      $i2f = 1;
      $fingerprints_specified++;
    }

    $i2f_string .= " -";
  }
  elsif ($opt eq "-DBF")
  {
    $dbf = 1;
    $need3d = 1;
    $fingerprints_specified++;
  }
  elsif ($opt eq "-SHD")
  {
    $shadow = 1;
    $need3d = 1;
    $fingerprints_specified++;
  }
  elsif ($opt eq "-JURS")
  {
    $jurs = 1;
    $need3d = 1;
    $fingerprints_specified++;
  }
  elsif ($opt eq "-tmpdir")
  {
    $tmpdir = $ARGV[$argptr++];
    die "Cannot write to tmpdir '$tmpdir'" unless (-w $tmpdir);
  }
  elsif ($opt eq "-SP2FB")
  {
    my $tmp = $ARGV[$argptr++];
    while ($tmp ne "-SP2FB" && $argptr < @ARGV)
    {
      $sp2fb .= " ${tmp}";
      $tmp = $ARGV[$argptr++];
    }
  }
  elsif ($opt eq '-')
  {
    $argptr--;
    last OPTION;
  }
  elsif ($opt eq "-bindir")
  {
    my $bindir = $ARGV[$argptr++];

    die "Missing or invalid binary directory '$bindir'" unless (-d $bindir);

    unshift(@bindir, $bindir);
  }
  elsif ($opt eq "-CMD")
  {
    $cmd_fname = $ARGV[$argptr++];
  }
  elsif ($opt eq "-JOIN")
  {
    $join_existing_file = $ARGV[$argptr++];
    die "Missing or empty file to join '$join_existing_file'" unless (-s $join_existing_file);
  }
  elsif ($opt eq "-CHIRAL")
  {
    $chirality_fingerprint_opts = $ARGV[$argptr++];
    $fingerprints_specified++;
  }
  elsif ($opt =~ /^-/)
  {
    die "Unrecognised option '$opt'";
  }
  else
  {
    $argptr--;
    last OPTION;
  }
}

$fingerprints_specified++ if length($reactions);

$fingerprints_specified++ if (@EZ);

#print STDERR "extra iwfp options '$extra_iwfp_options'\n" if ($verbose);

if ($argptr >= @ARGV)
{
  print STDERR "Insufficient arguments\n";
  usage(1);
}

# When flushing output, each generator may have a different set
# of argumements needed.

sub flush_args {
  use strict 'vars';

  # If flushing is not requested, return an empty string.
  if (! $flush) {
    return "";
  }

  my $caller = $_[0];
  print STDERR "From $caller\n";
  if ($caller eq "temperature") {
    return " -Y flush";
  }
  if ($caller eq "iwdescr") {
    return " -B flush";
  }
  if ($caller eq "svmfp_evaluate") {
    return " -Y flush";
  }
  if ($caller eq "assign_formal_charges") {
    return " -o flush -o smi";
  }
  if ($caller eq "iwfp") {
    return " -U flush";
  }
  if ($caller eq "tnass") {
    return " -Y flush";
  }
  if ($caller eq "pubchem_fingerprints") {
    return " -X flush";
  }

  return "";
}

my $dsctmp2;    # used later

if (length($descriptors))
{
  my $input_files;
  $input_files = "";
  my $i;
  for ($i = $argptr; $i < @ARGV; $i++)
  {
    $input_files .= " $ARGV[$i]"
  }

  my $dsctmp1 = "$tmpdir/gfpdsctmp1$$";

  my $cmd = "$descriptors $input_files > $dsctmp1";

  print STDERR "Executing '$cmd'\n" if ($verbose);

  system($cmd);

  die "'$cmd' failed" unless (-s $dsctmp1);

  my $just_columns_with_same_sign = find_executable("just_columns_with_same_sign");

  $dsctmp2 = "$tmpdir/gfpdsctmp2$$";

  $cmd = "$just_columns_with_same_sign -j -v -p $dsctmp1 > $dsctmp2";

  print STDERR "Executing '$cmd'\n" if ($verbose);

  system($cmd);

  die "'$cmd' failed" unless (-s $dsctmp2);

  unlink($dsctmp1);
}

$dash_g = "-g all" if (length($temperature_dash_m));

# If no fingerprints specified, take the default set.

if (0 == $fingerprints_specified)
{
  $mk = 1;
  $mk2 = 1;
  $iwfp = 1;

  if (length($temperature_dash_m))
  {
    $mpr = 1;
    $fingerprints_specified = 4;
  }
  else
  {
    $mpr = 1;
    $fingerprints_specified = 4;
  }
}

if ($mpr && $temperature)
{
  print STDERR "Cannot specify both -MPR and individual molecular properties\n";
  exit 2;
}

if (length($extra_iwfp_options) > 0 && 0 == $iwfp && 0 == length(@nciwfp))
{
  print STDERR "The -w option only makes sense with the -IWFP option\n";
  usage(32);
}

my $hbd_cmd_first;
my $hbd_cmd_pipe;

if ($hbd)
{
  my $tnass = find_executable ("tnass");

  my $charges = "$ENV{LILLYMOL_HOME}/queries/charges";
  my $hbonds = "$ENV{LILLYMOL_HOME}/queries/hbonds";
  my $hbondpatterns = "$ENV{LILLYMOL_HOME}/queries/hbondpatterns";

  die "Yipes, cannot access charges directory '$charges', see Ian" unless (-f $charges && -r $charges);
  die "Yipes, cannot access hbonds directory '$hbonds', see Ian" unless (-f $hbonds && -r $hbonds);
  die "Yipes, cannot access hbondpatterns directory '$hbondpatterns', see Ian" unless (-f $hbondpatterns && -r $hbondpatterns);

  my $hbd_cmd = "$tnass -A D -H a=F:${hbonds}/acceptor -H d=${hbonds}/donor.qry -H label -N F:${charges}/queries -q F:${hbondpatterns}/nass -J FPHB ";

  if ($work_as_filter)
  {
    $hbd_cmd_first = "$hbd_cmd $dash_g -F $aromatic_smiles $input_qualifiers FILE";
  }
  else
  {
    $hbd_cmd_first = "$hbd_cmd $dash_g $aromatic_smiles $input_qualifiers FILE";
  }
  $hbd_cmd_pipe  = "$hbd_cmd -F - ";
}

my $mk_cmd_first;
my $mk_cmd_pipe;

if ($mk || $mk2 || $ncmk)
{
  my $maccskeys =  find_executable("maccskeys");

  my $mk_cmd = "$maccskeys -E autocreate -A D -n ";

  $mk_cmd .= " -J FPMK" if ($mk);

  $mk_cmd .= " -J LEVEL2=FPMK2" if ($mk2);

  $mk_cmd .= " -J NC=NCMK" if ($ncmk);

  if ($work_as_filter)
  {
    $mk_cmd_first = "$mk_cmd $dash_g $aromatic_smiles -f $input_qualifiers FILE";
  }
  elsif ($work_as_tdt_filter)
  {
    # $mk_cmd_first = "$mk_cmd $dash_g -f $input_qualifiers FILE";
    $mk_cmd_first = "$mk_cmd $dash_g $aromatic_smiles -f FILE";
  }
  else
  {
    $mk_cmd_first = "$mk_cmd $dash_g $aromatic_smiles $input_qualifiers FILE";
  }

  $mk_cmd_pipe  = "$mk_cmd -f -";
}

my $temperature_cmd_first;
my $temperature_cmd_pipe;

if ($mpr)
{
  my $temperature = find_executable("temperature");

  my $mpr_cmd = "$temperature -A D -J MPR -E autocreate $temperature_dash_m" . flush_args("temperature");

  if ($work_as_filter)
  {
    $temperature_cmd_first = "$mpr_cmd -f $dash_g $aromatic_smiles $input_qualifiers FILE";
  }
  elsif ($work_as_tdt_filter)
  {
    $temperature_cmd_first = "$mpr_cmd -f $dash_g $aromatic_smiles FILE";
  }
  else
  {
    $temperature_cmd_first = "$mpr_cmd $dash_g $aromatic_smiles $input_qualifiers FILE";
  }

  $temperature_cmd_pipe  = "$mpr_cmd -f -";
}

if ($temperature)
{
  my $t = find_executable('temperature');

  my $b = "-d 4 -B ";
  $b .= 'n' if ($natoms);
  $b .= 'r' if ($nrings);
  $b .= 'a' if ($aroma);

  my $temperature_cmd = "${t} -A D -J NCMP -E autocreate ${temperature_dash_m} ${b}" . flush_args("temperature");

  $temperature_cmd_first = "$temperature_cmd $dash_g $aromatic_smiles $input_qualifiers FILE";
  $temperature_cmd_pipe  = "$temperature_cmd -f -";
  $mpr = 1;    # we can switch now
}

my $w_cmd_first;
my $w_cmd_pipe;

if ($w)
{
  my $t = find_executable("iwdescr.sh");

  my $iwdescr_cmd = "${t} ${w_specification} -B quiet" . flush_args("iwdescr");

  if ($work_as_filter)
  {
    $w_cmd_first = "${iwdescr_cmd} FILE";
  }
  elsif ($work_as_tdt_filter)
  {
    $w_cmd_first = "${iwdescr_cmd} -G FILTER -";
  }
  else
  {
    $w_cmd_first = "${iwdescr_cmd} FILE";
  }

  $w_cmd_pipe = "${iwdescr_cmd} -G FILTER -";
}

my $model_cmd_first;
my $model_cmd_pipe;

if (@model)
{
   my $svmfp_evaluate = find_executable('svmfp_evaluate.v2.sh');

   my $svmfp_evaluate_cmd = "${svmfp_evaluate}" . flush_args("svmfp_evaluate");

   $model_cmd_first = "${svmfp_evaluate_cmd}";
   $model_cmd_pipe = "${svmfp_evaluate_cmd} -J FILTER";
}

my $assign_formal_charges_first;

if ($chg)
{
  my $assign_formal_charges = find_executable("assign_formal_charges");

  $assign_formal_charges_first = "${assign_formal_charges} -f l -S - -o tdt FILE";
}

my $iwfp_cmd_first;
my $iwfp_cmd_pipe;

if ($iwfp)
{
  my $iwfp = find_executable("iwfp");

  my $iwfp_cmd = "$iwfp -A D -E autocreate -J FPIW $extra_iwfp_options" . flush_args("iwfp");

  if ($work_as_filter)
  {
    $iwfp_cmd_first = "$iwfp_cmd $dash_g $aromatic_smiles $input_qualifiers FILE";
  }
  elsif ($work_as_tdt_filter)
  {
    $iwfp_cmd_first = "$iwfp_cmd $dash_g $aromatic_smiles -f -";
  }
  else
  {
    $iwfp_cmd_first = "$iwfp_cmd $dash_g $aromatic_smiles $input_qualifiers FILE";
  }
  $iwfp_cmd_pipe  = "$iwfp_cmd -f -";
}

if ($nciwfp)
{
  my $iwfp = find_executable("iwfp");

  my $iwfp_cmd = "$iwfp -A D -E autocreate -J NCIW -c 200000 $extra_iwfp_options" . flush_args("iwfp");

  if ($work_as_filter)
  {
    $iwfp_cmd_first = "$iwfp_cmd $dash_g $aromatic_smiles $input_qualifiers FILE";
  }
  elsif ($work_as_tdt_filter)
  {
    $iwfp_cmd_first = "$iwfp_cmd $dash_g $aromatic_smiles -f -";
  }
  else
  {
    $iwfp_cmd_first = "$iwfp_cmd $dash_g $aromatic_smiles $input_qualifiers FILE";
  }
  $iwfp_cmd_pipe  = "$iwfp_cmd -f -";
}

my $nciwfp_cmd_first;
my $nciwfp_cmd_pipe;

if (@nciwfp)
{
  my $iwfp = find_executable("iwfp");

  my $iwfp_cmd = "${iwfp} -A D -E autocreate -c 200000" . flush_args("iwfp");

  if ($work_as_filter)
  {
    $nciwfp_cmd_first = "${iwfp_cmd} ${dash_g} ${aromatic_smiles} ${input_qualifiers} FILE";
  }
  elsif ($work_as_tdt_filter)
  {
    $nciwfp_cmd_first = "${iwfp_cmd} ${dash_g} ${aromatic_smiles} ${input_qualifiers} -f -";
  }
  else
  {
    $nciwfp_cmd_first = "${iwfp_cmd} ${dash_g} ${aromatic_smiles} ${input_qualifiers} FILE";
  }

  $nciwfp_cmd_pipe = "${iwfp_cmd} -f -";
}

if (@iwfpust)
{
  my $iwfp = find_executable('iwfp');

  my $iwfp_cmd = "$iwfp -A D -E autocreate" . flush_args("iwfp");

  if ($work_as_filter)
  {
    $iwfp_cmd_first = "$iwfp_cmd $dash_g $aromatic_smiles $input_qualifiers FILE";
  }
  elsif ($work_as_tdt_filter)
  {
    $iwfp_cmd_first = "$iwfp_cmd $dash_g $aromatic_smiles $input_qualifiers -f -";
  }
  else
  {
    $iwfp_cmd_first = "$iwfp_cmd $dash_g $aromatic_smiles $input_qualifiers FILE";
  }
  $iwfp_cmd_pipe  = "$iwfp_cmd -f -";
}

my $gc_cmd_first;
my $gc_cmd_pipe;

if ($gc)
{
  my $tnass = find_executable("tnass");

  my $queries = "$ianhome/lib/wildman_crippen.dat";
  die "Where are the GC queries '$queries'" unless (-s $queries && -r $queries);

  my $gc_cmd = "$tnass -A D -J NCGC -E autocreate -q S:$queries". flush_args("tnass");

  $gc_cmd_first = "$gc_cmd $dash_g $aromatic_smiles $input_qualifiers FILE";
  $gc_cmd_pipe  = "$gc_cmd -F -";
  
  print STDERR "Ignore the 'no available queries' message from tnass, it is harmless\n";
}

my $pubchem_fingerprints_first;
my $pubchem_fingerprints_pipe;

if ($pubchem)
{
  my $pubchem_fingerprints = find_executable('pubchem_fingerprints.sh');

  my $pubchem_fingerprints_cmd = "${pubchem_fingerprints} -E autocreate -l" . flush_args("pubchem_fingerprints");

  $pubchem_fingerprints_first = "${pubchem_fingerprints_cmd} ${dash_g} ${aromatic_smiles} ${input_qualifiers} FILE";
  $pubchem_fingerprints_pipe = "${pubchem_fingerprints_cmd} -f -";
}
my $ttust_cmd_first;
my $ttust_cmd_pipe;

if (@ttust)
{
  my $topotorsion_fingerprints = find_executable("topotorsion_fingerprints");

  my $tt_cmd = "${topotorsion_fingerprints} -A D -E autocreate";

  if ($work_as_filter)
  {
    $ttust_cmd_first = "${tt_cmd} ${dash_g} ${aromatic_smiles} ${input_qualifiers} FILE";
  }
  elsif ($work_as_tdt_filter)
  {
    $ttust_cmd_first = "${tt_cmd} ${dash_g} ${aromatic_smiles} ${input_qualifiers} -f FILE";
  }
  else
  {
    $ttust_cmd_first = "${tt_cmd} ${dash_g} ${aromatic_smiles} ${input_qualifiers} FILE";
  }

  $ttust_cmd_pipe = "${tt_cmd} -f -";
}

my $map_cmd_first;
my $map_cmd_pipe;

if (@mapust)
{
  my $map = find_executable("extended_atom_pairs");

  my $map_cmd = "$map -A D -E autocreate ${extra_map_options}";
  $map_cmd .= " -c $mapmin" if (length($mapmin));
  $map_cmd .= " -C $mapmax" if (length($mapmax));

  if ($work_as_filter)
  {
    $map_cmd_first = "$map_cmd $dash_g $aromatic_smiles $input_qualifiers FILE";
  }
  elsif ($work_as_tdt_filter)
  {
    $map_cmd_first = "$map_cmd $dash_g $aromatic_smiles $input_qualifiers -f FILE";
  }
  else
  {
    $map_cmd_first = "$map_cmd $dash_g $aromatic_smiles $input_qualifiers FILE";
  }
  $map_cmd_pipe  = "$map_cmd -f -";
}

my $mapO_cmd_first;
my $mapO_cmd_pipe;

if (length($mapO))
{
  my $map = find_executable("extended_atom_pairs");

  my $map_cmd = "$map -A D -E autocreate -J NCMAPO -O F:${mapO} ${extra_map_options}";

  if ($work_as_filter)
  {
    $mapO_cmd_first = "$map_cmd $dash_g $aromatic_smiles $input_qualifiers -f FILE";
  }
  else
  {
    $mapO_cmd_first = "$map_cmd $dash_g $aromatic_smiles $input_qualifiers FILE";
  }
  $mapO_cmd_pipe  = "$map_cmd -f -";
}

# All the EC fingerprints do not work properly running as TDT filters because
# of lack of standardisation. They are OK in the middle of a pipeline because
# something earlier will have standardised the molecule. They cannot be at
# the beginning of a tdt filter. For now, we just comment out those lines
# and let the computations fail. Next time we do an update, lift these
# restrictions

my $ec_cmd_first;
my $ec_cmd_pipe;

if (@ecust)
{
  my $extended_connectivity = find_executable("iwecfp");

  my $ec_cmd = "$extended_connectivity -E autocreate ";
  $ec_cmd .= "-R $max_ec_shell_length " if (length($max_ec_shell_length));
  $ec_cmd .= "-r $min_ec_shell_length" if (length($min_ec_shell_length));
  $ec_cmd .= $extra_ec_options if (length($extra_ec_options));

  $ec_cmd .= " -m" if ($multiplicative_ec);

  if ($work_as_filter)
  {
    $ec_cmd_first = "$ec_cmd $dash_g $aromatic_smiles $input_qualifiers FILE";
  }
  elsif ($work_as_tdt_filter)
  {
    $ec_cmd_first = "$ec_cmd $dash_g $aromatic_smiles -f $input_qualifiers FILE";
  }
  else
  {
    $ec_cmd_first = "$ec_cmd $dash_g $aromatic_smiles $input_qualifiers FILE";
  }

  $ec_cmd_pipe = "$ec_cmd -f -";
}

my $rand_cmd_first;
my $rand_cmd_pipe;

if ($rand)
{
  my $random_fingerprint = find_executable("random_fingerprint");

  if ($work_as_filter)
  {
    $rand_cmd_first = "$random_fingerprint $random_fp_options FILE";
  }
  elsif ($work_as_tdt_filter)
  {
    $rand_cmd_first = "$random_fingerprint -f $random_fp_options FILE";
  }
  else
  {
    $rand_cmd_first = "$random_fingerprint $random_fp_options FILE";
  }
  $rand_cmd_pipe = "$random_fingerprint $random_fp_options -f -";
}

my $es_cmd_first;
my $es_cmd_pipe;

if ($es > 1)
{
  print STDERR "Cannot have multiple -ES options\n";
  exit(1);
}

if ($es)
{
  my $estate = find_executable("jwestate");

  my $es_cmd = "$estate $es_options";

  if ($work_as_filter)
  {
    $es_cmd_first = "$es_cmd -f ${dash_g} $input_qualifiers FILE";
  }
  else
  {
    $es_cmd_first = "$es_cmd ${dash_g} $input_qualifiers FILE";
  }
  $es_cmd_pipe = "$es_cmd -f -";
}

my $molecular_abstraction_bfile_name;

if (length($molecular_abstraction_bfile_contents))
{
  $molecular_abstraction_bfile_name = "${tmpdir}/gfpmavbfile$$";
  open(my $bfile, ">$molecular_abstraction_bfile_name");
  die "Cannot open '${molecular_abstraction_bfile_name}'" unless $bfile;
  print $bfile $molecular_abstraction_bfile_contents;
  close($bfile);
  $molecular_abstraction_options .= " -B ${molecular_abstraction_bfile_name}";
  push(@files_to_be_deleted, $molecular_abstraction_bfile_name);
}

my $mabs_cmd_first;
my $mabs_cmd_pipe;

if (length($molecular_abstraction_options))
{
  my $molecular_abstraction = find_executable("molecular_abstraction");

  my $input_file_specification;
  if ($work_as_filter)
  {
    $input_file_specification = "$aromatic_smiles $input_qualifiers ";
  }
  elsif ($work_as_tdt_filter)
  {
    $input_file_specification = "$aromatic_smiles $input_qualifiers -f ";
  }
  else
  {
    $input_file_specification = "$aromatic_smiles $input_qualifiers ";
  }

  my $mabs_cmd_common = "$molecular_abstraction -E autocreate -A D -Y nbits=1024 -l $molecular_abstraction_options ";

  $mabs_cmd_first = "$mabs_cmd_common $dash_g $input_file_specification FILE";

  $mabs_cmd_pipe = "$mabs_cmd_common -f $input_file_specification -";
}

my $ring_fingerprint_first;
my $ring_fingerprint_pipe;

if ($ring_fingerprint)
{
  my $ring_fingerprint_exe = find_executable("ring_fingerprint");

  my $input_file_specification = "-p -d -E autocreate -l ${input_qualifiers} ";
  if ($work_as_tdt_filter)
  {
    $input_file_specification .= " -f"
  }
  elsif ($work_as_filter)
  {
    $input_file_specification .= " -i smi";
  }

  $ring_fingerprint_first = "${ring_fingerprint_exe} ${dash_g} ${input_file_specification} FILE";
  $ring_fingerprint_pipe  = "${ring_fingerprint_exe} ${input_file_specification} -f -";
}

my $similarity_to_first;
my $similarity_to_pipe;

if ($similarity_to)
{
  my $similarity_to_exe = find_executable("similarity_to");

  if ($similarity_to_options eq "")
  {
    $similarity_to_options = "-T 0.8";
  }

  my $input_file_specification = "${similarity_to_options} -E autocreate -l";
  if ($work_as_tdt_filter || $work_as_filter)
  {
    $input_file_specification .= " -f";
  }

  $similarity_to_first = "${similarity_to_exe} ${input_file_specification} ";     # note we do NOT insert FILE 
  $similarity_to_pipe = "${similarity_to_exe} ${input_file_specification} -f ";     # note we do NOT insert - 
}

my $fingerprint_substructure_cmd_first;
my $fingerprint_substructure_cmd_pipe;

if (length($fingerprint_substructure_options))
{
  my $fingerprint_substructure = find_executable("fingerprint_substructure");

  my $input_file_specification;
  if ($work_as_filter)
  {
    $input_file_specification = "$aromatic_smiles $input_qualifiers -";
  }
  elsif ($work_as_tdt_filter)
  {
    $input_file_specification = "$aromatic_smiles $input_qualifiers -f";
  }
  else
  {
    $input_file_specification = "$aromatic_smiles $input_qualifiers ";
  }

  my $fingerprint_substructure_cmd_common = "$fingerprint_substructure -z i -A D -Y nbits=1024 -q $fingerprint_substructure_options";

  $fingerprint_substructure_cmd_first = "$fingerprint_substructure_cmd_common $dash_g $input_file_specification FILE";
  $fingerprint_substructure_cmd_pipe  = "$fingerprint_substructure_cmd_common $input_file_specification -f -";
}

my $molecular_transformations_first;
my $molecular_transformations_pipe;

if(length ($reactions))
{
  my $molecular_transformations = find_executable("molecular_transformations");

  my $molecular_transformations_cmd = "$molecular_transformations -J FPISO -z i -z w -j -A D $reactions -F 'OUT=ISOSMI<'";

  if ($work_as_filter)
  {
    $molecular_transformations_first = "$molecular_transformations_cmd ${dash_g} -F 'IN=\$SMI<' $aromatic_smiles $input_qualifiers FILE ";
  }
  else
  {
    $molecular_transformations_first = "$molecular_transformations_cmd ${dash_g} $aromatic_smiles $input_qualifiers FILE ";
  }

  $molecular_transformations_pipe = "$molecular_transformations_cmd -F 'IN=\$SMI<' ";
}

my $ez_descriptor_first;
my $ez_descriptor_pipe;

if (@EZ)
{
  my $ez_descriptor = find_executable("ez_descriptor");
  my $ez_descriptor_cmd = "${ez_descriptor} -E autocreate -A D -J FPEZ -z i -z f";
  foreach $_ (@EZ)
  {
    $ez_descriptor_cmd .= " -s '$_'";
  }
  if ($work_as_filter)
  {
    $ez_descriptor_first = "$ez_descriptor_cmd -f FILE";
  }
  else
  {
    $ez_descriptor_first = "$ez_descriptor_cmd FILE";
  }

  $ez_descriptor_pipe = "$ez_descriptor_cmd -f -";
}

my $ezf_fingerprint_first;
my $ezf_fingerprint_pipe;

if ($ezf)
{
  my $ez_fingerprint = find_executable("ez_fingerprint");
  my $ez_fingerprint_cmd = "${ez_fingerprint} -v -r 10 -d 250 -E autocreate -l";

  $ez_fingerprint_cmd .= " -S ${ezf_write}"  if (length($ezf_write));
  $ez_fingerprint_cmd .= " -q ${ezf_read} -n"   if (length($ezf_read));

  if ($work_as_filter)
  {
    $ezf_fingerprint_first = "${ez_fingerprint_cmd} -f FILE";
  }
  else
  {
    $ezf_fingerprint_first = "${ez_fingerprint_cmd} FILE";
  }

  $ezf_fingerprint_pipe = "${ez_fingerprint_cmd} -f -";
}

my $ez_fingerprint_v2_first;
my $ez_fingerprint_v2_pipe;

if ($ct)
{
  my $ez_fingerprint_v2 = find_executable("ez_fingerprint_v2");
  my $ez_fingerprint_v2_cmd = "${ez_fingerprint_v2} -A D -P za -E autocreate ${dash_g} -l";

  if ($work_as_filter)
  {
    $ez_fingerprint_v2_first = "${ez_fingerprint_v2_cmd} -f FILE";
  }
  else
  {
    $ez_fingerprint_v2_first = "${ez_fingerprint_v2_cmd} FILE";
  }

  $ez_fingerprint_v2_pipe = "${ez_fingerprint_v2_cmd} -f -";
}

my $tsubstructure_first;
my $tsubstructure_pipe;

if (length($tsubstructure_fingerprint_options))
{
  my $tsubstructure = find_executable("tsubstructure");
  my $tsubstructure_command = "${tsubstructure} -E autocreate -A D -u ${tsubstructure_fingerprint_options} ";

  if ($tsubstructure_fingerprints_non_colliding)
  {
    $tsubstructure_command .= "-J NCTS ";
  }
  else
  {
    $tsubstructure_command .= "-J FPTS ";
  }

  if ($tsubstructure_bitrep)
  {
     $tsubstructure_command .= " -M bitrep=$tsubstructure_bitrep";
  }

  if ($work_as_filter)
  {
    $tsubstructure_first = "$tsubstructure_command FILE";
  }
  elsif ($work_as_tdt_filter)
  {
    $tsubstructure_first = "$tsubstructure_command -M tdt FILE";
  }
  else
  {
    $tsubstructure_first = "$tsubstructure_command FILE";
  }

  $tsubstructure_pipe = "$tsubstructure_command -M tdt -";
}

my $ring_substitution_first;
my $ring_substitution_pipe;

if ($rs)
{
  my $ring_substitution = find_executable("ring_substitution");
  my $ring_substitution_command = "$ring_substitution -A D -l -M full -M sfb -E autocreate";

  if ($work_as_filter)
  {
    $ring_substitution_first = "$ring_substitution_command ${dash_g} FILE";
  }
  elsif ($work_as_tdt_filter)
  {
    $ring_substitution_first = "$ring_substitution_command ${dash_g} -f FILE";
  }
  else
  {
    $ring_substitution_first = "$ring_substitution_command ${dash_g} FILE";
  }

  $ring_substitution_pipe = "$ring_substitution_command -f -";
}

my $ring_size_fingerprint_first;
my $ring_size_fingerprint_pipe;

if ($rze)
{
  my $ring_size_fingerprint = find_executable("ring_size_fingerprint");
  # my $ring_size_fingerprint_command = "$ring_size_fingerprint -E autocreate -x -e -t -r 8 ";
  # This is what Ian has, -x and -e is not compatible
  my $ring_size_fingerprint_command = "$ring_size_fingerprint -E autocreate -e -t ";

  if ($work_as_filter)
  {
    $ring_size_fingerprint_first = "$ring_size_fingerprint_command ${dash_g} FILE"
  }
  elsif ($work_as_tdt_filter)
  {
    $ring_size_fingerprint_first = "$ring_size_fingerprint_command ${dash_g} -f"
  }
  else
  {
    $ring_size_fingerprint_first = "$ring_size_fingerprint_command ${dash_g} FILE"
  }

  $ring_size_fingerprint_pipe = "$ring_size_fingerprint_command -f -";
}

my $iwpathd_first;
my $iwpathd_pipe;
my $iwpathdt_first;
my $iwpathdt_pipe;
my $iwpathdc_first;
my $iwpathdc_pipe;
my $iwpathdtc_first;
my $iwpathdtc_pipe;

if ($pathd || $pathdt || $pathdc || $pathdtc)
{
  my $iwpathd = find_executable("iwpathd");
  my $comm = "$iwpathd -A D ${dash_g} -E autocreate";
  my $iwpathd_command = "$comm -J NCPD";
  my $iwpathdt_command = "$comm -X tshape=5 -J NCPDT";
  my $iwpathdc_command = "$comm -P C -J NCPDC";
  my $iwpathdtc_command = "$comm -P C -X tshape=5 -J NCPDTC";

  if ($work_as_filter)
  {
    $iwpathd_first = "$iwpathd_command FILE";
    $iwpathdt_first = "$iwpathdt_command FILE";
    $iwpathdc_first = "$iwpathdc_command FILE";
    $iwpathdtc_first = "$iwpathdtc_command FILE";
  }
  elsif ($work_as_tdt_filter)
  {
    $iwpathd_first = "$iwpathd_command -f -";
    $iwpathdt_first = "$iwpathdt_command -f -";
    $iwpathdc_first = "$iwpathdc_command -f -";
    $iwpathdtc_first = "$iwpathdtc_command -f -";
  }
  else
  {
    $iwpathd_first = "$iwpathd_command FILE";
    $iwpathdt_first = "$iwpathdt_command FILE";
    $iwpathdc_first = "$iwpathdc_command FILE";
    $iwpathdtc_first = "$iwpathdtc_command FILE";
  }

  $iwpathd_pipe = "$iwpathd_command -f -";
  $iwpathdt_pipe = "$iwpathd_command -f -";
  $iwpathdc_pipe = "$iwpathdc_command -f -";
  $iwpathdtc_pipe = "$iwpathdc_command -f -";
}


my $dicer_first;
my $dicer_pipe;

if ($dicer)
{
  my $dicer_exe = find_executable("dicer");
  my $dicer_command = "$dicer_exe -X 500 -A D -J NCDC -k 3 -c -M 12"; 

  if ($work_as_filter)
  {
    $dicer_first = "$dicer_command FILE";
  }
  elsif ($work_as_tdt_filter)
  {
    $dicer_first = "$dicer_command -f FILE";
  }
  else
  {
    $dicer_first = "$dicer_command FILE";
  }

  $dicer_pipe = "$dicer_command -f -";
}

my $dbf_pipe;

if ($dbf)
{
  my $dbf_exe = find_executable("dbf");

  $dbf_pipe = "$dbf_exe -J NCDBF -";
}

my $shadow_pipe;

if ($shadow)
{
  my $shadow_exe = find_executable("tshadow");

  $shadow_pipe = "${shadow_exe} -J NCSHD -";
}

my $jurs_pipe;

if ($jurs)
{
  my $jurs_exe = find_executable("jursfp");

  $jurs_pipe = "${jurs_exe} -J NCJURS - 2> /dev/null";
}

my $corina_first;
my $corina_pipe;

if ($need3d)
{
  my $corina = find_executable("rcorina");
  my $corina_command = "$corina -x -r 4";

  if ($work_as_filter)
  {
    $corina_first = "${corina_command} -f FILE";
  }
  else
  {
    $corina_first = "${corina_command} -t FILE";
  }

  $corina_pipe = "${corina_command} -f -";
}

my $erg_cmd_first;
my $erg_cmd_pipe;

if ($erg)
{
  my $erg_exe = find_executable("ErG");

  my $erg_cmd = "${erg_exe} -J NCERG -j";

  if ($work_as_filter)
  {
    $erg_cmd_first = "${erg_cmd} FILE";
  }
  elsif ($work_as_tdt_filter)
  {
    $erg_cmd_first = "${erg_cmd} -f -";
  }
  else
  {
    $erg_cmd_first = "${erg_cmd} FILE";
  }

  $erg_cmd_pipe = "${erg_cmd} -f -";
}

my $ncerg_cmd_first;
my $ncerg_cmd_pipe;

if ($ncerg)
{
  my $erg_exe = find_executable("ErG");

  my $erg_cmd = "${erg_exe} -B all -M NCERG";

  if ($work_as_filter)
  {
    $ncerg_cmd_first = "${erg_cmd} FILE";
  }
  elsif ($work_as_tdt_filter)
  {
    $ncerg_cmd_first = "${erg_cmd} -f -";
  }
  else
  {
    $ncerg_cmd_first = "${erg_cmd} FILE";
  }

  $ncerg_cmd_pipe = "${erg_cmd} -f -";
}

my $fperg_cmd_first;
my $fperg_cmd_pipe;

if ($fperg)
{
  my $erg_exe = find_executable("ErG");

  my $erg_cmd = "${erg_exe} -B all -M FPERG";

  if ($work_as_filter)
  {
    $fperg_cmd_first = "${erg_cmd} FILE";
  }
  elsif ($work_as_tdt_filter)
  {
    $fperg_cmd_first = "${erg_cmd} -f -";
  }
  else
  {
    $fperg_cmd_first = "${erg_cmd} FILE";
  }

  $fperg_cmd_pipe = "${erg_cmd} -f -";
}


my $catsp_cmd_first;
my $catsp_cmd_pipe;

if ($catsp)
{
  my $catsp_exe = find_executable("jwcats");

  my $tag = "NCCATSP";

  my $catsp_cmd = "${catsp_exe} -E autocreate -p ";
  if ($catsp_min_path > 0)
  {
    $catsp_cmd .= " -z ${catsp_min_path}";
    $tag .= $catsp_min_path;
  }
  if ($catsp_max_path > 0)
  {
    $catsp_cmd .= " -m ${catsp_max_path}";
    $tag .= $catsp_max_path;
  }

  $catsp_cmd .= " -J ${tag}";

  if ($work_as_filter)
  {
    $catsp_cmd_first = "${catsp_cmd} ${dash_g} FILE";
  }
  elsif ($work_as_tdt_filter)
  {
    $catsp_cmd_first = "${catsp_cmd} ${dash_g} -f - ";
  }
  else
  {
    $catsp_cmd_first = "${catsp_cmd} ${dash_g} ${input_qualifiers} FILE";
  }

  $catsp_cmd_pipe = "${catsp_cmd} -f -";
}

my $cats_cmd_first;
my $cats_cmd_pipe;

if ($cats)
{
  my $cats_exe = find_executable("jwcats");

  my $tag = "NCCATS";
  my $cats_cmd = "${cats_exe} -E autocreate";
  if ($cats_min_path > 0)
  {
    $cats_cmd .= " -z ${cats_min_path}";
    $tag .= $cats_min_path;
  }
  if ($cats_max_path > 0)
  {
    $cats_cmd .= " -m ${cats_max_path}";
    $tag .= $cats_max_path;
  }

  $cats_cmd .= " -J ${tag}";

  if ($work_as_filter)
  {
    $cats_cmd_first = "${cats_cmd} ${dash_g} FILE";
  }
  elsif ($work_as_tdt_filter)
  {
    $cats_cmd_first = "${cats_cmd} ${dash_g} -f -";
  }
  else
  {
    $cats_cmd_first = "${cats_cmd} ${dash_g} ${input_qualifiers} FILE";
  }

  $cats_cmd_pipe = "${cats_cmd} -f -";
}

# With the non-first clogp invocations, one could argue that invoking fileconv a second time
# is not necessary. Probably...

my $clogp_cmd_first;
my $clogp_cmd_pipe;

if ($clogp)
{
  my $fileconv = find_executable("fileconv");
  my $clogp_exe = find_executable("bb_clogp-5.4_fp");
  my $clogp_cmd = "${clogp_exe}";
  my $fileconv_options = "-A I -E autocreate -f lod -g all -I 0 -V -S -";    # always invoke chemical standardisation

  $clogp_cmd .= " -p ${clogp_bit_replicates}" if ($clogp_bit_replicates > 0);

  $clogp_cmd .= " -k dd" if ($clogd);
  $clogp_cmd .= " -k pd" if ($clogpd);

  my $sed_cmd = "sed -e 's/^\$SMI<\\(..*\\)>/\\1/'";     # expose the smiles

  if ($work_as_tdt_filter)
  {
    $clogp_cmd_first = "${fileconv} -i tdt -i info -o tdt -o info ${fileconv_options} FILE | ${sed_cmd} | ${clogp_cmd}";
  }
  else
  {
    $clogp_cmd_first = "${fileconv} ${fileconv_options} FILE | ${clogp_cmd} ";
  }

  $clogp_cmd_pipe = "${fileconv} -i tdt -i info -o tdt -o info ${fileconv_options} - |${sed_cmd}| ${clogp_cmd}";

  if ($work_as_tdt_filter)
  {
    $clogp_cmd_first .= " -q";
    $clogp_cmd_pipe  .= " -q";
  }
}

my $marvin_cmd_first;

if ($marvin)
{
  my $marvin_exe = find_executable('iwmarvin.sh');
  my $m2gfp = find_executable('marvin2gfp.sh');

  if ($work_as_filter)
  {
    print STDERR "Marvin does not work as a filter\n";
    exit(2);
  }
  else
  {
    $marvin_cmd_first = "${marvin_exe} FILE | ${m2gfp} ${marvin_string} -";
  }
}

my $psa_cmd_first;
my $psa_cmd_pipe;

if ($psa)
{
  my $psa_exe = find_executable("psafp");
  my $psa_cmd = "${psa_exe}";

  if ($psa_bit_replicates > 0)
  {
    $psa_cmd .= " -p ${psa_bit_replicates} -J NCPSA${psa_bit_replicates}"
  }

  if ($work_as_filter)
  {
    $psa_cmd_first = "${psa_cmd} ${dash_g} ${aromatic_smiles} FILE";
  }
  elsif ($work_as_tdt_filter)
  {
    $psa_cmd_first = "${psa_cmd} ${dash_g} ${aromatic_smiles} -f -";
  }
  else
  {
    $psa_cmd_first = "${psa_cmd} ${dash_g} ${aromatic_smiles} FILE";
  }

  $psa_cmd_pipe = "${psa_cmd} -f -";
}

my $abr_cmd_first;
my $abr_cmd_pipe;

if ($abr)
{
  my $abr_exe = find_executable("abr.sh");
  my $abr_cmd = "${abr_exe}";

  if ($abr_bit_replicates > 0)
  {
    $abr_cmd .= " -p ${abr_bit_replicates} -J NCABR${abr_bit_replicates}"
  }
  else
  {
    $abr_cmd .= " -J NCABR";
  }

  if ($work_as_filter)
  {
    $abr_cmd_first = "${abr_cmd} ${dash_g} ${aromatic_smiles} FILE";
  }
  elsif ($work_as_tdt_filter)
  {
    $abr_cmd_first = "${abr_cmd} ${dash_g} ${aromatic_smiles} -f -";
  }
  else
  {
    $abr_cmd_first = "${abr_cmd} ${dash_g} ${aromatic_smiles} FILE";
  }

  $abr_cmd_pipe = "${abr_cmd} -f -";
}

my $abp_cmd_first;
my $abp_cmd_pipe;

if ($abp)
{
  my $abp_exe = find_executable("ap.sh");
  my $abp_cmd = "${abp_exe}";

  if ($abp_bit_replicates > 0)
  {
    $abp_cmd .= " -p ${abp_bit_replicates} -J NCABP${abp_bit_replicates}"
  }
  else
  {
    $abp_cmd .= " -J NCABP";
  }

  if ($work_as_filter)
  {
    $abp_cmd_first = "${abp_cmd} ${dash_g} ${aromatic_smiles} FILE";
  }
  elsif ($work_as_tdt_filter)
  {
    $abp_cmd_first = "${abp_cmd} ${dash_g} ${aromatic_smiles} -f -";
  }
  else
  {
    $abp_cmd_first = "${abp_cmd} ${dash_g} ${aromatic_smiles} FILE";
  }

  $abp_cmd_pipe = "${abp_cmd} -f -";
}

my $atp_first;
my $atp_cmd_pipe;

if ($atp)
{
  my $atp_exe = find_executable("atom_triples");
  my $atp_cmd = "${atp_exe} -A D -l ${atp_options}";

  if ($work_as_filter)
  {
    $atp_first = "${atp_cmd} ${dash_g} ${aromatic_smiles} FILE";
  }
  elsif ($work_as_tdt_filter)
  {
    $atp_first = "${atp_cmd} ${dash_g} ${aromatic_smiles} -f -";
  }
  else
  {
    $atp_first = "${atp_cmd} ${dash_g} ${aromatic_smiles} FILE"
  }

  $atp_cmd_pipe = "${atp_cmd} -f -";
}

my $atr_first;
my $atr_cmd_pipe;

if ($atr)
{
  my $atr_exe = find_executable("atom_triples");
  my $atr_cmd = "${atr_exe} -A D -l ${atr_options}";

  if ($work_as_filter)
  {
    $atr_first = "${atr_cmd} ${dash_g} ${aromatic_smiles} FILE";
  }
  elsif ($work_as_tdt_filter)
  {
    $atr_first = "${atr_cmd} ${dash_g} ${aromatic_smiles} -f -";
  }
  else
  {
    $atr_first = "${atr_cmd} ${dash_g} ${aromatic_smiles} FILE"
  }

  $atr_cmd_pipe = "${atr_cmd} -f -";
}

my $ata_first;
my $ata_cmd_pipe;

if ($ata)
{
  my $ata_exe = find_executable("atom_triples");
  my $ata_cmd = "${ata_exe} -A D -l ${ata_options}";

  if ($work_as_filter)
  {
    $ata_first = "${ata_cmd} ${dash_g} ${aromatic_smiles} FILE";
  }
  elsif ($work_as_tdt_filter)
  {
    $ata_first = "${ata_cmd} ${dash_g} ${aromatic_smiles} -f -";
  }
  else
  {
    $ata_first = "${ata_cmd} ${dash_g} ${aromatic_smiles} FILE"
  }

  $ata_cmd_pipe = "${ata_cmd} -f -";
}

my $flatten_counted_pipe;

if ($flatten_counted)
{
  if (1 == $fingerprints_specified)
  {
    print STDERR "In order to flatten fingerprints, there must be other fingerprints specified\n";
    usage(1);
  }

  my $flatten_counted_exe = find_executable("gfp_flatten_counted");

  $flatten_counted_pipe = "${flatten_counted_exe} -";
}

my $chirality_fingerprint_first;
my $chirality_fingerprint_pipe;

if (length($chirality_fingerprint_opts) > 0)
{
  my $chirality_fingerprint_exe = find_executable("chirality_fingerprint");
  if ("DEF" eq $chirality_fingerprint_opts)
  {
    $chirality_fingerprint_opts = "-e 1 -r 15 -y 1";
  }
  else
  {
    $chirality_fingerprint_opts =~ s/-/ -/g;
  }
  $chirality_fingerprint_exe .= " ${chirality_fingerprint_opts}";

  if ($work_as_filter)
  {
    $chirality_fingerprint_first = "${chirality_fingerprint_exe} ${dash_g} FILE";
  }
  elsif ($work_as_tdt_filter)
  {
    $chirality_fingerprint_first = "${chirality_fingerprint_exe} -f -";
  }
  else
  {
    $chirality_fingerprint_first = "${chirality_fingerprint_exe} ${dash_g} FILE";
  }

  $chirality_fingerprint_pipe = "${chirality_fingerprint_exe} -f -";
}

my $first = 1;

# Build the command once and then substitute the actual file name

my $cmd = "";

if ($convert_to_unique_smiles)
{
  my $fileconv = find_executable("fileconv");
  $cmd = "${fileconv} -c 1 -A D -E autocreate -o tdtnausmi -S -";
  $cmd .= " -g all -f lod" if (length($dash_g));    # fileconv does not recognise -l
  $cmd .= " FILE";
  $first = 0;
}
elsif (length($translate_elements))
{
  my $fileconv = find_executable("fileconv");
  $cmd = "${fileconv} -c 1 -A D -E autocreate -o tdt -S - ${translate_elements}";
  $cmd .= " -g all -f lod" if (length($dash_g));    # fileconv does not recognise -l
  $cmd .= " FILE";
  $first = 0;
}

# Built the command

while ($fingerprints_specified > 0)
{
  print STDERR "$fingerprints_specified fingerprints remaining\n" if ($verbose);
  my $file;

  if ($clogp)
  {
    $cmd = "${clogp_cmd_first}";
    $clogp = 0;
  }
  elsif ($chg)
  {
    $cmd = $assign_formal_charges_first;
    $chg = 0;
  }
  elsif ($mpr)
  {
    if ($first)
    {
      $cmd = $temperature_cmd_first;
    }
    else
    {
      $cmd .= "| $temperature_cmd_pipe";
    }

    $mpr = 0;
  }
  elsif ($mk || $mk2 || $ncmk)
  {
    if ($first)
    {
      $cmd = $mk_cmd_first;
    }
    else
    {
      $cmd .= "| $mk_cmd_pipe";
    }
    $fingerprints_specified -= ($mk + $mk2 + $ncmk) - 1;    # the number of fingerprints will be decremented at the bottom of the loop
    $mk = 0;
    $mk2 = 0;
    $ncmk = 0;
  }
  elsif ($iwfp)
  {
    if ($first)
    {
      $cmd = $iwfp_cmd_first;
    }
    else
    {
      $cmd .= "| $iwfp_cmd_pipe";
    }
    $iwfp = 0;
  }
  elsif ($nciwfp)
  {
    if ($first)
    {
      $cmd = $iwfp_cmd_first;
    }
    else
    {
      $cmd .= "| $iwfp_cmd_pipe";
    }
    $nciwfp = 0;
  }
  elsif (@nciwfp)
  {
    my $options = shift(@nciwfp);

    if ($first)
    {
      $cmd = "$nciwfp_cmd_first";
      $cmd =~ s/ / $options /;
#     print STDERR "Expanded to '${cmd}'\n";
    }
    else
    {
      my $tmp = $nciwfp_cmd_pipe;
      $tmp =~ s/ / $options /;
      $cmd .= "| $tmp";
    }
  }
  elsif (@iwfpust)
  {
    my $options = shift(@iwfpust);

#   print STDERR "Options are '${options}'\n";

    if ($first)
    {
      $cmd = "$iwfp_cmd_first";
      $cmd =~ s/ / $options /;
#     print STDERR "Expanded to '${cmd}'\n";
    }
    else
    {
      my $tmp = $iwfp_cmd_pipe;
      $tmp =~ s/ / $options /;
      $cmd .= "| $tmp";
    }
  }
  elsif ($hbd)
  {
    if ($first)
    {
      $cmd = $hbd_cmd_first;
    }
    else
    {
      $cmd .= "| $hbd_cmd_pipe";
    }
    $hbd = 0;
  }
  elsif ($w)
  {
    if ($first)
    {
      $cmd = $w_cmd_first;
    }
    else
    {
      $cmd .= "| $w_cmd_pipe";
    }
    $w = 0;
  }
  elsif ($gc)
  {
    if ($first)
    {
      $cmd = $gc_cmd_first;
    }
    else
    {
      $cmd .= "| $gc_cmd_pipe"
    }

    $gc = 0;
  }
  elsif ($pubchem)
  {
    if ($first)
    {
      $cmd = $pubchem_fingerprints_first;
    }
    else
    {
      $cmd .= "| $pubchem_fingerprints_pipe"
    }

    $pubchem = 0;
  }
  elsif ($rand)
  {
    if ($first)
    {
      $cmd = $rand_cmd_first;
    }
    else
    {
      $cmd .= "| $rand_cmd_pipe"
    }

    $rand = 0;
  }
  elsif (@ttust)
  {
    my $options = shift(@ttust);
    if ($first)
    {
      $cmd = "${ttust_cmd_first}";
      $cmd =~ s/ / $options /;
    }
    else
    {
      my $tmp = $ttust_cmd_pipe;
      $tmp =~ s/ / $options /;
      $cmd .= "| $tmp";
    }
  }
  elsif (@ecust)
  {
    my $options = shift(@ecust);

#   print STDERR "Options are '${options}'\n";

    if ($first)
    {
      $cmd = "$ec_cmd_first";
      $cmd =~ s/ / $options /;
#     print STDERR "Expanded to '${cmd}'\n";
    }
    else
    {
      my $tmp = $ec_cmd_pipe;
      $tmp =~ s/ / $options /;
      $cmd .= "| $tmp";
    }
  }
  elsif ($es)
  {
    if ($first)
    {
      $cmd = $es_cmd_first;
    }
    else
    {
      $cmd .= "| $es_cmd_pipe";
    }
    $es = 0;
  }
  elsif (@mapust)
  {
    my $options = shift(@mapust);

    if ($first)
    {
      $cmd = $map_cmd_first;
      $cmd =~ s/ / $options /;
    }
    else
    {
      my $tmp = $map_cmd_pipe;
      $tmp =~ s/ / $options /;
      $cmd .= "| $tmp";
    }
  }
  elsif (length ($molecular_abstraction_options))
  {
    if ($first)
    {
      $cmd = $mabs_cmd_first;
    }
    else
    {
      $cmd .= "| $mabs_cmd_pipe";
    }
    $molecular_abstraction_options = "";
  }
  elsif ($ring_fingerprint)
  {
    if ($first)
    {
      $cmd = $ring_fingerprint_first;
    }
    else
    {
      $cmd .= "| $ring_fingerprint_pipe";
    }
    $ring_fingerprint = 0;
  }
  elsif ($similarity_to)
  {
    if ($first)
    {
      my $tmp = shift(@similarity_to_targets);

      $cmd = "${similarity_to_first} -p ${tmp} FILE"
    }

    my $fname;
    foreach $fname (@similarity_to_targets) 
    {
      $cmd .= "|${similarity_to_pipe} -p ${fname} -"
    }

    $similarity_to = 0;
  }
  elsif ($psa)
  {
    if ($first)
    {
      $cmd = $psa_cmd_first;
    }
    else
    {
      $cmd .= "| $psa_cmd_pipe";
    }

    $psa = 0;
  }
  elsif ($marvin)
  {
    if ($first)
    {
      $cmd = $marvin_cmd_first;
      $marvin = 0;
    }
    else {
       print STDERR "Can only run Marvin FP on their own, NCMV* FP omitted\n";
    }
  }
  elsif ($abr)
  {
    if ($first)
    {
      $cmd = $abr_cmd_first;
    }
    else
    {
      $cmd .= "| $abr_cmd_pipe";
    }

    $abr = 0;
  }
  elsif ($abp)
  {
    if ($first)
    {
      $cmd = $abp_cmd_first;
    }
    else
    {
      $cmd .= "| $abp_cmd_pipe";
    }

    $abp = 0;
  }
  elsif (length($fingerprint_substructure_options))
  {
    if ($first)
    {
      $cmd = $fingerprint_substructure_cmd_first;
    }
    else
    {
      $cmd .= "| $fingerprint_substructure_cmd_pipe";
    }
    $fingerprint_substructure_options = "";
  }
  elsif (length($reactions))
  {
    if ($first)
    {
      $cmd = "$molecular_transformations_first";
    }
    else
    {
      $cmd .= "| $molecular_transformations_pipe";
    }
    $reactions = "";
  }
  elsif (@EZ)
  {
    if ($first)
    {
      $cmd = "$ez_descriptor_first";
    }
    else
    {
      $cmd .= "| $ez_descriptor_pipe";
    }
  }
  elsif ($ezf)
  {
    if ($first)
    {
      $cmd = "${ezf_fingerprint_first}";
    }
    else
    {
      $cmd .= "| ${ezf_fingerprint_pipe}";
    }
  }
  elsif ($ct)
  {
    if ($first)
    {
      $cmd = "${ez_fingerprint_v2_first}";
    }
    else
    {
      $cmd .= "| ${ez_fingerprint_v2_pipe}";
    }
    $ct = 0;
  }
  elsif ($tsubstructure_fingerprint_options)
  {
    if ($first)
    {
      $cmd = "$tsubstructure_first";
    }
    else
    {
      $cmd .= "|$tsubstructure_pipe";
    }
    $tsubstructure_fingerprint_options = "";
  }
  elsif ($rs)
  {
    if ($first)
    {
      $cmd = "$ring_substitution_first";
    }
    else
    {
      $cmd .= "|$ring_substitution_pipe";
    }

    $rs = 0;
  }
  elsif ($rze)
  {
    if ($first)
    {
      $cmd = "$ring_size_fingerprint_first";
    }
    else
    {
      $cmd .= "|$ring_size_fingerprint_pipe";
    }

    $rze = 0;
  }
  elsif ($pathd)
  {
    if ($first)
    {
      $cmd = "$iwpathd_first";
    }
    else
    {
      $cmd .= "|$iwpathd_pipe";
    }

    $pathd = 0;
  }
  elsif ($pathdt)
  {
    if ($first)
    {
      $cmd = "$iwpathdt_first";
    }
    else
    {
      $cmd .= "|$iwpathdt_pipe";
    }

    $pathdt = 0;
  }
  elsif ($pathdc)
  {
    if ($first)
    {
      $cmd = "$iwpathdc_first";
    }
    else
    {
      $cmd .= "|$iwpathdc_pipe";
    }

    $pathdc = 0;
  }
  elsif ($pathdtc)
  {
    if ($first)
    {
      $cmd = "$iwpathdtc_first";
    }
    else
    {
      $cmd .= "|$iwpathdtc_pipe";
    }

    $pathdtc = 0;
  }
  elsif ($dicer)
  {
    if ($first)
    {
      $cmd = "$dicer_first";
    }
    else
    {
      $cmd .= "|$dicer_pipe";
    }

    $dicer = 0;
  }
  elsif ($shadow)
  {
    if ($first)
    {
      if ($need3d)
      {
        $cmd = "${corina_first} | $shadow_pipe ";
      }
    }
    else
    {
      if ($need3d)
      {
         $cmd .= "|${corina_pipe} | ${shadow_pipe} ";
      }
      else
      {
        $cmd .= "|${shadow_pipe}";
      }
    }

    $shadow = 0;
    $need3d = 0;
  }
  elsif ($dbf)
  {
    if ($first)
    {
      $cmd = "${corina_first} | $dbf_pipe ";
    }
    else
    {
      if ($need3d)
      {
        $cmd .= "|${corina_pipe}|${dbf_pipe} ";
      }
      else
      {
          $cmd .= "|${dbf_pipe}";
      }
    }

    $dbf = 0;
    $need3d = 0;
  }
  elsif ($jurs)
  {
    if ($first)
    {
      $cmd = "${corina_first} | $jurs_pipe ";
    }
    else
    {
      if ($need3d)
      {
        $cmd .= "|${corina_pipe} | ${jurs_pipe}";
      }
      else
      {
        $cmd .= "|${jurs_pipe}";
      }
    }

    $jurs = 0;
    $need3d = 0;
  }
  elsif ($erg)
  {
    if ($first)
    {
      $cmd = "${erg_cmd_first}";
    }
    else
    {
      $cmd .= "|${erg_cmd_pipe}";
    }
    $erg = 0;
  }
  elsif ($ncerg)
  {
    if ($first)
    {
      $cmd = "${ncerg_cmd_first}";
    }
    else
    {
      $cmd .= "|${ncerg_cmd_pipe}";
    }
    $ncerg = 0;
  }
  elsif ($fperg)
  {
    if ($first)
    {
      $cmd = "${fperg_cmd_first}";
    }
    else
    {
      $cmd .= "|${fperg_cmd_pipe}";
    }
    $fperg = 0;
  }
  elsif ($cats)
  {
    if ($first)
    {
      $cmd = "${cats_cmd_first}";
    }
    else
    {
      $cmd .= "|${cats_cmd_pipe}";
    }
    $cats = 0;
  }
  elsif ($catsp)
  {
    if ($first)
    {
      $cmd = "${catsp_cmd_first}";
    }
    else
    {
      $cmd .= "|${catsp_cmd_pipe}";
    }
    $catsp = 0;
  }
  elsif ($mapO)
  {
    if ($first)
    {
      $cmd = "${mapO_cmd_first}";
    }
    else
    {
      $cmd .= "|${mapO_cmd_pipe}";
    }
    $mapO = 0;
  }
  elsif ($atr)
  {
    if ($first)
    {
      $cmd = $atr_first;
    }
    else
    {
      $cmd .= "|${atr_cmd_pipe}";
    }
    $atr = 0;
  }
  elsif ($ata)
  {
    if ($first)
    {
      $cmd = $ata_first;
    }
    else
    {
      $cmd .= "|${ata_cmd_pipe}";
    }
    $ata = 0;
  }
  elsif ($atp)
  {
    if ($first)
    {
      $cmd = $atp_first;
    }
    else
    {
      $cmd .= "|${atp_cmd_pipe}";
    }
    $atp = 0;
  }
  elsif (length($chirality_fingerprint_first) > 0)
  {
    if ($first)
    {
      $cmd = $chirality_fingerprint_first;
    }
    else
    {
      $cmd .= "|${chirality_fingerprint_pipe}";
    }
    $fingerprints_specified--;
  }
  elsif ($user_insert)
  {
    $cmd .= $user_insert_string;
    $user_insert = 0;
  }
  elsif ($d2f)
  {
    $cmd .= $d2f_string;
    $d2f = 0;
  }
  elsif ($i2f)
  {
    $cmd .= $i2f_string;
    $i2f = 0;
  }

# Models are complicated by the fact that we don't want the fingerprints of the model
# joining what we are asked to produce. So, we temporarily protect any fingerprints 
# already in the stream, then after the model has been evaluated, we remove anything unprotected
# see the sed commands

  elsif (@model)
  {
    $cmd .= "|sed -e 's/^FP/QFP/' -e 's/^NC/QNC/'";    # protect fingerprints already present
    my $m = shift(@model);

    my $buckets = 10;
    my $replicates = 1;
    my $mdir;

    if (-d $m)
    {
      $mdir = $m;
    }
    elsif ($m =~ /(\S+):(\d+):(\d+)$/)
    {
      $replicates = $2;
      $buckets = $3;
      $mdir = $1;
    }
    elsif ($m =~ /(\S+):(\d+)$/)
    {
      $replicates = $2;
      $mdir = $1;
    }
    else
    {
      print STDERR "Unrecognised model specification '${m}'\n";
      exit 2;
    }

    if (! -d $mdir)
    {
       print STDERR "Missing model directory '${mdir}'\n";
       exit(1);
    }

    if ($first)
    {
      $cmd = "${model_cmd_first} -mdir ${mdir} ";
    }
    else
    {
      $cmd .= "|${model_cmd_pipe} -mdir ${mdir} ";
    }

    my $ndx = @model;

    $cmd .= " -J TAG=NCMDL${ndx}";    # create unique tag for each model

    $cmd .= " -J RPL=${replicates}" if ($replicates > 1);
    $cmd .= " -J BKT=${buckets}" if (10 != $buckets);
    if ($first)
    {
      $cmd .= " FILE";
    }
    else
    {
      $cmd .= " -";
    }
    $cmd .= "|sed -e 's/^NCMDL/QNCMDL/' -e '/^FP/d' -e '/^NC/d' -e 's/^QFP/FP/' -e 's/^QNC/NC/'";  # protect result just produced, remove model fingerprints, unprotect previously protected
  }
  elsif ($flatten_counted)
  {
    $cmd .= "| ${flatten_counted_pipe}";
    $flatten_counted = 0;
  }
  elsif (@yoyo)
  {
    if ($first)     # will never happen
    {
      $cmd = shift(@yoyo) . " FILE";
    }
    else
    {
      $cmd .= "|" . shift(@yoyo) . " -";
    }
  }
  else
  {
    print STDERR "Huh, what fingerprints did you specify?\n";
    exit 2;
  }

  $fingerprints_specified--;

  print STDERR "cmd is '$cmd'\n" if ($verbose);

  $first = 0;
}

if (length($sp2fb))
{
  my $sp2fb_cmd = find_executable('gfp_sparse_to_fixed');

  $cmd .= "| ${sp2fb_cmd} ${sp2fb} -";
}

if (length($join_existing_file))
{
  my $tdt_join = find_executable('tdt_join');
  $cmd .= "| ${tdt_join} -d -t PCN - -c 1 -C 1 ${join_existing_file}";
}

# quick check to see that all command line arguments are valid files

my $i = $argptr;
while ($i < @ARGV)
{
  my $fname = $ARGV[$i++];
  next if ($fname eq '-');
  if (! -s $fname)
  {
    print STDERR "Missing or empty input file '${fname}'\n";
    exit(1);
  }
}

print STDERR "Command is '$cmd'\n" if ($verbose);

my $cmd_file;
if (length($cmd_fname))
{
  open($cmd_file, ">$cmd_fname");
  die "Cannot open command file '${cmd_fname}'" unless $cmd_file;
}

while ($argptr < @ARGV)
{
  my $fname = $ARGV[$argptr++];

  if ($fname eq '-')
  {
    if (length($input_qualifiers))
    {
    }
    elsif ($work_as_tdt_filter)
    {
    }
    else
    {
      $fname = "-i smi -";
    }
    die "STDIN input does not work with descriptors" if (length($descriptors));
  }
  else
  {
    die "Missing or empty input file '$fname'" unless (-s $fname);
  }

  my $cmdcopy = $cmd;

  $cmdcopy =~ s/\bFILE\b/$fname/;

  if (length($descriptors))
  {
    my $gfp_add_descriptors;
    $gfp_add_descriptors = find_executable("gfp_add_descriptors");

    $cmdcopy .= "|$gfp_add_descriptors -D DDAT - $dsctmp2";
  }

  print STDERR "Changed to '$cmdcopy'\n" if ($verbose);

  system($cmdcopy);

  print $cmd_file "${cmdcopy}\n" if (length($cmd_fname));
}

close($cmd_file) if (length($cmd_fname));

if (defined($dsctmp2))
{
  unlink($dsctmp2);
}

foreach (@files_to_be_deleted)
{
  unlink($_);
}

sub check_executable 
{
  use strict 'vars';

  my $fname = $_[0];

  die "Missing or inaccessible programme '$fname'" unless (-s $fname && -x $fname);
}

sub find_executable 
{
  use strict 'vars';

  use vars qw (@bindir);

  my $fname = $_[0];

  my $dir;

  foreach $dir (@bindir)
  {
    return "$dir/$fname" if (-x "$dir/$fname");
  }

  foreach $dir (@bindir)
  {
    return "${dir}/${fname}" if (-x "${dir}/${fname}.sh");
  }

  die "Cannot find executable for '$fname'\n";
}

sub gather_isostere_reactions 
{
  use strict 'vars';

  my $dir = $_[0];

  opendir(D, $dir) || die "Cannot open isostere directory '$dir': $!";

  my $fname;

  my $rc = "";

  while (defined($fname = readdir(D)))
  {
    next unless ($fname =~ /\.rxn$/);

    $rc .= " -R $dir/$fname";
  }

  closedir(D);

  return $rc;
}

