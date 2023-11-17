/*
  Reads the output of Marvin computations and produces a fingerprint
*/

#include <stdlib.h>
#include <iostream>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"

#include "Foundational/iwmisc/sparse_fp_creator.h"
#include "Foundational/accumulator/accumulator.h"
#include "Foundational/iwmisc/set_or_unset.h"

using std::cerr;
using std::endl;

const char * prog_name = nullptr;

static int verbose = 0;

static int reading_pKx_data = 0;

static int process_clogp = 0;
static int process_clogD = 0;
static int process_apKa1 = 0;
static int process_bpKa1 = 0;

static int clogp_column = -1;
static int clogD_column = -1;
static int apKa1_column = -1;
static int bpKa1_column = -1;

static int clogp_offset = -1;
static int clogd_offset = -1;
static int apKa1_offset = -1;
static int bpKa1_offset = -1;

static Accumulator<float> acc_clogp, acc_clogD, acc_apKa1, acc_bpKa1;

static int reading_clogD_data = 0;

static int records_read = 0;

static IWString smiles_tag("$SMI<");

static IWString identifier_tag("PCN<");

static IWString tag("NCMPKLP<");

static int fault_tolerant_mode = 0;

static int bit_replicates = 9;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Produces fingerprints from Marvin logp and pK. descriptors\n";
  cerr << " -F <tag>       tag for fingerprints\n";
  cerr << " -p <n>         number of bit replicates in fingerprints\n";
  cerr << " -O <type>      what type of input, 'pKx', 'logp', 'apK', 'bpK'\n";
  cerr << " -f             run in a fault tolerant mode\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

static int
convert_clogx_to_positive_int (float logx)
{
  int rc = static_cast<int>(logx + 5.4999F);

  if (rc <= 0)
    return 1;
  else
    return rc;
}

static int
convert_pkx_to_positive_int (float pkx)
{
  if (pkx <= -0.2f)
    return 1;

  return static_cast<int>(pkx + 2.49999f);
}

static int
assign_if_not_empty (const const_IWSubstring & s,
                     Set_or_Unset<float> & v)
{
  if (0 == s.length())
    return 1;

  float f;

  if (! s.numeric_value(f))
  {
    cerr << "Invalid numeric '" << s << "'\n";
    return 0;
  }

  v = f;

  return 1;
}

static int
marvin2gfp_pKx (const const_IWSubstring & buffer,
                IWString_and_File_Descriptor & output)
{
  const_IWSubstring token;
  int i = 0;

  const_IWSubstring smiles;
  const_IWSubstring id;

  if (! buffer.nextword(smiles, i, '\t') || ! buffer.nextword(id, i, '\t'))
  {
    cerr << "Cannot extract smiles and/or id\n";
    return 0;
  }

  const_IWSubstring clogp_s, apKa1_s, bpKa1_s;

  for (int col = 2; buffer.nextword_single_delimiter(token, i, '\t'); col++)
  {
    if (col == clogp_column)
      clogp_s = token;
    else if (col == apKa1_column)
      apKa1_s = token;
    else if (col == bpKa1_column)
      bpKa1_s = token;

//  cerr << "Checking column " << col << endl;
  }

  Set_or_Unset<float> clogp, apKa1, bpKa1;

  float f;

  if (! assign_if_not_empty(clogp_s, clogp))
    return 0;
  if (! assign_if_not_empty(apKa1_s, apKa1))
    return 0;
  if (! assign_if_not_empty(bpKa1_s, bpKa1))
    return 0;

  if (verbose)
  {
    if (clogp.value(f))
      acc_clogp.extra(f);
    if (apKa1.value(f))
      acc_apKa1.extra(f);
    if (bpKa1.value(f))
      acc_bpKa1.extra(f);
  }

  Sparse_Fingerprint_Creator sfc;

  if (clogp_offset >= 0 && clogp.value(f))
  {
    int clogp_int  = convert_clogx_to_positive_int(f);
//  cerr << "Convert " << f << " to " << clogp_int << endl;
    for (int i = 0; i < bit_replicates; i++)
    {
      sfc.hit_bit(clogp_offset + i, clogp_int);
    }
  }

  if (apKa1_offset >= 0 && apKa1.value(f))
  {
    int apKa1_int = convert_pkx_to_positive_int(f);
    for (int i = 0; i < bit_replicates; i++)
    {
      sfc.hit_bit(apKa1_offset + i, apKa1_int);
    }
  }

  if (bpKa1_offset >= 0 && bpKa1.value(f))
  {
    int bpKa1_int = convert_pkx_to_positive_int(f);
    for (int i = 0; i < bit_replicates; i++)
    {
      sfc.hit_bit(bpKa1_offset + i, bpKa1_int);
    }
  }

  output << smiles_tag << smiles << ">\n";
  output << identifier_tag << id << ">\n";

  IWString tmp;
  sfc.daylight_ascii_form_with_counts_encoded(tag, tmp);
  output << tmp << '\n';

  output << "|\n";

  return 1;
}

static int
marvin2gfp_clogD (const const_IWSubstring & buffer,
                 IWString_and_File_Descriptor & output)
{
  const_IWSubstring token;

  const_IWSubstring smiles, id, string_clogd;

  int i = 0;
  if (! buffer.nextword(smiles, i, '\t') ||
      ! buffer.nextword(id, i, '\t') ||
      ! buffer.nextword(string_clogd, i, '\t'))
  {
    cerr << "Cannot split clogd data '" << buffer << "'\n";
    return 0;
  }

  float d;

  if (! string_clogd.numeric_value(d))
  {
    cerr << "Invalid logd value '" << buffer << "'\n";
    return 0;
  }

  acc_clogD.extra(d);

  output << smiles_tag << smiles << ">\n";
  output << identifier_tag << id << ">\n";

  Sparse_Fingerprint_Creator sfc;

  int clogd_int = convert_clogx_to_positive_int(d);
  for (int i = 0; i < bit_replicates; i++)
  {
    sfc.hit_bit(clogd_offset + i, clogd_int);
  }

  IWString tmp;
  sfc.daylight_ascii_form_with_counts_encoded(tag, tmp);
  output << tmp << '\n';

  output << "|\n";

  return 1;
}

/*
  Typical pK? header record looks like

smiles:n	pH=0.00	pH=1.00	pH=2.00	pH=3.00	pH=4.00	pH=5.00	pH=6.00	pH=7.00	pH=8.00	pH=9.00	pH=10.00	pH=11.00	pH=12.00	pH=13.00	pH=14.00	logP	apKa1	bpKa1	atoms

  We need to work out the columns for 'logP    apKa1   bpKa1'
*/

static int
analyse_header_pKx (const const_IWSubstring & header)
{
  int i = 0;
  const_IWSubstring token;

  for (int col = 0; header.nextword(token, i, '\t'); col++)
  {
    if ("logP" == token)
      clogp_column = col;
    else if ("apKa1" == token)
      apKa1_column = col;
    else if ("bpKa1" == token)
      bpKa1_column = col;
  }

  if (clogp_column < 0 || apKa1_column < 0 || bpKa1_column < 0)
  {
    cerr << "Cannot identify needed header records\n";
    cerr << "'" << header << "'\n";
    return 0;
  }

  if (verbose > 1)
    cerr << "clogp in " << clogp_column << " apKa1 in " << apKa1_column << " bpKa1 in " << bpKa1_column << endl;

  return 1;
}

static int
analyse_header_clogD (const const_IWSubstring & header)
{
  return 1;
}

static int
marvin2gfp_pKx (const const_IWSubstring & header,
                iwstring_data_source & input,
                IWString_and_File_Descriptor & output)
{
  if (! analyse_header_pKx(header))
    return 0;

  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    records_read++;

    if (! marvin2gfp_pKx (buffer, output))
    {
      cerr << "Cannot process '" << buffer << "'\n";
      return 0;
    }

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
marvin2gfp_clogD (const const_IWSubstring & header,
                  iwstring_data_source & input,
                  IWString_and_File_Descriptor & output)
{
  if (! analyse_header_clogD(header))
    return 0;

  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    records_read++;

    if (! marvin2gfp_clogD (buffer, output))
    {
      cerr << "Cannot process '" << buffer << "'\n";
      return 0;
    }

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
marvin2gfp (iwstring_data_source & input,
            IWString_and_File_Descriptor & output)
{
//input.set_translate_tabs(1);

  const_IWSubstring header;

  if (! input.next_record(header))
  {
    cerr << "Cannot read header\n";
    return 0;
  }

  if (reading_clogD_data)
    return marvin2gfp_clogD(header, input, output);

  if (reading_pKx_data)
    return marvin2gfp_pKx(header, input, output);

  cerr << "What kind of input is this?\n";
  return 0;
}

static int
marvin2gfp (const char * fname,
            IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return marvin2gfp(input, output);
}


static int
marvin2gfp (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vF:O:fp:");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('F'))
  {
    cl.value('F', tag);

    if (! tag.ends_with('<'))
      tag << '<';

    if (verbose)
      cerr << "Fingerprint tag '" << tag << "'\n";
  }

// The -p option must be before the -O option

  if (cl.option_present('p'))
  {
    if (! cl.value('p', bit_replicates) || bit_replicates <= 0)
    {
      cerr << "The bit replicates value (-p) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose)
      cerr << "Fingerprint bit replicates set to " << bit_replicates << endl;
  }

  if (! cl.option_present('O'))
  {
    cerr << "Must specify what fingerprints to produce via the -O option\n";
    usage(2);
  }

  if (cl.option_present('O'))
  {
    int i = 0;
    const_IWSubstring o;
    while (cl.value('O', o, i++))
    {
      if ("pKx" == o)
      {
        process_apKa1 = 1;
        process_bpKa1 = 1;
        reading_pKx_data = 1;
      }
      else if ("clogp" == o || "logp" == o)
      {
        process_clogp = 1;
        reading_pKx_data = 1;
      }
      else if ("apK" == o)
      {
        process_apKa1 = 1;
        reading_pKx_data = 1;
      }
      else if ("bpK" == o)
      {
        process_bpKa1 = 1;
        reading_pKx_data = 1;
      }
      else if ("clogd" == o || "logd" == o)
      {
        process_clogD = 1;
        reading_clogD_data = 1;
      }
      else
      {
        cerr << "Unrecognised -O qualifier '" << o << "'\n";
        usage(1);
      }
    }
  }

  int next_offset = 0;

  if (process_clogp)
  {
    clogp_offset = next_offset;
    next_offset += bit_replicates;
  }

  if (process_apKa1)
  {
    apKa1_offset = next_offset;
    next_offset += bit_replicates;
  }

  if (process_bpKa1)
  {
    bpKa1_offset = next_offset;
    next_offset += bit_replicates;
  }

  if (process_clogD)
  {
    clogd_offset = next_offset;
    next_offset += bit_replicates;
  }

  if (0 == next_offset)
  {
    cerr << "Huh, nothing being processed, specify output with -O option\n";
    usage(2);
  }

  if (cl.option_present('F'))
  {
    cl.value('F', tag);

    if (! tag.ends_with('<'))
      tag << '<';

    if (verbose)
      cerr << "Fingerprint tag '" << tag << "'\n";
  }
  else
  {
    tag = "NCMV";
    if (process_clogp)
      tag << 'P';
    if (process_clogD)
      tag << 'D';
    if (process_apKa1)
      tag << 'A';
    if (process_bpKa1)
      tag << 'B';

    tag << '<';
  }

  if (reading_clogD_data && reading_pKx_data)
  {
    cerr << "Cannot process both Marvin logD and PKx data\n";
    return 3;
  }

  if (! reading_clogD_data && ! reading_pKx_data)
  {
    process_clogp = 1;
    process_apKa1 = 1;
    process_bpKa1 = 1;
    reading_pKx_data = 1;
  }

  if (cl.option_present('f'))
  {
    fault_tolerant_mode = 1;

    if (verbose)
      cerr << "Will skip over some kinds of errors\n";
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! marvin2gfp(cl[i], output))
    {
      rc = i + 1;
      break;
    }
  }

  output.flush();

  if (verbose)
  {
    cerr << "Processed " << records_read << " molecules\n";
    if (acc_clogp.n() > 1)
      cerr << "clogp values between " << acc_clogp.minval() << " and " << acc_clogp.maxval() << " ave " << static_cast<float>(acc_clogp.average()) << endl;
    if (acc_apKa1.n() > 1)
      cerr << "apK values between " << acc_apKa1.minval() << " and " << acc_apKa1.maxval() << " ave " << static_cast<float>(acc_apKa1.average()) << endl;
    if (acc_bpKa1.n() > 1)
      cerr << "bpK values between " << acc_bpKa1.minval() << " and " << acc_bpKa1.maxval() << " ave " << static_cast<float>(acc_bpKa1.average()) << endl;
    if (acc_clogD.n() > 1)
      cerr << "clogD values between " << acc_clogD.minval() << " and " << acc_clogD.maxval() << " ave " << static_cast<float>(acc_clogD.average()) << endl;
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = marvin2gfp(argc, argv);

  return rc;
}
