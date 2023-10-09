/*
  Interface to Biobyte clogp
*/

#include <iostream>
#include <limits>
#include <memory>

extern "C" void stdSMILES(char *);
extern "C" int clogp_(int *, float *, int *);
extern "C" int smilin_(char *buf, long int *warn, long int *ok, long buf_len);

extern "C" int initialise_BioByte_Arrays();

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/charge_assigner.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

using std::cerr;
using std::endl;

const char * prog_name = nullptr;

static int verbose = 0;

static int molecules_read = 0;

static Chemical_Standardisation chemical_standardisation;

static int compute_clogd = 0;

static Charge_Assigner charge_assigner;

#define CHARGE_ASSIGNER_MULTIPLIER 100

static int reduce_to_largest_fragment = 1;

static IWString tag;

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");

static int bit_replicates = 9;  // 1 is not enough, 20 reduces performance, but this is largely a guess

static extending_resizable_array<int> clogp_fp_histogram;

static int function_as_filter = 0;

static int output_like_biobyte = 0;

static int fault_tolerant = 0;

// A value that will be used if 2 == fault_tolerant.
static float value_to_print_for_failed_calculations = 0.0f;

static int failed_calculations = 0;

static int smilin_failed = 0;

static int write_clogp_error_level = 0;

static Accumulator<float> acc_clogp, acc_clogd;

static int read_smiles_as_text = 0;

static int ntest = 0;
static int test_failures = 0;

static int apply_clogp_to_charged_molecule = 0;

#define NEUTRAL_UNCHARGED 0.34
#define QUAT -0.29
#define MULTI_MINUS 3.09
#define POSITIVE3 3.38
#define POSITIVE2 2.14
#define ZWIT_1_MINUS 0.71
#define ZWIT_2_MINUS 2.10

class LogD
{
  public:
    int quats_encountered;
    int no_charge_encountered;
    int multi_minus;
    int positive3;
    int positive2;
    int zwit_1_minus;
    int zwit_2_minus;

    LogD();
};

LogD::LogD()
{
   quats_encountered = 0;
   no_charge_encountered = 0;
   multi_minus = 0;
   positive3 = 0;
   positive2 = 0;
   zwit_1_minus = 0;
   zwit_2_minus = 0;

   return;
}

static LogD logd_stats;

static char output_separator = ' ';

static int clogp_max_error_level = 60;

static int molecules_above_clogp_max_error_level = 0;

/*
  Never implemented the clogp filtering, finish if anyone ever wants it
*/

static int filter_clogp_values = 0;

static float min_clogp = -std::numeric_limits<float>::max();
static float max_clogp =  std::numeric_limits<float>::max();

static int discarded_for_clogp_too_low = 0;
static int discarded_for_clogp_too_high = 0;

static IWString_and_File_Descriptor stream_for_discarded_by_filter;

static void
usage (int rc)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << "Interface to Biobyte clogp\n";
  cerr << "  -J <tag>      produce fingerprint\n";
  cerr << "  -p <n>        number of bit replicates in fingerprints\n";
  cerr << "  -f            work as a TDT filter\n";
  cerr << "  -N ...        charge assigner specification\n";
  cerr << "  -Q            do clogp calculation on charged molecule\n";
  cerr << "  -D            also compute clogd\n";
  cerr << "  -e <level>    max valid error level (default 60 - which means ingored)\n";
  cerr << "  -q            include clogp error level as descriptor\n";
  cerr << "  -c <min>      discard molecules with clogp values less than <min>\n";
  cerr << "  -C <max>      discard molecules with clogp values more than <min>\n";
  cerr << "  -B <fname>    write discarded molecules to <fname>\n";
  cerr << "  -Y <prefix>   insert <prefix> before all descriptor names\n";
  cerr << "  -U .          default fault tolerant mode - does not quit on failed calculation\n";
  cerr << "  -U <number>   write molecules that fail, but use <number> as the clogp/clogd values\n";
//cerr << "  -i <type>     input specification\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -v            verbose output\n";
  // clang-format off

  exit(rc);
}

static int
handle_rejected_molecule (Molecule & m,
                          const float clogp)
{
  if (stream_for_discarded_by_filter.is_open())
  {
    stream_for_discarded_by_filter << m.smiles() << ' ' << m.name() << ' ' << clogp << '\n';
    stream_for_discarded_by_filter.write_if_buffer_holds_more_than(4096);
    return 1;
  }

  return 1;
}

static int
do_output_like_biobyte(Molecule & m,
                       const float clogp,
                       const int errlvl,
                       IWString_and_File_Descriptor & output)
{
  output << clogp << ' ' << errlvl << ' ' << m.smiles() << ' ' << m.name() << '\n';

  return 1;
}

static void
preprocess (Molecule & m)
{
  if (reduce_to_largest_fragment)
    m.reduce_to_largest_fragment();

  m.remove_all_chiral_centres();
  m.revert_all_directional_bonds_to_non_directional();

  if (chemical_standardisation.active())
    chemical_standardisation.process(m);

  return;
}

static int
do_descriptor_output (const Molecule & m,
                      const float clogp,
                      const int   errlvl,
                      const float clogd_offset,
                      IWString_and_File_Descriptor & output)
{
  append_first_token_of_name(m.name(), output);
  output << output_separator << clogp;
  if (write_clogp_error_level)
    output << output_separator << errlvl;
  if (compute_clogd)
    output << output_separator << (clogp - clogd_offset) << output_separator << clogd_offset;

  output << '\n';

  return 1;
}

static int
convert_computed_to_positive_int (const float f)
{
  int rc = static_cast<int>(f + 5.4999F);

  if (rc <= 0)
    return 1;
  else
    return rc;
}

static int
write_fingerprint (const Sparse_Fingerprint_Creator & sfc,
                   IWString_and_File_Descriptor & output)
{
  IWString tmp;
  sfc.daylight_ascii_form_with_counts_encoded(tag, tmp);

  output << tmp << '\n';

  return 1;
}

static int
do_fingerprint_output (Molecule & m,
                       const float clogp,
                       const float clogd,
                       IWString_and_File_Descriptor & output)
                       
{
  const int int_logp = convert_computed_to_positive_int(clogp);
  if (verbose)
    clogp_fp_histogram[int_logp]++;

  Sparse_Fingerprint_Creator sfc;

  if (compute_clogd)
  {
    const int int_logd = convert_computed_to_positive_int(clogd);

    for (int i = 0; i < bit_replicates; i++)
    {
      sfc.hit_bit(i * 2, int_logp);
      sfc.hit_bit(i * 2 + 1, int_logd);
    }
  }
  else
  {
    for (int i = 0; i < bit_replicates; i++)
    {
      sfc.hit_bit(i * 2, int_logp);
    }
  }

  if (function_as_filter)
    write_fingerprint(sfc, output);
  else
  {
    output << smiles_tag << m.smiles() << ">\n";
    output << identifier_tag << m.name() << ">\n";
    write_fingerprint(sfc, output);
    output << "|\n";
  }

  return 1;
}


static int
do_output (Molecule & m,
           const float clogp,
           const int   errlvl,
           const float clogd_offset,
           IWString_and_File_Descriptor & output)
{
  if (verbose)
  {
    acc_clogp.extra(clogp);
    if (compute_clogd)
      acc_clogd.extra(clogp - clogd_offset);
  }

  if (tag.length())
    return do_fingerprint_output(m, clogp, clogp - clogd_offset, output);
  else if (output_like_biobyte)
    return do_output_like_biobyte(m, clogp, errlvl, output);
  else
    return do_descriptor_output(m, clogp, errlvl, clogd_offset, output);
}

#ifdef DO_CHARGE_ASSIGNER
static double
do_charge_assigner (Molecule & m,
                    Charge_Assigner & charge_assigner,
                    formal_charge_t * fc)
{
  const int matoms = m.natoms();

  std::fill_n(fc, matoms, static_cast<formal_charge_t>(0));

  const int n = charge_assigner.number_elements();

  double rc = 0.0;

  Molecule_to_Match target(&m);

  for (int i = 0; i < n; i++)
  {
    Substructure_Results sresults;

    Substructure_Hit_Statistics * q = charge_assigner[i];

    int nhits = q->substructure_search(target, sresults);

    if (0 == nhits)
      continue;

    double logd_offset;
    if (! q->numeric_value(logd_offset))   // hmmm, no logd offset!!
      continue;

//  cerr << nhits << " to query '" << q->comment() << "' offset " << logd_offset << '\n';
    
    for (int j = 0; j < nhits; j++)
    {
      const Query_Atoms_Matched * qam = sresults.query_atoms_matching(j);

      const Set_of_Atoms * e = sresults.embedding(j);

      for (int k = 0; k < qam->number_elements(); k++)
      {
        atom_number_t l = e->item(k);

        if (0 != fc[l])    // already hit by this or something else
          continue;

        const Substructure_Atom * a = qam->item(k);

        double charge_specification;
        if (! a->numeric_value(charge_specification))  // no charge placed here
          continue;

        rc += logd_offset;
        if (charge_specification < 0.0)
          fc[l] = -1;
        else
          fc[l] = 1;
      }
    }
  }

  return rc;
}
#endif

static int
assign_formal_charges(Molecule & m,
                      Charge_Assigner & charge_assigner,
                      formal_charge_t * & fc)
{
  if (nullptr == fc)
  {
    const int matoms = m.natoms();
    fc = new formal_charge_t[matoms];
    std::fill_n(fc, matoms, 0);
  }

  charge_assigner.process(m, fc);

  return 1;    // always
}

static double
compute_clogd_offset(Molecule & m,
                     Charge_Assigner & charge_assigner,
                     formal_charge_t * & fc)
{
  if (nullptr == fc)
    assign_formal_charges(m, charge_assigner, fc);

  const int matoms = m.natoms();

  int npos = 0;
  int nneg = 0;

  for (int i = 0; i < matoms; i++)
  {
    const auto fci = fc[i];

    if (0 == fci)
      ;
    else if (fci > 0)
      npos++;
    else
      nneg++;
  }

  if (0 == npos && 0 == nneg)
  {
    logd_stats.no_charge_encountered++;
    return NEUTRAL_UNCHARGED;
  }

  if (nneg > 1 && 0 == npos)
  {
    logd_stats.multi_minus++;
    return MULTI_MINUS;
  }

  if (npos >= 3 && 0 == nneg)
  {
    logd_stats.positive3++;
    return POSITIVE3;
  }

  if (2 == npos && 0 == nneg)
  {
    logd_stats.positive2++;
    return POSITIVE2;
  }

  if (npos > 0 && 1 == nneg)
  {
    logd_stats.zwit_1_minus++;
    return ZWIT_1_MINUS;
  }

  if (npos > 0 && nneg > 1)
  {
    logd_stats.zwit_2_minus++;
    return ZWIT_2_MINUS;
  }

  if (npos > 0 && 0 == nneg)
    ;
  else if (0 == npos && nneg > 0)
    ;
  else
  {
    cerr << "compute_clogd_offset:unusual charge state " << npos << "+ and " << nneg << "-\n";
    return 0;
  }

  double offset = 0.0;

  for (int i = 0; i < matoms; i++)
  {
    auto fci = fc[i];

    if (0 == fci)
      continue;

    if (fci < 0)
      fci = -fci;

    const int query_number = fci / CHARGE_ASSIGNER_MULTIPLIER;

//  cerr << "Atom " << i << " matched to query number " << query_number << '\n';

    double d;
    if (! charge_assigner[query_number]->numeric_value(d))
      continue;

    offset += static_cast<double>(d);
  }

  return offset;
}

// 

static int
do_clogp_calculation (IWString & smiles,
                      float & clogp,
                      int & errlvl)
{
  clogp = 0.0f;
  errlvl = 0;

//IWString initial_smiles(smiles);

  smiles.null_terminated_chars();

  char * s = const_cast<char *>(smiles.rawchars());

  stdSMILES(s);

//smiles.resize(strlen(s));

//if (smiles != initial_smiles)
//  cerr << "stdSMILES update " << initial_smiles << " to " << smiles << '\n';

  long int warn = 0;
  long int ok = 0;
//cerr << "Long int " << sizeof(long int) << '\n';

  smilin_(s, &warn, &ok, strlen(s));

  smiles.resize_keep_storage(strlen(s));

  if (! ok)
  {
    cerr << "smilin_ failed " << s << " warn " << warn << " ok " << ok << '\n';
    smilin_failed++;
    return 0;
  }

  int outlev = 0;
  clogp_(&outlev, &clogp, &errlvl);

  return 1;
}

static int
close_enough (const float c1, const float c2)
{
  if (c1 == c2)
    return 1;

  if (abs(c1 - c2) < 1.0e-04)
    return 1;

  return 0;
}

static int
do_tests (Molecule & m,
          IWString_and_File_Descriptor & output)
{
  IWString s0 = m.smiles();

  float clogp;
  int errlvl;
  if (! do_clogp_calculation(s0, clogp, errlvl))
  {
    cerr << "Cannot compute logP first attempt " << s0 << ' ' << m.name() << '\n';
    return 0;
  }

  float x;
  IWString s = m.unique_smiles();
  if (! do_clogp_calculation(s, x, errlvl))
  {
    cerr << "Cannot compute logP unique smiles " << s << ' ' << m.name() << '\n';
    return 0;
  }

  if (! close_enough(clogp, x))
  {
    cerr << "Test mismatch. " << s0 << ' ' << clogp << ' ' << s << ' ' << x << '\n';
    test_failures++;
    return 0;
  }

  for (int i = 0; i < ntest; ++i)
  {
    s = m.random_smiles();
    Molecule q;
    if (! q.build_from_smiles(s))
    {
      cerr << "Cannot interpret random smiles variant " << s << ' ' << m.name() << '\n';
      test_failures++;
      return 0;
    }

    IWString s2 = q.smiles();

    int myerrlvl;
    if (! do_clogp_calculation(s2, x, myerrlvl))
    {
      cerr << "clogP failed for variant " << s2 << ' ' << m.name() << '\n';
      return 0;
    }

    if (! close_enough(clogp, x))
    {
      cerr << "Clogp mismatch " << s0 << ' ' << clogp << ' ' << s2 << ' ' << x << '\n';
      test_failures++;
      return 0;
    }

    if (myerrlvl != errlvl)
    {
      cerr << "Error Level mismatch " << s0 << ' ' << errlvl << ' ' << s2 << ' ' << myerrlvl << '\n';
    }
  }

  return 1;
}

static int
print_fault_tolerant_output (Molecule & m,
                             IWString_and_File_Descriptor & output)
{
  return do_output(m, value_to_print_for_failed_calculations, 99, value_to_print_for_failed_calculations, output);
}

static int
bb_clogp (Molecule & m,
          IWString_and_File_Descriptor & output)
{
  preprocess(m);

  const int matoms = m.natoms();

  if (matoms > 256)
  {
    cerr << "Molecule too large for Biobyte fixed size arrays " << m.name() << '\n';
    return 1;
  }

  if (ntest > 0)
    return do_tests(m, output);

  const int has_formal_charges = m.has_formal_charges();    // before we possibly apply charges

  formal_charge_t * fc = nullptr;
  int formal_charges_assigned = 0;

  if (apply_clogp_to_charged_molecule)
  {
    if (! assign_formal_charges(m, charge_assigner, fc))
    {
      cerr << "iwclogp:cannot assign formal charges " << m.smiles() << ' ' << m.name() << '\n';
      return 0;
    }
    for (int i = 0; i < matoms; ++i) 
    {
      if (0 != fc[i])
      {
        m.set_formal_charge(i, fc[i] % CHARGE_ASSIGNER_MULTIPLIER);
        formal_charges_assigned++;
      }
    }
  }

//IWString s = m.unique_smiles();
  IWString s = m.smiles();

  float clogp;
  int errlvl;
  if (! do_clogp_calculation(s, clogp, errlvl))
  {
    cerr << "clogP failed for " << s << ' ' << m.name() << '\n';
    failed_calculations++;
    if (2 == fault_tolerant)
      return print_fault_tolerant_output(m, output);
    return fault_tolerant;
  }

  if (errlvl >= clogp_max_error_level)
  {
    if (verbose > 1)
      cerr << "CLOGP failed " << m.smiles() << ' ' << m.name() << " errlvl " << errlvl << '\n';

    molecules_above_clogp_max_error_level++;

    if (2 == fault_tolerant)
      return print_fault_tolerant_output(m, output);

    return handle_rejected_molecule(m, 99.0);   // kind of confusing since they will be mixed in with things that are filtered...
  }

  if (! filter_clogp_values)
    ;
  else if (clogp < min_clogp || clogp > max_clogp)
  {
    if (clogp < min_clogp)
      discarded_for_clogp_too_low++;
    else
      discarded_for_clogp_too_high++;

    return handle_rejected_molecule(m, clogp);
  }

  if (! compute_clogd)
    return do_output(m, clogp, 0, 0.0f, output);

  if (has_formal_charges)    // was determined before we possibly applied charges
  {
    logd_stats.quats_encountered++;
    return do_output(m, clogp, errlvl, QUAT, output);
  }

  float clogd_offset;

  if (! compute_clogd)
    clogd_offset = 0.0f;
  else if (apply_clogp_to_charged_molecule && 0 == formal_charges_assigned)    // already computed
    clogd_offset = NEUTRAL_UNCHARGED;
  else
    clogd_offset = compute_clogd_offset(m, charge_assigner, fc);

  if (nullptr != fc)
    delete [] fc;

  return do_output(m, clogp, errlvl, clogd_offset, output);
}

static int
bb_clogp (data_source_and_type<Molecule> & input,
          IWString_and_File_Descriptor & output)
{
  Molecule * m;
  while (NULL != (m = input.next_molecule()))
  {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    if (! bb_clogp(*m, output))
      return 0;

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
bb_clogp_filter_record (const_IWSubstring buffer,    // note local copy
                        IWString_and_File_Descriptor & output)
{
  assert (buffer.starts_with(smiles_tag));
  assert (buffer.ends_with('>'));

  buffer.remove_leading_chars(smiles_tag.length());
  buffer.chop();

  Molecule m;

  if (! m.build_from_smiles(buffer))
  {
    cerr << "Cannot interpret smiles '" << buffer << "'\n";
    return 0;
  }

  return bb_clogp(m, output);
}

static int
bb_clogp_filter (iwstring_data_source & input,
                 IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    output << buffer << '\n';
    output.write_if_buffer_holds_more_than(4096);
//  cerr << "Processed '" << buffer << "'\n";

    if (! buffer.starts_with(smiles_tag))
      continue;

    if (! bb_clogp_filter_record(buffer, output))
    {
      cerr << "Fatal error processing '" << buffer << "', read " << input.lines_read() << " lines\n";
      return 0;
    }
  }

  return 1;
}

static int
bb_clogp_filter (const char * fname,
                 IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);
  if (! input.good())
  {
    cerr << "bb_clogp_filter:cannot open '" << fname << "'\n";
    return 0;
  }

  return bb_clogp_filter(input, output);
}

static int
bb_clogp (const char * fname, FileType input_type, 
          IWString_and_File_Descriptor & output)
{
  assert(NULL != fname);

  if (input_type == FILE_TYPE_INVALID)
  {
    input_type = discern_file_type_from_name(fname);
    assert(0 != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good())
  {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1)
    input.set_verbose(1);

  return bb_clogp(input, output);
}


static int
text_processing_record (IWString & smiles,
                        const IWString & id,
                        IWString_and_File_Descriptor & output)
{
  float clogp;
  int errlvl;

  if (! do_clogp_calculation(smiles, clogp, errlvl))
  {
    cerr << "Fatal error processing " << smiles << " '" << id << "'\n";
    return fault_tolerant;
  }

  if (output_like_biobyte)
  {
    output << clogp << ' ' << errlvl << ' ' << smiles << ' ' << id << '\n';
    return 1;
  }

  output << id << output_separator << clogp;
  if (write_clogp_error_level)
    output << output_separator << errlvl;
  output << '\n';

  return 1;
}

static int
text_processing_record (const const_IWSubstring & buffer,
                        IWString_and_File_Descriptor & output)
{
  IWString smiles, id;

  if (buffer.split(smiles, ' ', id))
    ;
  else if (1 == buffer.nwords())
    smiles = buffer;
  else
  {
    cerr << "Invalid input\n";
    return 0;
  }

  return text_processing_record(smiles, id, output);
}

static int
text_processing (iwstring_data_source & input,
                 IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    if (! text_processing_record(buffer, output))
    {
      cerr << "Fatal error processing '" << buffer << "', read " << input.lines_read() << " lines\n";
      return 0;
    }

    output.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

static int
text_processing (const char * fname,
                 IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return text_processing(input, output);
}

/*
  input will be $SMI<....>
*/

static int
filter_with_text_processing_record (const_IWSubstring buffer,   // note our own copy
                            IWString_and_File_Descriptor & output)
{
  buffer.remove_leading_chars(smiles_tag.length());
  buffer.chop();

  IWString smiles(buffer);

  float clogp;
  int errlvl;
  if (! do_clogp_calculation(smiles, clogp, errlvl))
  {
    cerr << "clogP failed " << smiles << '\n';
    return fault_tolerant;
  }

  return 1;
}

static int
filter_with_text_processing(iwstring_data_source & input,
                            IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    output << buffer << '\n';

    output.write_if_buffer_holds_more_than(4096);

    if (! output.starts_with(smiles_tag))
      continue;

    if (! filter_with_text_processing_record(buffer, output))
    {
      cerr << "Fatal error processing '" << buffer << "', read " << input.lines_read() << " lines\n";
      return 0;
    }
  }

  return 1;
}

static int
filter_with_text_processing(const char * fname,
                            IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return filter_with_text_processing(input, output);
}

static int
parse_range (const Command_Line & cl,
             char lflag,
             char uflag,
             float & zmin,
             float & zmax,
             const char * clogp_or_clogd)
{
  if (cl.option_present(lflag))
  {
    if (! cl.value(lflag, zmin))
    {
      cerr << "The min " << clogp_or_clogd << " value(-" << lflag << ") must be a valid floating point number\n";
      return 0;
    }

    if (verbose)
      cerr << "Will discard molecules with " << clogp_or_clogd << " values < " << min_clogp << '\n';
  }

  if (cl.option_present(uflag))
  {
    if (! cl.value(uflag, zmax))
    {
      cerr << "The max " << clogp_or_clogd << " value(-C) must be a valid floating point number\n";
      usage(3);
    }

    if (cl.option_present(lflag) && zmax < zmin)
    {
      cerr << "Inconsistent min " << zmin << " and max " << zmax << " " << clogp_or_clogd << " values\n";
      return 3;
    }

    if (verbose)
      cerr << "Will discard molecules with " << clogp_or_clogd << " values > " << zmax <<'\n';
  }

  return 1;
}

static int
bb_clogp (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vA:E:i:g:lJ:p:fB:Y:e:qN:st:c:C:k:U:IQD");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('A'))
  {
    if (! process_standard_aromaticity_options(cl, verbose, 'A'))
    {
      cerr << "Cannot initialise aromaticity specifications\n";
      usage(5);
    }
  }
  else 
    set_global_aromaticity_type(Daylight);

  if (cl.option_present('E'))
  {
    if (! process_elements(cl, verbose, 'E'))
    {
      cerr << "Cannot initialise elements\n";
      return 6;
    }
  }

  if (cl.option_present('g'))
  {
    if (! chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g'))
    {
      cerr << "Cannot process chemical standardisation options (-g)\n";
      usage(32);
    }
  }

  if (cl.option_present('f'))
  {
    function_as_filter = 1;

    if (verbose)
      cerr << "Will work as a TDT filter\n";
  }

  if (cl.option_present('U'))
  {
    const_IWSubstring u = cl.string_value('U');

    if ('.' == u)
      fault_tolerant = 1;
    else if (u.numeric_value(value_to_print_for_failed_calculations))
      fault_tolerant = 2;
    else
    {
      cerr << "Unrecognised fault tolerance directive '" << u << "'\n";
      usage(1);
    }

    if (verbose)
      cerr << "Will run in fault tolerant mode\n";
  }

  if (cl.option_present('I'))
  {
    output_like_biobyte = 1;
    if (verbose)
      cerr << "Will write output like BioByte\n";
  }

  if (cl.option_present('e'))
  {
    if (! cl.value('e', clogp_max_error_level) || clogp_max_error_level < 0)
    {
      cerr << "The maximum valid error value(-e option) must be a whole positive number\n";
      usage(13);
    }

    if (verbose)
      cerr << "Error levels above " << clogp_max_error_level << " will be considered invalid\n";
  }

  if (cl.option_present('J'))
  {
    cl.value('J', tag);

    if (! tag.ends_with('<'))
      tag << '<';

    if (verbose)
      cerr << "logP and logD written as fingerprints\n";

    if (cl.option_present('p'))
    {
      if (! cl.value('p', bit_replicates) || bit_replicates <= 0)
      {
        cerr << "The bit replicates value(-p) must be a whole +ve number\n";
        usage(2);
      }

      if (verbose)
        cerr << "Fingerprint bit replicates set to " << bit_replicates << '\n';
    }
  }

  if (cl.option_present('Q'))
  {
    apply_clogp_to_charged_molecule = 1;

    if (verbose)
      cerr << "Will do clogp calculation on charged molecules\n";
  }

  if (cl.option_present('D'))
  {
    compute_clogd = 1;

    if (verbose)
      cerr << "Will compute clogd\n";
  }

  if (apply_clogp_to_charged_molecule || compute_clogd)
  {
    if (! cl.option_present('N'))
    {
      cerr << "Request charged molecule and/or clogd calculation, but no charge assigner specified\n";
      usage(1);
    }

    if (! charge_assigner.construct_from_command_line(cl, verbose, 'N'))
    {
      cerr << "Cannot initialise charge assigner\n";
      return 1;
    }

    charge_assigner.set_apply_charges_to_molecule(0);
    charge_assigner.set_assigned_charge_multiplier(CHARGE_ASSIGNER_MULTIPLIER);
  }
  else if (cl.option_present('N'))
    cerr << "Charge assigner specified (-N) but no need, ignored\n";

  if (cl.option_present('c') || cl.option_present('C'))
  {
    if (function_as_filter)
    {
      cerr << "Sorry, the -c and -C options cannot be used with the -f option\n";
      return 1;
    }

    if (! parse_range(cl, 'c', 'C', min_clogp, max_clogp, "clogp"))
    {
      cerr << "Cannot parse min and max clogp specification(s)\n";
      usage(3);
    }

    filter_clogp_values = 1;
  }

  if (cl.option_present('k'))
  {
    int k;
    if (! cl.value('k', k) || k < 1)
    {
      cerr << "The default float precision option (-k) must be a whole +ve integer\n";
      usage(1);
    }

    set_default_iwstring_float_concatenation_precision(k);

    if (verbose)
      cerr << "Default float string concatenation " << k << '\n';
  }
  else
    set_default_iwstring_float_concatenation_precision(4);

  if (cl.option_present('s'))
  {
    if (compute_clogd)
    {
      cerr << "Cannot use the -s option (text input) when clogD calculation requested\n";
      return 1;
    }

    read_smiles_as_text = 1;

    if (verbose)
      cerr << "Input files not interpreted as molecules, just text\n";
  }

  if (cl.option_present('t'))
  {
    if (! cl.value('t', ntest) || ntest < 1)
    {
      cerr << "The number of tests done (-t) must be a whole +ve number\n";
      usage(1);
    }

    if (verbose)
      cerr << "Will perform " << ntest << " tests on each input molecule\n";
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i'))
  {
    if (! process_input_type(cl, input_type))
    {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  }
  else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-"))
    input_type = FILE_TYPE_SMI;
  else if (! all_files_recognised_by_suffix(cl))
    return 4;

  if (cl.empty())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.option_present('B'))
  {
    IWString b = cl.option_value('B');

    if (! b.ends_with(".smi"))
      b << ".smi";

    if (! stream_for_discarded_by_filter.open(b.null_terminated_chars()))
    {
      cerr << "Cannot open stream for discarded molecules '" << b << "'\n";
      return 3;
    }

    if (verbose)
      cerr << "Discarded molecules written to '" << b << "'\n";
  }

  initialise_BioByte_Arrays();

  IWString_and_File_Descriptor output(1);

  if (output_like_biobyte)
//  output << "clogP Error Smiles Id\n";
    ;
  else if (0 == tag.length())
  {
    IWString prefix;
    if (cl.option_present('Y'))
    {
      cl.value('Y', prefix);
      if (verbose)
        cerr << "Will prepend '" << prefix << "' to descriptor names\n";
    }

    output << "Name" << output_separator << prefix << "clogp";
    if (cl.option_present('q'))
    {
      output << output_separator << prefix << "clogp_errlvl";
      write_clogp_error_level = 1;
      if (verbose)
        cerr << "Will write clogp error estimate\n";
    }

    if (compute_clogd)
      output << output_separator << prefix << "brnsclogD" << output_separator << prefix << "brnsclogDoffset";

    output << "\n";
  }

  int rc = 0;
  if (read_smiles_as_text && function_as_filter)
  {
    for (int i = 0; i < cl.number_elements(); ++i)
    {
      if (! filter_with_text_processing(cl[i], output))
      {
        rc = i + 1;
        break;
      }
    }
  }
  else if (read_smiles_as_text)
  {
    for (int i = 0; i < cl.number_elements(); ++i)
    {
      if (! text_processing(cl[i], output))
      {
        rc = i + 1;
        break;
      }
    }
  }
  else if (function_as_filter)
  {
    if (! bb_clogp_filter(cl[0], output))
      rc = 1;
  }
  else
  {
    for (int i = 0; i < cl.number_elements(); i++)
    {
      if (! bb_clogp(cl[i], input_type, output))
      {
        rc = i + 1;
        break;
      }
    }
  }

  output.flush();

  if (verbose)
  {
    cerr << "Read " << molecules_read << " molecules\n";
    if (smilin_failed)
      cerr << smilin_failed << " molecules failed BioByte smiles interpretation\n";
    if (failed_calculations)
      cerr << failed_calculations << " molecules failed clogp\n";
    if (molecules_above_clogp_max_error_level)
      cerr << molecules_above_clogp_max_error_level << " above clogp max\n";
    if (acc_clogp.n())
    {
      cerr << "clogP values btw " << acc_clogp.minval() << " and " << acc_clogp.maxval() << " ave " << static_cast<float>(acc_clogp.average()) << '\n';
      if (compute_clogd)
        cerr << "clogD values btw " << acc_clogd.minval() << " and " << acc_clogd.maxval() << " ave " << static_cast<float>(acc_clogd.average()) << '\n';
    }

    for (int i = 0; i < clogp_fp_histogram.number_elements(); ++i)
    {
      if (clogp_fp_histogram[i] > 0)
        cerr << clogp_fp_histogram[i] << " molecules in range " << i << '\n';
    }
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = bb_clogp(argc, argv);

  return rc;
}
