/*
  Convert output from clogp to descriptor form
  Tips for mapping file descriptors
  http://unixwiz.net/techtips/remap-pipe-fds.html
*/

#include <ctype.h>
#include <stdlib.h>
#include <unistd.h>

#include <algorithm>
#include <limits>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/charge_assigner.h"
#include "Molecule_Lib/element.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

#include <sys/types.h>
#include <sys/wait.h>

using std::cerr;
using std::endl;

static int verbose = 0;

static int molecules_read = 0;

/*
  If we have a "missing" clogp value, the missing_values string
  will contain one or three missing value symbols
*/

static IWString missing_values('.');

static Charge_Assigner charge_assigner;

#define CHARGE_ASSIGNER_MULTIPLIER 100

static IWString descriptor_prefix;

static int rerun_charge_assigner_queries = 0;

static int clogp_max_error_level = 99;

static int molecules_above_clogp_max_error_level = 0;

static Chemical_Standardisation chemical_standardisation;

static int fault_tolerant = 0;

/*
  Never implemented the clogp filtering, finish if anyone ever wants it
*/

static int filter_clogp_values = 0;

static float min_clogp = -std::numeric_limits<float>::max();
static float max_clogp = std::numeric_limits<float>::max();

static int discarded_for_clogp_too_low = 0;
static int discarded_for_clogp_too_high = 0;

static int filter_clogd_values = 0;

static float min_clogd = -std::numeric_limits<float>::max();
static float max_clogd = std::numeric_limits<float>::max();

static int discarded_for_clogd_too_low = 0;
static int discarded_for_clogd_too_high = 0;

static IWString_and_File_Descriptor stream_for_discarded_by_filter;

static int append_clogp_to_filter_output = 0;
static int append_clogp_to_rejected_output = 0;

static int truncate_to_first_token_of_name = 1;

// Write output that looks like
//   smiles id result
static int write_smiles_file = 0;

/*
  Bug. My initial implementation was supposed to have clogp+clogd, but I
  made a mistake and put out clogp twice.
  Therefore we can optionally control what gets output
*/

#define CLOGPCLOGP 1
#define CLOGPCLOGD 2
#define CLOGDCLOGD 3

static int to_fingerprint = CLOGPCLOGP;

/*
  We can produce a fingerprint if needed
*/

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");

static IWString tag;

// Adding just one bit to a large fingerprint does not do much, so we set multiple bits
// 1 is not enough, 20 reduces performance, but this is largely a guess.
// But there are plenty of models that do benefit from more bit replicates, this
// needs to be explored on a case by case basis.

static int bit_replicates = 9; 

static extending_resizable_array<int> clogp_fp_histogram;

static int work_as_tdt_filter = 0;

static int flush_after_every_molecule = 0;

class LogD {
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

LogD::LogD() {
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

static void
usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << "Converts output from Biobyte clogp to descriptors, incl Bruns logD\n";
  cerr << "  -M <string>   missing value string (default '.')\n";
  cerr << "  -e <level>    max valid error level (default 0)\n";
  cerr << "  -c <min>      work as filter, discard if clogp .lt. <min>\n";
  cerr << "  -C <max>      work as filter, discard if clogp .gt. <min>\n";
  cerr << "  -d <min>      work as filter, discard if clogd .lt. <min>\n";
  cerr << "  -D <max>      work as filter, discard if clogd .gt. <min>\n";
  cerr << "  -a            append clogp/clogd to filtered output\n";
  cerr << "  -B <fname>    write filtered molecules to <fname>\n";
  cerr << "  -b            append clogp/clogd to rejected output\n";
  cerr << "  -Y <prefix>   insert <prefix> before all descriptor names\n";
  cerr << "  -N ...        charge assigner specifications - logD\n";
  cerr << "  -r            re-run charge assigner queries individually\n";
  cerr << "  -u            fault tolerant mode\n";
  cerr << "  -U            do NOT truncate to the first token of the name\n";
  cerr << "  -k <string>   control what gets fingerprinted (pp,pd,dd)\n";
  cerr << "  -J <tag>      tag for creating fingerprints\n";
  cerr << "  -p <n>        number of bit replicates in fingerprints\n";
  cerr << "  -q            a TDT has been fed through clogp\n";
  cerr << "  -X ...        more options\n";
  cerr << "  -E autocreate autocreate elemnents\n";
  display_standard_aromaticity_options(cerr);
  display_standard_chemical_standardisation_options(cerr, 'g');
  cerr << "  -v            verbose output\n";
  // clang-format on

  exit(rc);
}

static void
preprocess(Molecule& m) {
  m.reduce_to_largest_fragment();

  if (chemical_standardisation.active()) {
    chemical_standardisation.process(m);
  }

  return;
}

static int
write_first_token_of_identifier(const_IWSubstring id,  // local copy
                                IWString_and_File_Descriptor& output) {
  id.truncate_at_first(' ');

  output << id << ' ';

  return 1;
}

static int
rejected_by_logd(Molecule& m, float clogp, float clogd) {
  if (!stream_for_discarded_by_filter.is_open()) {
    return 1;
  }

  stream_for_discarded_by_filter << m.smiles() << ' ' << m.name();

  if (append_clogp_to_rejected_output) {
    stream_for_discarded_by_filter << ' ' << clogp << ' ' << clogd;
  }

  stream_for_discarded_by_filter << '\n';

  stream_for_discarded_by_filter.write_if_buffer_holds_more_than(32768);

  return 1;
}

static int
common_write_fingerprint(const const_IWSubstring& smiles, const const_IWSubstring& id,
                         const Sparse_Fingerprint_Creator& sfc,
                         IWString_and_File_Descriptor& output) {
  output << smiles_tag << smiles << ">\n";

  if (id.length() > 0) {
    output << identifier_tag << id << ">\n";
  }

  IWString tmp;
  sfc.daylight_ascii_form_with_counts_encoded(tag, tmp);

  output << tmp << '\n';

  if (!work_as_tdt_filter) {
    output << "|\n";
  }

  if (flush_after_every_molecule) {
    output.flush();
  }

  return 1;
}

static int
convert_computed_to_positive_int(float f) {
  int rc = static_cast<int>(f + 5.4999F);

  if (rc <= 0) {
    return 1;
  } else {
    return rc;
  }
}

/*
  To convert logP and logD to fingerprints, we get lazy and
  use a non-colliding fingerprint. We should use a counted
  fixed width type, but those are less well supported in
  the svmfp framework
*/

static int
write_fingerprint_data(Molecule& m, float clogp, float offset,
                       IWString_and_File_Descriptor& output) {
  // if (m.name().length() > 0)
  //   output << identifier_tag << m.name() << ">\n";

  int int_logp = convert_computed_to_positive_int(clogp);
  int int_logd = convert_computed_to_positive_int(clogp - offset);

  // cerr << " Convert " << clogp << " to " << int_logp << '\n';

  if (verbose) {
    clogp_fp_histogram[int_logp]++;
  }

  Sparse_Fingerprint_Creator sfc;

  if (CLOGPCLOGP == to_fingerprint) {
    for (int i = 0; i < bit_replicates; i++) {
      sfc.hit_bit(i * 2, int_logp);
      sfc.hit_bit(i * 2 + 1, int_logp);
    }
  } else if (CLOGPCLOGD == to_fingerprint) {
    for (int i = 0; i < bit_replicates; i++) {
      sfc.hit_bit(i * 2, int_logp);
      sfc.hit_bit(i * 2 + 1, int_logd);
    }
  } else if (CLOGDCLOGD == to_fingerprint) {
    for (int i = 0; i < bit_replicates; i++) {
      sfc.hit_bit(i * 2, int_logd);
      sfc.hit_bit(i * 2 + 1, int_logd);
    }
  } else {
    cerr << "Not sure what to fingerprint, cannot continue\n";
    return 0;
  }

  // cerr << "SFC contains " << sfc.nbits() << " bits\n";

  return common_write_fingerprint(m.smiles(), m.name(), sfc, output);
}

static int
write_fingerprint_data(const const_IWSubstring& smiles, const const_IWSubstring& id,
                       const const_IWSubstring& slogp,
                       IWString_and_File_Descriptor& output) {
  float clogp;

  if (!slogp.numeric_value(clogp)) {
    cerr << "Invalid clogp '" << slogp << "', set to zero\n";
    clogp = 0.0F;
  }

  int int_logp = convert_computed_to_positive_int(clogp);

  // cerr << " Convert " << clogp << " to " << int_logp << '\n';

  if (verbose) {
    clogp_fp_histogram[int_logp]++;
  }

  Sparse_Fingerprint_Creator sfc;

  for (int i = 0; i < bit_replicates; i++) {
    sfc.hit_bit(2 * i, int_logp);
  }

  return common_write_fingerprint(smiles, id, sfc, output);
}

static int
logd_final_processing(Molecule& m, float clogp, float offset,
                      IWString_and_File_Descriptor& output) {
  if (tag.length()) {
    return write_fingerprint_data(m, clogp, offset, output);
  }

  float clogd = clogp - offset;

  if (!filter_clogd_values) {
    if (truncate_to_first_token_of_name) {
      write_first_token_of_identifier(m.name(), output);
    } else {
      output << m.name();
    }

    output << clogp << ' ' << clogd << ' ' << offset << '\n';

    if (flush_after_every_molecule) {
      output.flush();
    }

    return 1;
  }

  if (clogd < min_clogd) {
    discarded_for_clogd_too_low++;
    return rejected_by_logd(m, clogp, clogd);
  } else if (clogd > max_clogd) {
    discarded_for_clogd_too_high++;
    return rejected_by_logd(m, clogp, clogd);
  }

  output << m.smiles() << ' ' << m.name();

  if (append_clogp_to_filter_output) {
    output << ' ' << clogd;
  }

  output << '\n';

  return 1;
}

static double
do_rerun_charge_assigner_queries(Molecule& m, Charge_Assigner& charge_assigner,
                                 formal_charge_t* fc) {
  int matoms = m.natoms();

  std::fill_n(fc, matoms, 0);

  int n = charge_assigner.number_elements();

  double rc = 0.0;

  Molecule_to_Match target(&m);

  for (int i = 0; i < n; i++) {
    Substructure_Results sresults;

    Substructure_Hit_Statistics* q = charge_assigner[i];

    int nhits = q->substructure_search(target, sresults);

    if (0 == nhits) {
      continue;
    }

    double logd_offset;
    if (!q->numeric_value(logd_offset)) {  // hmmm, no logd offset!!
      continue;
    }

    //  cerr << nhits << " to query '" << q->comment() << "' offset " << logd_offset <<
    //  '\n';

    for (int j = 0; j < nhits; j++) {
      const Query_Atoms_Matched* qam = sresults.query_atoms_matching(j);

      const Set_of_Atoms* e = sresults.embedding(j);

      for (int k = 0; k < qam->number_elements(); k++) {
        atom_number_t l = e->item(k);

        if (0 != fc[l]) {  // already hit by this or something else
          continue;
        }

        const Substructure_Atom* a = qam->item(k);

        double charge_specification;
        if (!a->numeric_value(charge_specification)) {  // no charge placed here
          continue;
        }

        rc += logd_offset;
        if (charge_specification < 0.0) {
          fc[l] = -1;
        } else {
          fc[l] = 1;
        }
      }
    }
  }

  return rc;
}

#define NEUTRAL_UNCHARGED 0.34
#define QUAT -0.29
#define MULTI_MINUS 3.09
#define POSITIVE3 3.38
#define POSITIVE2 2.14
#define ZWIT_1_MINUS 0.71
#define ZWIT_2_MINUS 2.10

static int
clogp2descriptors(Molecule& m, float clogp, IWString_and_File_Descriptor& output) {
  preprocess(m);

  // First any existing formal charge. We are assuming that the only
  // kind of formal charge that could be present is +ve

  if (m.has_formal_charges()) {
    //  append_result(output, clogp - QUAT, QUAT);
    logd_stats.quats_encountered++;
    return logd_final_processing(m, clogp, QUAT, output);
  }

  int matoms = m.natoms();

  formal_charge_t* fc = new formal_charge_t[matoms];
  std::unique_ptr<formal_charge_t[]> free_fc(fc);
  set_vector(fc, matoms, 0);

  if (0 == charge_assigner.process(m, fc)) {
    //  append_result(output, (clogp - NEUTRAL_UNCHARGED), NEUTRAL_UNCHARGED);
    logd_stats.no_charge_encountered++;
    return logd_final_processing(m, clogp, NEUTRAL_UNCHARGED, output);
  }

  int npos = 0;
  int nneg = 0;

  for (int i = 0; i < matoms; i++) {
    formal_charge_t fci = fc[i];

    if (0 == fci) {
      ;
    } else if (fci > 0) {
      npos++;
    } else {
      nneg++;
    }
  }

  if (0 == npos && 0 == nneg)  // really should have been picked up before
  {
    //  append_result(output, (clogp - NEUTRAL_UNCHARGED), NEUTRAL_UNCHARGED);
    logd_stats.no_charge_encountered++;
    return logd_final_processing(m, clogp, NEUTRAL_UNCHARGED, output);
  }

  if (nneg > 1 && 0 == npos) {
    //  append_result(output, (clogp - MULTI_MINUS), MULTI_MINUS);
    logd_stats.multi_minus++;
    return logd_final_processing(m, clogp, MULTI_MINUS, output);
  }

  if (npos >= 3 && 0 == nneg) {
    //  append_result(output, (clogp - POSITIVE3), POSITIVE3);
    logd_stats.positive3++;
    return logd_final_processing(m, clogp, POSITIVE3, output);
  }

  if (2 == npos && 0 == nneg) {
    //  append_result(output, (clogp - POSITIVE2), POSITIVE2);
    logd_stats.positive2++;
    return logd_final_processing(m, clogp, POSITIVE2, output);
  }

  if (npos > 0 && 1 == nneg) {
    //  append_result(output, (clogp - ZWIT_1_MINUS), ZWIT_1_MINUS);
    logd_stats.zwit_1_minus++;
    return logd_final_processing(m, clogp, ZWIT_1_MINUS, output);
  }

  if (npos > 0 && nneg > 1) {
    //  append_result(output, (clogp - ZWIT_2_MINUS), ZWIT_2_MINUS);
    logd_stats.zwit_2_minus++;
    return logd_final_processing(m, clogp, ZWIT_2_MINUS, output);
  }

  if (npos > 0 && 0 == nneg) {
    ;
  } else if (0 == npos && nneg > 0) {
    ;
  } else {
    cerr << "clogp2descriptors:unusual charge state " << npos << "+ and " << nneg
         << "-\n";
    return 0;
  }

  // Doesn't quality as any of the special cases. Compute the offset from
  // the queries that match

  if (rerun_charge_assigner_queries) {
    float offset = do_rerun_charge_assigner_queries(m, charge_assigner, fc);
    cerr << "Offset from rerun " << offset << '\n';
    //  append_result(output, (clogp - offset), offset);
    return logd_final_processing(m, clogp, offset, output);
  }

  double offset = static_cast<double>(0.0);

  for (int i = 0; i < matoms; i++) {
    formal_charge_t fci = fc[i];
    if (0 == fci) {
      continue;
    }

    if (fci < 0) {
      fci = -fci;
    }

    int query_number = fci / CHARGE_ASSIGNER_MULTIPLIER;

    //  cerr << "Atom " << i << ' ' << m.smarts_equivalent_for_atom(i) <<  " matched to
    //  query number " << query_number << '\n';

    double d;
    if (!charge_assigner[query_number]->numeric_value(d)) {
      continue;
    }

    offset += static_cast<double>(d);
  }

  // append_result(output, (clogp - offset), offset);

  return logd_final_processing(m, clogp, offset, output);
}

static int
common_write_smiles(const const_IWSubstring& smiles, const const_IWSubstring& id,
                    const const_IWSubstring& sclogp, int append_sclogp,
                    IWString_and_File_Descriptor& output) {
  output << smiles << ' ' << id;

  if (append_sclogp) {
    output << ' ' << sclogp;
  }

  output << '\n';

  output.write_if_buffer_holds_more_than(32768);

  return 1;
}

static int
write_to_rejected_file_if_needed(const const_IWSubstring& smiles,
                                 const const_IWSubstring& id,
                                 const const_IWSubstring& sclogp) {
  if (!stream_for_discarded_by_filter.is_open()) {
    return 1;
  }

  return common_write_smiles(smiles, id, sclogp, append_clogp_to_rejected_output,
                             stream_for_discarded_by_filter);
}

static int
common_rejection_handling(const const_IWSubstring& smiles, const const_IWSubstring& id,
                          const const_IWSubstring& sclogp) {
  if (filter_clogp_values) {
    return write_to_rejected_file_if_needed(smiles, id, sclogp);
  }

  return 1;
}

static int
_clogp2descriptors_record(const const_IWSubstring& buffer, int id_must_be_present,
                          IWString_and_File_Descriptor& output) {
  const_IWSubstring sclogp, serrlvl, smiles;
  IWString id;
  int i = 0;

  if (!buffer.nextword(sclogp, i) || !buffer.nextword(serrlvl, i) ||
      !buffer.nextword(smiles, i)) {
    cerr << "clogp2descriptors_record:not enough tokens\n";
    return 0;
  }

  if (id_must_be_present && !buffer.nextword(id, i)) {
    cerr << "clogp2descriptors_record:identifier missing\n";
    return 0;
  }

  if (!truncate_to_first_token_of_name) {
    const_IWSubstring tmp;
    while (buffer.nextword(tmp, i)) {
      id << ' ' << tmp;
    }
  }

  // A negative error level indicates total failure

  if ('-' == serrlvl[0]) {
    cerr << "clogp failed '" << id << "'\n";
    return common_rejection_handling(smiles, id, missing_values);
  }

  if (clogp_max_error_level <= 60)  // need to interpret serrlvl and check
  {
    int errlvl;

    if (!serrlvl.numeric_value(errlvl) || errlvl < 0) {
      cerr << "clogp2descriptors_record:invalid errlvl\n";
      return 0;
    }

    if (errlvl > clogp_max_error_level) {
      molecules_above_clogp_max_error_level++;

      return common_rejection_handling(smiles, id, sclogp);
    }
  }

  // errlvl OK, what next

  if (filter_clogp_values) {
    float clogp;
    if (!sclogp.numeric_value(clogp)) {
      cerr << "Non numeric logp value '" << sclogp << "'\n";
      return 0;
    }

    //  cerr << "Within filter, value " << clogp << '\n';

    if (clogp < min_clogp) {
      discarded_for_clogp_too_low++;
      return common_rejection_handling(smiles, id, sclogp);
    } else if (clogp > max_clogp) {
      discarded_for_clogp_too_high++;
      return common_rejection_handling(smiles, id, sclogp);
    } else {
      return common_write_smiles(smiles, id, sclogp, append_clogp_to_filter_output,
                                 output);
    }
  }

  if (!charge_assigner.active()) {
    if (tag.length()) {
      return write_fingerprint_data(smiles, id, sclogp, output);
    }

    if (write_smiles_file) {
      output << smiles << ' ';
    }

    write_first_token_of_identifier(id, output);
    output << sclogp << '\n';

    return 1;
  }

  Molecule m;

  if (!m.build_from_smiles(smiles)) {
    cerr << "clogp2descriptors_record:invalid smiles '" << smiles << "'\n";
    return 0;
  }

  m.set_name(id);

  float clogp;

  if (!sclogp.numeric_value(clogp)) {
    cerr << "clogp2descriptors_record:invalid numeric for clogp\n";
    return 0;
  }

  return clogp2descriptors(m, clogp, output);
}

static int
clogp2descriptors_record(const const_IWSubstring& buffer, int id_must_be_present,
                         IWString_and_File_Descriptor& output) {
  if (_clogp2descriptors_record(buffer, id_must_be_present, output)) {
    return 1;
  }

  return fault_tolerant;
}

static int
clogp2descriptors(iwstring_data_source& input, IWString_and_File_Descriptor& output) {
  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    molecules_read++;

    if (!clogp2descriptors_record(buffer, 1, output)) {
      cerr << "clogp2descriptors:fatal error '" << buffer << "'\n";
      return 0;
    }

    output.write_if_buffer_holds_more_than(32768);
  }

  output.flush();

  return 1;
}

// Biobyte seems to have an internal maximum record lenght, so it will wrap records.
// Probably need to check that behaviour when new versions come out

static int
clogp2descriptors_fingerprint_sent_to_clogp(iwstring_data_source& input,
                                            IWString_and_File_Descriptor& output) {
  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    if (buffer.starts_with("   0.00 -1 "))  // not processed by clogp
    {
      buffer.remove_leading_chars(11);
      output << buffer;
      output.write_if_buffer_holds_more_than(4096);

      if ('|' == buffer) {
        ;
      } else {
        while (!buffer.ends_with('>')) {
          //        cerr << "Reading continuation record " << input.lines_read() << '\n';
          if (!input.next_record(buffer)) {
            cerr << "GACK, unterminated TDT completion not found\n";
            return 0;
          }

          buffer.remove_leading_chars(11);
          output << buffer;
          output.write_if_buffer_holds_more_than(4096);
        }
      }

      output << '\n';

      continue;
    }

    if (!clogp2descriptors_record(buffer, 0, output)) {
      cerr << "clogp2descriptors::fatal error '" << buffer << "'\n";
      return 0;
    }

    output.write_if_buffer_holds_more_than(4096);
  }

  output.flush();

  return 1;
}

static int
clogp2descriptors_fingerprint_sent_to_clogp(const char* fname,
                                            IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return clogp2descriptors_fingerprint_sent_to_clogp(input, output);
}

static int
do_parent_process(int pipe_to_subprocess[2], int pipe_from_subprocess[2], pid_t child,
                  iwstring_data_source& input, IWString_and_File_Descriptor& output) {
  close(pipe_to_subprocess[0]);    // we will be writing to the sub process
  close(pipe_from_subprocess[1]);  // we read from the sub process via this pipe

  char readbuf[4096];

  IWString buffer;

  cerr << "Parent working with child " << child << '\n';

  while (input.next_record(buffer)) {
    output << buffer << '\n';
    output.write_if_buffer_holds_more_than(4096);

    if (!buffer.starts_with(smiles_tag)) {
      continue;
    }

    buffer.remove_up_to_first('<');
    buffer.chop();
    buffer << " Q\n";

    int bytes_written = write(pipe_to_subprocess[1], buffer.rawchars(), buffer.length());
    if (bytes_written != buffer.length()) {
      cerr << "Yipes, cannot write to biobyte pipe\n";
      return 0;
    }
    cerr << "Parent wrote " << bytes_written << " bytes to pipe\n";

    size_t bytes_read = read(pipe_from_subprocess[0], readbuf, sizeof(readbuf));

    output.strncat(readbuf, static_cast<int>(bytes_read));
  }

  output.flush();

  close(pipe_to_subprocess[0]);
  close(pipe_to_subprocess[1]);
  close(pipe_from_subprocess[0]);
  close(pipe_from_subprocess[1]);

  wait(NULL);  // only one child created here

  return 1;
}

static int
do_child_process(const char* biobyte_cmd, int pipe_from_parent[2],
                 int pipe_to_parent[2]) {
  close(pipe_from_parent[1]);  // close write end of pipe from parent
  close(pipe_to_parent[0]);    // close read  end of pipe to   parent

  int d = dup2(pipe_from_parent[0], 0);  // my stdin is coming from parent
  if (0 != d) {
    cerr << "Dup2 of stdin failed in child\n";
    exit(1);
  }

  d = dup2(pipe_to_parent[1], 1);  // my stdout goes to parent
  if (1 != d) {
    cerr << "Dup2 of stdout failed in child\n";
    exit(1);
  }

  close(pipe_from_parent[0]);  // now duplicated, close original
  close(pipe_to_parent[1]);    // now duplicated, close original

#ifdef CHILD_DOES_READS
  char buffer[4096];

  while (1) {
    int bytes_read = read(0, buffer, sizeof(buffer));
    if (0 == bytes_read) {
      cerr << "Zero byte read\n";
      break;
    }
    cerr << "Read " << bytes_read << " bytes\n";
    write(1, buffer, bytes_read);
  }

  exit(0);
#endif

  // cerr << "Child going to biobyte '" << biobyte_cmd << "'\n";

  const char* argv[] = {"HELLO", NULL};

  execv(biobyte_cmd, const_cast<char**>(argv));

  cerr << "Child process '" << biobyte_cmd << "' not started, catastrophic failure\n";

  exit(1);
}

static int
clogp2descriptors_tdt(iwstring_data_source& input, const char* biobyte_clogp_cmd,
                      IWString_and_File_Descriptor& output) {
  IWString cmd = biobyte_clogp_cmd;
  int pipe_to_subprocess[2], pipe_from_subprocess[2];

  if (0 != pipe(pipe_to_subprocess) || 0 != pipe(pipe_from_subprocess))  // create pipes
  {
    cerr << "Cannot initialise pipe\n";
    return 0;
  }

  pid_t child = fork();
  if (child < 0) {
    cerr << "Cannot fork sub process for logp computation\n";
    return 0;
  }

  if (0 == child) {
    do_child_process(biobyte_clogp_cmd, pipe_to_subprocess, pipe_from_subprocess);
  } else {
    do_parent_process(pipe_to_subprocess, pipe_from_subprocess, child, input, output);
  }

  return 1;
}

static int
clogp2descriptors_tdt(const char* biobyte_clogp_cmd,
                      IWString_and_File_Descriptor& output) {
  iwstring_data_source input(0);

  return clogp2descriptors_tdt(input, biobyte_clogp_cmd, output);
}

static int
clogp2descriptors(const char* fname, IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return clogp2descriptors(input, output);
}

static int
parse_range(const Command_Line& cl, char lflag, char uflag, float& zmin, float& zmax,
            const char* clogp_or_clogd) {
  if (cl.option_present(lflag)) {
    if (!cl.value(lflag, zmin)) {
      cerr << "The min " << clogp_or_clogd << " value(-" << lflag
           << ") must be a valid floating point number\n";
      return 0;
    }

    if (verbose) {
      cerr << "Will discard molecules with " << clogp_or_clogd << " values < "
           << min_clogp << '\n';
    }
  }

  if (cl.option_present(uflag)) {
    if (!cl.value(uflag, zmax)) {
      cerr << "The max " << clogp_or_clogd
           << " value(-C) must be a valid floating point number\n";
      usage(3);
    }

    if (cl.option_present(lflag) && zmax < zmin) {
      cerr << "Inconsistent min " << zmin << " and max " << zmax << " " << clogp_or_clogd
           << " values\n";
      return 3;
    }

    if (verbose) {
      cerr << "Will discard molecules with " << clogp_or_clogd << " values > " << zmax
           << '\n';
    }
  }

  return 1;
}

static void
DisplayDashXOptions(std::ostream& output) {
  output << " -X flush      flush output after each molecule\n";
  output << " -X smiles     for clogp output, write 'smiles id clogp'\n";

  ::exit(0);
}

static int
clogp2descriptors(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:M:e:N:g:Y:ruc:C:B:abd:D:J:p:k:F:qUX:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('Y')) {
    cl.value('Y', descriptor_prefix);

    if (verbose) {
      cerr << "Will prepend '" << descriptor_prefix
           << "' to the names of all descriptors\n";
    }
  }

  if (cl.option_present('e')) {
    if (!cl.value('e', clogp_max_error_level) || clogp_max_error_level < 0) {
      cerr
          << "The maximum valid error value(-e option) must be a whole positive number\n";
      usage(13);
    }

    if (verbose) {
      cerr << "Error levels above " << clogp_max_error_level
           << " will be considered invalid\n";
    }
  }

  if (cl.option_present('u')) {
    fault_tolerant = 1;

    if (verbose) {
      cerr << "Will run in fault tolerant mode\n";
    }
  }

  if (cl.option_present('c') || cl.option_present('C')) {
    if (!parse_range(cl, 'c', 'C', min_clogp, max_clogp, "clogp")) {
      cerr << "Cannot parse min and max clogp specification(s)\n";
      usage(3);
    }

    filter_clogp_values = 1;
  }

  if (cl.option_present('d') || cl.option_present('D')) {
    if (!parse_range(cl, 'd', 'D', min_clogd, max_clogd, "clogd")) {
      cerr << "Cannot parse min and max clogd specification(s)\n";
      usage(3);
    }

    filter_clogd_values = 1;
  }

  if (filter_clogp_values && filter_clogd_values) {
    cerr << "Sorry, cannot filter both clogp and clogd\n";
    usage(3);
  }

  if (cl.option_present('B')) {
    IWString b = cl.option_value('B');

    if (!b.ends_with(".smi")) {
      b << ".smi";
    }

    if (!stream_for_discarded_by_filter.open(b.null_terminated_chars())) {
      cerr << "Cannot open stream for discarded molecules '" << b << "'\n";
      return 3;
    }

    if (verbose) {
      cerr << "Discarded molecules written to '" << b << "'\n";
    }
  }

  if (cl.option_present('a')) {
    append_clogp_to_filter_output = 1;

    if (verbose) {
      cerr << "Will append clogp to filtered output\n";
    }
  }

  if (cl.option_present('b')) {
    append_clogp_to_rejected_output = 1;

    if (verbose) {
      cerr << "Will append clogp to rejected output\n";
    }
  }

  if (cl.option_present('J')) {
    cl.value('J', tag);

    if (!tag.ends_with('<')) {
      tag << '<';
    }

    if (verbose) {
      cerr << "logP and logD written as fingerprints\n";
    }

    if (cl.option_present('p')) {
      if (!cl.value('p', bit_replicates) || bit_replicates <= 0) {
        cerr << "The bit replicates value(-p) must be a whole +ve number\n";
        usage(2);
      }

      if (verbose) {
        cerr << "Fingerprint bit replicates set to " << bit_replicates << '\n';
      }
    }

    if (cl.option_present('k')) {
      const_IWSubstring k = cl.string_value('k');

      if ("pp" == k) {
        to_fingerprint = CLOGPCLOGP;
      } else if ("pd" == k) {
        if (!cl.option_present('N')) {
          cerr
              << "Request for 'pd' output, but no charge assigner for logD computation\n";
          return 3;
        }
        to_fingerprint = CLOGPCLOGD;
        if ("NCLP<" == tag) {  // too difficult to be properly flexible within gfp_make
          tag = "NCLPD<";
        }
      } else if ("dd" == k) {
        to_fingerprint = CLOGDCLOGD;
        if (!cl.option_present('N')) {
          cerr
              << "Request for 'pd' output, but no charge assigner for logD computation\n";
          return 3;
        }
        if ("NCLP<" == tag) {  // too difficult to be properly flexible within gfp_make
          tag = "NCLD<";
        }
      } else {
        cerr << "Unrecognised -k qualifier '" << k << "'\n";
        usage(2);
      }
    }
  }

  if (tag.length() &&
      CLOGPCLOGP == to_fingerprint) {  // no need for clogD, so charge assigner not needed
    ;
  } else if (cl.option_present('N'))  // doing chemically meaningful things
  {
    if (!process_standard_aromaticity_options(cl)) {
      cerr << "Cannot process -A options\n";
      usage(8);
    }

    if (!charge_assigner.construct_from_command_line(cl, verbose, 'N')) {
      cerr << "Cannot initialise charge assigner\n";
      usage(7);
    }

    charge_assigner.set_apply_charges_to_molecule(0);
    charge_assigner.set_assigned_charge_multiplier(CHARGE_ASSIGNER_MULTIPLIER);

    if (cl.option_present('g')) {
      if (!chemical_standardisation.construct_from_command_line(cl, (verbose > 1), 'g')) {
        cerr << "Cannot initialise chemical standardisation (-g)\n";
        usage(6);
      }
    }

    if (cl.option_present('r')) {
      rerun_charge_assigner_queries = 1;
      if (verbose) {
        cerr << "Will re-run charge assigner queries\n";
      }
    }
  }

  if (cl.option_present('q')) {
    work_as_tdt_filter = 1;

    if (verbose) {
      cerr << "Will work as a TDT filter\n";
    }
  }

  if (cl.option_present('U')) {
    truncate_to_first_token_of_name = 0;

    if (verbose) {
      cerr << "Will NOT truncate to the first token of the name\n";
    }
  }

  // Must do missing values after charge assigner

  IWString m = '.';

  if (cl.option_present('M')) {
    IWString m = cl.string_value('M');

    if (verbose) {
      cerr << "Missing values written as '" << m << "'\n";
    }
  }

  missing_values = m;
  if (charge_assigner.active()) {
    missing_values << ' ' << m << ' ' << m;
  }

  if (cl.option_present('X')) {
    const_IWSubstring x;
    for (int i = 0; cl.value('X', x, i); ++i) {
      if (x == "flush") {
        flush_after_every_molecule = 1;
        if (verbose) {
          cerr << "Will flush after every molecule\n";
        }
      } else if (x == "smiles") {
        write_smiles_file = 1;
        if (verbose) {
          cerr << "Will write a smiles file with the computed value appended to the "
                  "name\n";
        }
      } else if (x == "help") {
        DisplayDashXOptions(cerr);
      } else {
        cerr << "Unrecognised -X qualifier '" << x << "'\n";
        DisplayDashXOptions(cerr);
      }
    }
  }

  IWString_and_File_Descriptor output(1);

  if (filter_clogd_values || filter_clogp_values) {
    ;
  } else if (tag.length()) {
    ;
  } else if (cl.option_present('F')) {
    ;
  } else if (cl.option_present('q')) {
    ;
  } else if (write_smiles_file) {
  } else {
    output << "Name " << descriptor_prefix << "clogp";

    if (charge_assigner.active()) {
      output << ' ' << descriptor_prefix << "brnsclogD " << descriptor_prefix
             << "brnsclogDoffset";
    }

    output << '\n';
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(1);
  }

  int rc = 0;

  if (cl.option_present(
          'F'))  // never could get this to work, not sure why, it is incomplete..
  {
    const char* biobyte_clogp_cmd = cl.option_value('F');

    if (!clogp2descriptors_tdt(biobyte_clogp_cmd, output)) {
      rc = 1;
    }
  } else if (cl.option_present('q')) {
    rc = clogp2descriptors_fingerprint_sent_to_clogp(cl[0], output);
  } else {
    for (int i = 0; i < cl.number_elements(); i++) {
      if (!clogp2descriptors(cl[i], output)) {
        rc = i + 1;
        break;
      }
    }
  }

  if (stream_for_discarded_by_filter.is_open()) {
    stream_for_discarded_by_filter.flush();
  }

  if (verbose) {
    cerr << "Processed " << molecules_read << " molecules's\n";
    cerr << molecules_above_clogp_max_error_level
         << " molecules above max clogp error level " << clogp_max_error_level << '\n';
  }

  if (verbose && filter_clogp_values) {
    if (-std::numeric_limits<float>::max() != min_clogp) {
      cerr << "Discarded " << discarded_for_clogp_too_low << " molecules for clogp < "
           << min_clogp << '\n';
    }
    if (std::numeric_limits<float>::max() != max_clogp) {
      cerr << "Discarded " << discarded_for_clogp_too_high << " molecules for clogp > "
           << max_clogp << '\n';
    }
  }

  if (verbose && filter_clogd_values) {
    if (-std::numeric_limits<float>::max() != min_clogd) {
      cerr << "Discarded " << discarded_for_clogd_too_low << " molecules for clogd < "
           << min_clogd << '\n';
    }
    if (std::numeric_limits<float>::max() != max_clogd) {
      cerr << "Discarded " << discarded_for_clogd_too_high << " molecules for clogd > "
           << max_clogd << '\n';
    }
  }

  if (verbose && tag.length()) {
    cerr << clogp_fp_histogram.number_elements()
         << " different clogp fingerprint values\n";
    for (int i = 0; i < clogp_fp_histogram.number_elements(); i++) {
      if (clogp_fp_histogram[i]) {
        cerr << clogp_fp_histogram[i] << " molecules have value " << i << '\n';
      }
    }
  }

  return rc;
}

int
main(int argc, char** argv) {
  int rc = clogp2descriptors(argc, argv);

  return rc;
}
