#include <algorithm>
#include <iostream>
#include <memory>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/element.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"

#include "Molecule_Tools/iwecfp_lib.h"

namespace iwecfp {

using std::cerr;

class Options {
  private:
    // Each call to the generator needs an array of Sparse_Fingerprint_Creator.
    // Allocate once and then clear before each call.
    int _nsfc;
    std::unique_ptr<Sparse_Fingerprint_Creator[]> _sfc;

    int _verbose;

    int _function_as_tdt_filter;

    int _flush_after_each_molecule;

    int _fixed_width_fingerprint;

    int _bit_count_multiplier;

    int _bit_replicates;
    int _bit_replicate_offset;

    int _descriptor_file_output;

    int _max_shell_radius;

    int _each_shell_gets_own_fingerprint;

    // For fingerprint output.
    IWString _tag;

    IWString _smiles_tag;
    IWString _identifier_tag;

    Accumulator_Int<uint32_t> _nbits_acc;

    // If writing a descriptor file, a write buffer.
    std::unique_ptr<int[]> _count;

    // Preprocessing.
    int _reduce_to_largest_fragment;

    Chemical_Standardisation _chemical_standardisation;

    Atom_Typing_Specification _atom_typing;
    
    uint64_t _molecules_read = 0;

  // private functions
    int ParseReplicatesDirective(const const_IWSubstring& y);

    int WriteFixedWidthFingerprint(Molecule& m, uint32_t nbits,
                                    const Sparse_Fingerprint_Creator& sfc,
                                    IWString_and_File_Descriptor& output);
    int WriteDescriptorFileRow(Molecule& m, uint32_t ncols,
                                const Sparse_Fingerprint_Creator& sfc,
                                IWString_and_File_Descriptor& output);
    int WriteArrayOfFingerprintsInner(Sparse_Fingerprint_Creator* sfc,
                                       IWString_and_File_Descriptor& output);
    int WriteArrayOfFingerprints(Sparse_Fingerprint_Creator* sfc,
                                 IWString_and_File_Descriptor& output);
  public:
    Options();
    ~Options();

    int function_as_tdt_filter() const {
      return _function_as_tdt_filter;
    }

    int descriptor_file_output() const {
      return _descriptor_file_output;
    }

    const IWString& smiles_tag() const {
      return _smiles_tag;
    }

    int Initialise(Command_Line& cl);

    int Fingerprint(Molecule& m, iwecfp::Iwecfp& generator,
                     IWString_and_File_Descriptor& output);
    // Only does this if descriptor file output is requested.
    int WriteDescriptorFileHeader(IWString_and_File_Descriptor& output) const;

    int Preprocess(Molecule& m);

    // If queries are present, what do we do withb
    int HandleNoStartAtoms(Molecule& m, IWString_and_File_Descriptor& output);

    int DoOutput(Molecule& m, Sparse_Fingerprint_Creator* sfc,
                 IWString_and_File_Descriptor& output);

    void MaybeFlush(IWString_and_File_Descriptor& output);

    int Report(std::ostream& output) const;
};

Options::Options() {
  _nsfc = 0;
  _verbose = 0;
  _function_as_tdt_filter = 0;
  _flush_after_each_molecule = 0;
  _fixed_width_fingerprint = 0;
  _descriptor_file_output = 0;
  _max_shell_radius = 0;
  _bit_count_multiplier = 0;
  _bit_replicates = 0;
  _bit_replicate_offset = 741;
  _each_shell_gets_own_fingerprint = 0;
  _reduce_to_largest_fragment = 0;
  _smiles_tag = "$SMI<";
  _identifier_tag = "PCN<";
}

Options::~Options() {
}

void
DisplayDashYOptions(std::ostream& output) {
  output << " -Y flush flush output after each molecule\n";
  output << " -Y fixed=<nbits> generate fixed width fingerprints with <nbits> bits\n";
  output << " -Y desc=<ncols> generate descriptor file output with <ncols> columns\n";
  output << " -Y scale=<n> multiple all bit counts by <n>. Beware truncation at 255.\n";
  output << " -Y replicates=<n>[,<offset>] replicate bits with optional offset\n";
  output << " -Y pchiral include central atom isotope offset in generated bits\n";
  ::exit(0);
}

int
Options::ParseReplicatesDirective(const const_IWSubstring& y) {
  if (y.empty()) {
    return 0;
  }

  int i = 0;
  const_IWSubstring r, o;
  y.nextword(r, i, ',');
  y.nextword(o, i, ',');

  if (r.empty()) {
    return 0;
  }

  if (!r.numeric_value(_bit_replicates) || _bit_replicates < 1) {
    cerr << "The number of bit replicates must be a positive whole number\n";
    return 0;
  }

  if (o.empty()) {
    return 1;
  }

  if (!o.numeric_value(_bit_replicate_offset) || _bit_replicate_offset == 0) {
    cerr << "Invalid bit replicate offset '" << o << "'\n";
    return 0;
  }

  return 1;
}


int
Options::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');

  if (cl.option_present('f')) {
    _function_as_tdt_filter = 1;
    if (_verbose) {
      cerr << "Will function as a tdt filter\n";
    }
  }

  if (cl.option_present('P')) {
    const_IWSubstring p = cl.string_value('P');
    if (!_atom_typing.build(p)) {
      cerr << "Cannot discern atom type '" << p << "'\n";
      return 0;
    }
  }

  if (cl.option_present('R')) {
    cl.value('R', _max_shell_radius);
  }

  if (cl.option_present('Y')) {
    const_IWSubstring y;
    for (int i = 0; cl.value('Y', y, i); ++i) {
      if (y == "flush") {
        _flush_after_each_molecule = 1;
        if (_verbose) {
          cerr << "Will flush output after each molecule\n";
        }
      } else if (y.starts_with("desc=")) {
        y.remove_leading_chars(5);
        if (!y.numeric_value(_descriptor_file_output) || _descriptor_file_output < 2) {
          cerr << "The number of columns for descriptor output must be a whole +ve number\n";
          return 0;
        }
        if (_verbose) {
          cerr << "Will write a descriptor file with " << _descriptor_file_output << " columns\n";
        }
        _count.reset(new int[_descriptor_file_output]);
      } else if (y.starts_with("fixed=")) {
        y.remove_leading_chars(6);
        if (!y.numeric_value(_fixed_width_fingerprint) || _fixed_width_fingerprint < 2) {
          cerr << "The number of bits for a fixed width fingerprint must be a whole +ve number\n";
          return 0;
        }
        if (_verbose) {
          cerr << "Will write a fixed width fingerprint with "
               << _fixed_width_fingerprint << " bits\n";
        }
      } else if (y.starts_with("scale=")) {
        y.remove_leading_chars(6);
        if (!y.numeric_value(_bit_count_multiplier) || _bit_count_multiplier < 2) {
          cerr << "The bit count scale= directive must be a whole +ve number > 1, '"
               << y << "' invalid\n";
          return 0;
        }
        if (_verbose) {
          cerr << "Bit counts multipled by " << _bit_count_multiplier << '\n';
        }
      } else if (y.starts_with("replicates=")) {
        y.remove_leading_chars(11);
        if (!ParseReplicatesDirective(y)) {
          cerr << "Invalid replicates directive '" << y << "'\n";
          return 0;
        }
//    } else if (y == "pchiral") {
//      _central_atom_possible_chiral = 1;
//      if (_verbose) {
//        cerr << "Central atoms will only be CD4 or CD3H types\n";
//      }
      } else if (y == "help") {
        DisplayDashYOptions(cerr);
      } else {
        cerr << "Unrecongised -Y qualifier '" << y << "'\n";
        DisplayDashYOptions(cerr);
      }
    }
  }

  if (cl.option_present('s')) {
    if (!cl.option_present('R')) {
      cerr << "Sorry, the -s option only works when the -R option is specified\n";
      Usage(1);
    }

    _each_shell_gets_own_fingerprint = 1;
  }

  if (_each_shell_gets_own_fingerprint) {
    _sfc.reset(new Sparse_Fingerprint_Creator[_max_shell_radius + 1]);
    _nsfc = _max_shell_radius + 1;
  } else {
    _sfc.reset(new Sparse_Fingerprint_Creator[1]);
    _nsfc = 1;
  }

  if (_descriptor_file_output && _function_as_tdt_filter) {
    cerr << "Cannot generate descriptors while working as a TDT filter\n";
    return 0;
  }

  if (cl.option_present('l')) {
    _reduce_to_largest_fragment = 1;
    if (_verbose) {
      cerr << "Will reduce to largest fragment\n";
    }
  }

  if (cl.option_present('g')) {
    if (!_chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
      cerr << "Cannot initialise chemical standardisation (-g)\n";
      Usage(14);
    }
  }

  if (_descriptor_file_output) {
    // Descriptor output does not need a fingerprint tag.
  } else if (!cl.option_present('J') && !cl.option_present('P')) {
    _atom_typing.set_atom_type(IWATTYPE_COMPLEX);
    if (_fixed_width_fingerprint) {
      _tag = "FPECC<";
    } else {
      _tag = "NCECC<";
    }
  } else if (!cl.option_present('P') && cl.option_present('J')) {
    cerr << "Must specify atom typing to use with the -P option\n";
    Usage(3);
  } else {
    _atom_typing.swap_atomic_number_atom_type_to_atomic_number_prime();

    if (cl.option_present('J')) {
      _tag = cl.string_value('J');
    } else {
      if (_fixed_width_fingerprint) {
        _tag = "FPEC";
      } else {
        _tag = "NCEC";
      }
      _atom_typing.append_to_tag(_tag);
    }

    if (!_tag.ends_with('<')) {
      _tag += '<';
    }

    if (_verbose) {
      cerr << "Extended connectivity index written as non-colliding sparse fingerprints, tag '"
           << _tag << "'\n";
    }
  }

  if (_each_shell_gets_own_fingerprint) {
    _tag.chop();
  }

  return 1;
}

int
Options::WriteDescriptorFileHeader(IWString_and_File_Descriptor& output) const {
  if (!_descriptor_file_output) {
    return 1;
  }

  output << "Id";
  for (int i = 0; i < _descriptor_file_output; ++i) {
    output << " EC" << _max_shell_radius << '_' << i;
  }
  output << '\n';

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Generated fingerprints for " << _molecules_read << " molecules\n";
  if (_nbits_acc.n() > 0) {
    output << "Fingerprints had between " << _nbits_acc.minval() << " and "
           << _nbits_acc.maxval() << " ave "
           << static_cast<float>(_nbits_acc.average_if_available_minval_if_not())
           << " bits set\n";
  }

  return 1;
}

int
Options::Preprocess(Molecule& m) {
  m.remove_all(1);

  if (_reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment();
  }

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }

  m.compute_aromaticity_if_needed();

  return 1;
}

void
Options::MaybeFlush(IWString_and_File_Descriptor& output) {
  if (_flush_after_each_molecule) {
    output.flush();
  } else {
    output.write_if_buffer_holds_more_than(32768);
  }
}

int
Options::WriteArrayOfFingerprintsInner(Sparse_Fingerprint_Creator* sfc,
                                       IWString_and_File_Descriptor& output) {
  IWString tmp;

  if (_bit_count_multiplier > 0) {
    sfc[0].Scale(_bit_count_multiplier);
  }

  if (_each_shell_gets_own_fingerprint) {
    for (int i = 0; i <= _max_shell_radius; ++i) {
      IWString tag_this_fp(_tag);
      tag_this_fp << i;
      tmp.resize_keep_storage(0);
      sfc[i].daylight_ascii_form_with_counts_encoded(tag_this_fp, tmp);
      output << tmp << '\n';
    }
  } else {
    if (_tag.starts_with("FP")) {
      static constexpr int kNbits = 2048;
      IW_Bits_Base bits(kNbits);
      for (const auto& [bit, count] : sfc[0].bits_found()) {
        (void) count;
        bits.set_bit(bit % kNbits);
      }
      bits.daylight_ascii_representation_including_nset_info(tmp);
      output << _tag << tmp << ">\n";
    } else {
      sfc[0].daylight_ascii_form_with_counts_encoded(_tag, tmp);
      output << tmp << '\n';
    }
  }

  return output.size();
}

int
Options::WriteArrayOfFingerprints(Sparse_Fingerprint_Creator* sfc,
                                  IWString_and_File_Descriptor& output) {
  if (_bit_replicates == 0) {
    WriteArrayOfFingerprintsInner(sfc, output);
    return 1;
  }

  if (_each_shell_gets_own_fingerprint) {
    return 1;  // Preserve original lazy behaviour.
  }

  Sparse_Fingerprint_Creator tmp;
  for (const auto& [b, c] : sfc[0].bits_found()) {
    tmp.hit_bit(b, c);
  }

  if (_bit_replicates > 0) {
    for (const auto& [b, c] : sfc[0].bits_found()) {
      for (int i = 1; i <= _bit_replicates; ++i) {
        const uint32_t newbit = b + i * _bit_replicate_offset;
        tmp.hit_bit(newbit, c);
      }
    }
  }

  WriteArrayOfFingerprintsInner(&tmp, output);

  return 1;
}

int
Options::WriteFixedWidthFingerprint(Molecule& m, uint32_t nbits,
                                    const Sparse_Fingerprint_Creator& sfc,
                                    IWString_and_File_Descriptor& output) {
  output << _smiles_tag << m.smiles() << ">\n";
  sfc.write_constant_width_fingerprint(nbits, _tag, output);
  output << _identifier_tag << m.name() << ">\n";
  output << "|\n";

  return 1;
}

int
Options::WriteDescriptorFileRow(Molecule& m, uint32_t ncols,
                                const Sparse_Fingerprint_Creator& sfc,
                                IWString_and_File_Descriptor& output) {
  append_first_token_of_name(m.name(), output);

  std::fill_n(_count.get(), ncols, 0);

  sfc.WriteAsDescriptors(ncols, _count.get(), output);
  output << '\n';

  return 1;
}

int
Options::DoOutput(Molecule& m, Sparse_Fingerprint_Creator* sfc,
                  IWString_and_File_Descriptor& output) {
  if (_verbose) {
    _nbits_acc.extra(sfc[0].nbits());
  }

  if (_function_as_tdt_filter) {
    WriteArrayOfFingerprints(sfc, output);
  } else if (_fixed_width_fingerprint) {
    WriteFixedWidthFingerprint(m, _fixed_width_fingerprint, sfc[0], output);
  } else if (_descriptor_file_output) {
    WriteDescriptorFileRow(m, _descriptor_file_output, sfc[0], output);
  } else {
    output << _smiles_tag << m.smiles() << ">\n";
    WriteArrayOfFingerprints(sfc, output);
    output << _identifier_tag << m.name() << ">\n";
    output << "|\n";
  }

  MaybeFlush(output);

  return output.good();
}

int
Options::HandleNoStartAtoms(Molecule& m, IWString_and_File_Descriptor& output) {
  if (!_function_as_tdt_filter) {
    output << _smiles_tag << m.smiles() << ">\n";
    output << _identifier_tag << m.name() << ">\n";
  }

  output << _tag << ">\n";

  if (!_function_as_tdt_filter) {
    output << "|\n";
  }

  output.write_if_buffer_holds_more_than(4096);
  return 1;
}

int
Options::Fingerprint(Molecule& m, iwecfp::Iwecfp& generator,
                     IWString_and_File_Descriptor& output) {
  ++_molecules_read;

  Preprocess(m);

  const int matoms = m.natoms();

  std::unique_ptr<atype_t[]> atom_constant = std::make_unique<atype_t[]>(matoms);

  if (!_atom_typing.assign_atom_types(m, atom_constant.get(), nullptr)) {
    cerr << "Cannot assign atom types '" << m.name() << "'\n";
    return 0;
  }

  for (int i = 0; i < _nsfc; ++i) {
    _sfc[i].clear();
  }

  const iwecfp::FingerprintResult result =
      generator.Fingerprint(m, atom_constant.get(), _sfc.get());

  switch (result) {
    case iwecfp::FingerprintResult::kOk:
      return DoOutput(m, _sfc.get(), output);
    case iwecfp::FingerprintResult::kNoStartAtoms:
      return HandleNoStartAtoms(m, output);
    case iwecfp::FingerprintResult::kFatal:
      return 0;
  }

  return DoOutput(m, _sfc.get(), output);
}

int
Process(Options& options,
        iwecfp::Iwecfp& generator,
        Molecule& m,
        IWString_and_File_Descriptor& output) {
  return options.Fingerprint(m, generator, output);
}

int
Process(Options& options,
        iwecfp::Iwecfp& generator,
        data_source_and_type<Molecule>& input,
        IWString_and_File_Descriptor& output) {
  Molecule* m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (!Process(options, generator, *m, output)) {
      return 0;
    }
  }

  return 1;
}

int
ProcessFilterRecord(Options& options,
                    iwecfp::Iwecfp& generator,
                    const const_IWSubstring& buffer,
                    IWString_and_File_Descriptor& output) {
  assert(buffer.ends_with('>'));

  const_IWSubstring smiles(buffer);
  smiles.remove_up_to_first('<');
  smiles.chop();

  Molecule m;
  if (!m.build_from_smiles(smiles)) {
    cerr << "Cannot parse smiles '" << smiles << "'\n";
    return 0;
  }

  return Process(options, generator, m, output);
}

int
ProcessAsFilter(Options& options,
                iwecfp::Iwecfp& generator,
                iwstring_data_source& input,
                IWString_and_File_Descriptor& output) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    output << buffer << '\n';

    if (!buffer.starts_with(options.smiles_tag())) {
      continue;
    }

    if (!ProcessFilterRecord(options, generator, buffer, output)) {
      cerr << "Fatal error on line " << input.lines_read() << '\n';
      cerr << buffer << '\n';
      return 0;
    }

    options.MaybeFlush(output);
  }

  return output.good();
}

int
Process(Options& options,
        iwecfp::Iwecfp& generator,
        const char* fname, FileType input_type,
        IWString_and_File_Descriptor& output) {
  if (options.function_as_tdt_filter()) {
    iwstring_data_source input(fname);
    if (!input.good()) {
      cerr << "Cannot open '" << fname << "'\n";
      return 0;
    }

    return ProcessAsFilter(options, generator, input, output);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return Process(options, generator, input, output);
}


int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:g:J:i:fr:R:P:mlsM:bB:Q:Y:q:X:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');

  if (!process_elements(cl)) {
    cerr << "Cannot initialise elements (-E)\n";
    return 1;
  }

  if (!process_standard_aromaticity_options(cl, verbose)) {
    cerr << "Cannot process aromaticity options (-A)\n";
    return 1;
  }

  iwecfp::Iwecfp generator;
  if (!generator.Initialise(cl)) {
    std::cerr << "Cannot initialise iwecfp\n";
    return 1;
  }

  Options options;
  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise output handling\n";
    return 1;
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('f')) {
    // TDT filter, ignore input type.
  } else if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      Usage(6);
    }
  } else if (cl.size() == 1 && ::strcmp(cl[0], "-") == 0) {
    input_type = FILE_TYPE_SMI;
  } else if (!all_files_recognised_by_suffix(cl)) {
    cerr << "Cannot determine input type(s)\n";
    return 7;
  }

  IWString_and_File_Descriptor output(1);

  if (options.descriptor_file_output()) {
    options.WriteDescriptorFileHeader(output);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(2);
  }

  int rc = 0;
  for (const char* fname : cl) {
    if (!Process(options, generator, fname, input_type, output)) {
      rc = 1;
      break;
    }
  }

  output.flush();

  if (generator.verbose()) {
    generator.Report(std::cerr);
  }

  generator.Finalise();

  return rc;
}

}  // namespace iwecfp

int
main(int argc, char** argv) {
  int rc = iwecfp::Main(argc, argv);

  return rc;
}
