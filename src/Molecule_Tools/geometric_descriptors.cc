// Generate descriptors based on 3d distances.
// Could perhaps be made part of dbf, but this seemed cleaner.

#include <iostream>
#include <memory>
#include <optional>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

namespace geometric_descriptors {

using std::cerr;
using std::endl;

const char * prog_name = nullptr;

char output_separator = ' ';

int molecules_read = 0;

int molecules_not_matching_queries = 0;

struct JobParameters {
  resizable_array_p<Substructure_Query> queries;

  bool take_first_of_multiple_matches = false;

  bool skip_molecules_not_matching = false;

  // The number of atoms in each embedding
  int atoms_in_embedding = 0;

  // Given atoms_in_embedding atoms, how many pair-wise distances does that
  // genrate : n*(n-1)/2. Could be recomputed each time...
  int number_features = 0;

  IWString feature_name_stem = "d3d";
};

// Number of molecules with insufficient dimensionality.
int no_coordinates = 0;

void
Usage(int rc) {
  cerr << "Generates 2D distance features based on query matched atoms\n";
  cerr << " -q <queries>   query specification(s)\n";
  cerr << " -s <smarts>    query specification(s) as smarts\n";
  cerr << " -z ...         what to do with query matching problems\n";
  cerr << " -z i           ignore molecules not matching any query\n";
  cerr << " -z f           if there are multiple query matches, take the first one\n";
  cerr << " -v           verbose output\n";
  exit(rc);
}

void
Preprocess(Molecule& m) {
  return;
}

int
WriteHeader(const JobParameters& params,
            IWString_and_File_Descriptor& output) {
  output << "ID";
  for (int i = 0; i < params.atoms_in_embedding; ++i) {
    for (int j = i + 1; j < params.atoms_in_embedding; ++j) {
      output << output_separator << params.feature_name_stem << '.' << i << '.' << j;
    }
  }
  output << '\n';
  output.write_if_buffer_holds_more_than(2048);

  return 1;
}

std::optional<Set_of_Atoms>
DiscernIncludeAtoms(Molecule & m,
                    JobParameters& params,
                    IWString_and_File_Descriptor& output) {
  Molecule_to_Match target(&m);
  int max_atoms_matched = 0;
  for (Substructure_Query * q : params.queries) {
    Substructure_Results sresults;
    const int nhits = q->substructure_search(target, sresults);
    if (nhits == 0) {
      if (sresults.max_query_atoms_matched_in_search() > max_atoms_matched)
        max_atoms_matched = sresults.max_query_atoms_matched_in_search();
      continue;
    }
    if (nhits == 1) {
    } else if (params.take_first_of_multiple_matches) {
    } else {
      cerr << "Got " << nhits << " hits to query in " << m.name() << "\n";
      return std::nullopt;
    }
    const Set_of_Atoms * e = sresults.embedding(0);
    const int n = e->number_elements();
    // If we are the first query match, set the number of atoms expected in each embedding.
    if (params.atoms_in_embedding == 0) {
      params.atoms_in_embedding = n;
      params.number_features = n * (n - 1) / 2;
      WriteHeader(params, output);
    } else if (n != params.atoms_in_embedding) {
      cerr << "Inconsistent atom count in embeddings, got " << n << " expected " << params.atoms_in_embedding << '\n';
      return std::nullopt;
    }

    Set_of_Atoms result(*e);
    return result;
  }

  cerr << "No hits to " << m.name() << ", matched as many as " << max_atoms_matched << " query atoms\n";
  return std::nullopt;
}

int
GeometricDescriptors(Molecule & m,
                     const Set_of_Atoms& matched_atoms,
                     JobParameters& params,
                     IWString_and_File_Descriptor& output) {
  std::unique_ptr<float[]> result(new float[params.number_features]);

  const int natoms = matched_atoms.number_elements();

  int ndx = 0;
  for (int i = 0; i < natoms; ++i) {
    const atom_number_t a1 = matched_atoms[i];
    for (int j = i + 1; j < natoms; ++j, ++ndx) {
      const atom_number_t a2 = matched_atoms[j];
      float d = m.distance_between_atoms(a1, a2);
      result[ndx] = d;
    }
  }

  append_first_token_of_name(m.name(), output);
  for (int i = 0; i <  params.number_features; ++i) {
    output << output_separator << result[i];
  }
  output << '\n';

  output.write_if_buffer_holds_more_than(2048);

  return output.good();
}

// comment
int
GeometricDescriptors(Molecule & m,
                     JobParameters& params,
                     IWString_and_File_Descriptor& output) {
  if (m.highest_coordinate_dimensionality() == 1) {   // huh
    cerr << "No coodinates " << m.name() << " ignored\n";
    no_coordinates++;
    return 1;
  }

  const int matoms = m.natoms();
  if (matoms < 2) {
    cerr << "Too few atoms " << m.name() << " ignored\n";
    return 1;
  }

  std::optional<Set_of_Atoms> matched_atoms = DiscernIncludeAtoms(m, params, output);
  if (matched_atoms == std::nullopt) {
    molecules_not_matching_queries++;
    if (params.skip_molecules_not_matching) {
      return 1;
    }
    return 0;
  }

  return GeometricDescriptors(m, matched_atoms.value(), params, output);
}

int
GeometricDescriptors(data_source_and_type<Molecule> & input,
                     JobParameters& params,
                     IWString_and_File_Descriptor& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);
    molecules_read++;
    Preprocess(*m);
    if (! GeometricDescriptors(*m, params, output)) {
      cerr << "Error processing " << m->name() << '\n';
      return 0;
    }
  }

  return 1;
}

int
GeometricDescriptors(const char * fname,
                     FileType input_type,
                     JobParameters& params,
                     IWString_and_File_Descriptor& output) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
    assert(input_type != FILE_TYPE_INVALID);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good())
  {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }
  
  return GeometricDescriptors(input, params, output);
}
                     

int
GeometricDescriptors(int argc, char ** argv) {
  Command_Line cl(argc, argv, "A:E:vq:s:z:i:");
  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');

  if (! process_standard_aromaticity_options(cl, verbose))
    Usage(6);

  JobParameters params;

  if (cl.option_present('q')) {
    if (! process_queries(cl, params.queries, 'q', verbose)) {
      cerr << "Cannot read queries (-q)\n";
      return 1;
    }
  } else if (cl.option_present('s')) {
    IWString smarts;
    for (int i = 0; cl.value('s', smarts, i); ++i) {
      std::unique_ptr<Substructure_Query> q(new Substructure_Query);
      if (! q->create_from_smarts(smarts)) {
        cerr << "Cannot parse smarts '" << smarts << "'\n";
        return 1;
      }
      params.queries.add(q.release());
    }
  } else {
    cerr << "Must specify queries via the -q or -s options\n";
    Usage(1);
  }

  if (cl.option_present('z')) {
    IWString z;
    for (int i = 0; cl.value('z', z, i); ++i) {
      z.to_lowercase();
      if (z == 'i' || z.starts_with("ign")) {
        params.skip_molecules_not_matching = true;
        if (verbose)
          cerr << "Will skip molecules not matching queries\n";
      } else if (z == 'f' || z == "first") {
        params.take_first_of_multiple_matches = true;
        if (verbose)
          cerr << "Will take the first of any multiple matches\n";
      } else {
        cerr << "Unrecognised -z qualifier '" << z << "'\n";
        Usage(1);
      }
    }
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (! cl.option_present('i'))
  {
    if (1 == cl.number_elements() && 0 == strcmp("-", cl[0]))   // reading a pipe, assume smiles
      input_type = FILE_TYPE_SMI;
    else if (! all_files_recognised_by_suffix(cl))
    {
      cerr << "Cannot discern all file types, use the -i option\n";
      return 4;
    }
  }
  else if (! process_input_type(cl, input_type)) {
    Usage(3);
  }

  IWString_and_File_Descriptor output(1);

  for (const char * fname : cl) {
    if (! GeometricDescriptors(fname, input_type, params, output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }
  output.flush();

  if (verbose) {
    cerr << "processed " << molecules_read << " molecules\n";
    if (molecules_not_matching_queries)
      cerr << molecules_not_matching_queries << " molecules did not match the queries\n";
  }
  return 0;
}

}  // geometric_descriptors


int
main (int argc, char ** argv)
{
  geometric_descriptors::prog_name = argv[0];

  int rc = geometric_descriptors::GeometricDescriptors(argc, argv);

  return rc;
}
