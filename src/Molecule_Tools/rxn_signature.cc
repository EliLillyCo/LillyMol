/*
  main purpose is to generate signature information around the raction core, defined by the
  changing atoms
*/

#include <memory>

#include "cmdline.h"
#include "misc.h"
#include "accumulator.h"
#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION
#include "iwqsort.h"

#include "molecule_to_query.h"
#include "rxn_file.h"
#include "aromatic.h"
#include "atom_typing.h"
#include "iwstandard.h"
#include <vector>


const char * prog_name = NULL;

static int verbose = 0;

static int reactions_read = 0;

static Chemical_Standardisation chemical_standardisation;

static extending_resizable_array<int> acc_changing_atoms;

static int nr = 0;
static int * radius = nullptr;

static Atom_Typing_Specification atom_typing_specification;

static int reactions_discarded_for_incomplete_atom_map = 0;

static int remove_duplicate_reagents_atom_maps_scrambled = 0;

static extending_resizable_array<int> unmapped_atom_count;

static Accumulator<double> acc_unmapped_atom_fraction;

static int reactions_discarded_for_isotopic_atoms = 0;

static int reactions_with_no_reagent_atoms_in_products = 0;

static int reactions_containing_duplicate_atom_map_numbers = 0;

static IWString_and_File_Descriptor stream_for_changing_atoms;

static int ignore_bad_reactions = 0;

static int bad_reactions_ignored = 0;

static IWString_and_File_Descriptor stream_for_bad_reactions;

static int outputLeftSignature = true;

static int outputRightSignature = false;

static int outputAllParts = false;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Produces signatures around the atoms that change in a reaction\n";
  cerr << "  -r <rad>      radius from changing atoms to signature\n";
  cerr << "  -C <fname>    write changed atom counts to <fname>\n";
  cerr << "  -F <fname>    ignore otherwise bad reactions and write them to <fname>\n";
  cerr << "  -P <type>     atom typing specification\n";
  cerr << "  -O <type>     one of 'l', 'r', 'b' to output the signatures for the left-side parts of the reaction, the right-side parts, or both, respectively.  (Default 'l')\n";
  cerr << "  -a            all parts are reported in the signature (otherwise only the first reactant or product is reported \n";
  cerr << "  -v            verbose output\n";

  exit(rc);
}

typedef unsigned int atype_t;

static int
identify_atoms_in_signature(Molecule & m,
                              const int rad,
                              const int * changed,
                              int * in_signature)
{
  const int matoms = m.natoms();

  std::copy_n(changed, matoms, in_signature);    // by default, all changing atoms in fp

  if (rad > 0)
  {
    for (int i = 0; i < matoms; ++i)
    {
      if (0 == changed[i])    // we look from each of the changing atoms
        continue;

      const int ifrag = m.fragment_membership(i);

      for (int j = 0; j < matoms; ++j)   // check all the others
      {
        if (in_signature[j])
          continue;

        if (ifrag != m.fragment_membership(j))
          continue;

        if (m.bonds_between(i, j) <= rad)
          in_signature[j] = 1;
      }
    }
  }


  m.compute_aromaticity_if_needed();

  Set_of_Atoms extras;

  for (int i = 0; i < matoms; ++i)
  {
    if (0 == in_signature[i])
      continue;

    const Atom * a = m.atomi(i);

    const int acon = a->ncon();

    if (acon == m.nbonds(i))    // fully saturated, no double bonds here
      continue;

    for (int j = 0; j < acon; ++j)
    {
      const Bond * b = a->item(j);

      if (b->is_single_bond() || b->is_aromatic())
        continue;

      const int k = b->other(i);

      if (in_signature[k])
        continue;

      extras.add(k);
    }
  }

//  cerr << "Identified " << extras.size() << " doubly bonded extra in " << m.name() << ", radius " << rad << endl;
  extras.set_vector(in_signature, 1);
  

  return 1;
}

typedef class BondAtoms
{
    public:
    BondAtoms(int first,int second)
    {
        if (first > second)
        {
            atoms[0] = second;
            atoms[1] = first;              
        }
        else
        {
            atoms[0] = first;              
            atoms[1] = second;
        }
    }
    ~BondAtoms()
    {           
    }
    
    bool operator <(const class BondAtoms &other)
    {
        if(this->atoms[0] < other.atoms[0] || (this->atoms[0] == other.atoms[0] &&  this->atoms[1] < other.atoms[1]))
            return true;
        
        return false;
    }
    
    void swap(class BondAtoms &other)
    {
        int temp = this->atoms[0];
        this->atoms[0] = other.atoms[0];
        other.atoms[0] = temp;
        temp = this->atoms[1];
        this->atoms[1] = other.atoms[1];
        other.atoms[1] = temp;
    }
    
    int atoms[2];
} BondAtoms;
  
  
  
typedef class SignatureListItem
{
  private:
    std::vector<IWString> reactantList;
    std::vector<IWString> productList;
      
  
  public:
    SignatureListItem()
    {    
    }
    ~SignatureListItem()
    {           
    }
    
    void addReactantSignature(IWString sigToAdd)
    {
      reactantList.push_back(sigToAdd);
    }
    
    void addProductSignature(IWString sigToAdd)
    {
      productList.push_back(sigToAdd);
    }
    
    IWString toString()
    {
      // first sort the lists
      
      IWString returnVal("");
      
      if (!reactantList.empty())
      {
        std::sort(reactantList.begin(), reactantList.end());
        IWString thisSep("");
        for (std::vector<IWString>::const_iterator thisIter = reactantList.begin(); thisIter != reactantList.end() ; ++thisIter)
        {
          returnVal += thisSep + *thisIter;
          thisSep = "."; 
        }
      }
      if (!productList.empty())
      {
        if (!reactantList.empty())
          returnVal += ">>";  // separator if both reactants and products are to be included
          

        std::sort(productList.begin(), productList.end()); 
        IWString thisSep("");
        for (std::vector<IWString>::const_iterator thisIter = productList.begin(); thisIter != productList.end() ; ++thisIter)
        {
          returnVal += thisSep + *thisIter;
          thisSep = "."; 
        }
      }

      return returnVal;   
          
    }
} SignatureListItem;
  
  
  
static IWString
get_unique_smiles_core(ISIS_RXN_FILE_Molecule & r,
                         int radius,
                         atype_t * atype,
                         const int *changed)
{

  const int matoms = r.natoms();
  int * xref = new int[matoms]; std::unique_ptr<int[]> free_xref(xref);

  Molecule subset;
  
 
  //need to mark the changed atoms
  
  int *inFrag = new int[matoms+1]; std::unique_ptr<int[]> free_inFrag(inFrag);
  
  identify_atoms_in_signature(r, radius, changed, inFrag);
   
  for (int i = 0; i < matoms; ++i)
  {
    if (0 == changed[i])
      continue;

    atype[i] += 1702;
  }
  
  r.create_subset(subset, inFrag, 1, xref);

  const int nb = r.nedges();
 
  r.compute_aromaticity_if_needed();
  std::vector<BondAtoms> ringBonds;

  for (int i = 0; i < nb; ++i)
  {
    const Bond * b = r.bondi(i);

    const auto a1 = b->a1();
    const auto a2 = b->a2();

    const auto x1 = xref[a1];
    const auto x2 = xref[a2];

    //cerr << "Atom " << a1 << " mapped to " << x1 << " a2 " << a2 << " to " << x2 << endl;

    if (x1 >= 0)    // in the subset
      subset.set_isotope(x1, atype[a1]);

    if (x2 >= 0)
      subset.set_isotope(x2, atype[a2]);

    if (x1 < 0 || x2 < 0)    // not in the subset
      continue;

    if (b->nrings() > 0)
    {
        //printf("\nBond between atoms %i and %i is in a ring\n", x1,x2);
        Bond * subsetBond = const_cast<Bond *>(subset.bond_between_atoms(x1, x2));
        ringBonds.push_back(BondAtoms(x1,x2));
        
        if (b->is_aromatic())
        {
          //Bond * b = const_cast<Bond *>(subset.bond_between_atoms(x1, x2));
          subsetBond->set_bond_type(AROMATIC_BOND);
          subsetBond->set_permanent_aromatic(1);
        }
    }
  }

// now clear all implicit Hs from the core structure.
  
  int newAtomCount = subset.natoms();
  for (int i=0;  i < newAtomCount ; ++i)
  {
    subset.set_implicit_hydrogens(i, 0);
    
    //cerr << "i=" << i << "   ImplicitHs=" << subset.implicit_hydrogens(i) << "\n";
  }
  
  subset.setDoNotComputeAromaticity(true);  
  const IWString & s = subset.unique_smiles();
  subset.setDoNotComputeAromaticity(false);  
  
  int *smilesAtomOrder = new int[newAtomCount]; std::unique_ptr<int[]> free_imilesAtomOrder(smilesAtomOrder);
  
  subset.smiles_atom_order (smilesAtomOrder);

  // add the information about any ring bonds to the fragment "smiles".  This looks like this
  //  [4982C][11741O][6684C](:[4979C]):[4979C]::3,4:3,5::
  // indicating the the bonds between smiles atoms 3 and 4 and between 3 and 5 are ring bonds

  // first, replace the atoms indices that are relative to the subset with the indices related to position in the smiles
  
  std::ostringstream ss;
  ss <<  s;
  if (ringBonds.size())
  {
    for(std::vector<BondAtoms>::iterator thisOne = ringBonds.begin() ; thisOne != ringBonds.end() ; ++thisOne)  
      *thisOne = BondAtoms(smilesAtomOrder[thisOne->atoms[0]], smilesAtomOrder[thisOne->atoms[1]]);
      
    ss << "::";
 
    // sort them so that the result is still cannonical
 
    for(std::vector<BondAtoms>::iterator firstOne = ringBonds.begin() ; firstOne != ringBonds.end() ; ++firstOne)   
      for(std::vector<BondAtoms>::iterator secondOne = firstOne + 1 ; secondOne != ringBonds.end() ; ++secondOne)
        if(*secondOne < *firstOne)
            secondOne->swap(*firstOne);
    
    
    for(std::vector<BondAtoms>::iterator firstOne = ringBonds.begin() ; firstOne != ringBonds.end() ; ++firstOne)
        ss << firstOne->atoms[0] <<","  << firstOne->atoms[1] << ":";
  
  }

  return IWString(ss.str());

}




static int
write_changing_atom_count(RXN_File & rxn,
                          const int c,
                          IWString_and_File_Descriptor & output)
{
  IWString s;
  s = rxn.name();

  if (s.ends_with(".rxn"))
    s.chop(4);

  output << s;
  output << ' ' << c << '\n';
   
  output.write_if_buffer_holds_more_than(4096);

  return 1;
}


static int
ensure_locator_arrays_filled(RXN_File & rxn, const int highestMapNumber)
{
  Molecule_to_Query_Specifications mqs;
  mqs.set_make_embedding(1);
  mqs.set_interpret_atom_alias_as_smarts(0);

  IWReaction notused;
  RXN_File_Create_Reaction_Options rxcfcro;
  rxcfcro.set_only_create_query_from_first_reagent(1);

  int *tmp = new_int(highestMapNumber+1, 1); std::unique_ptr<int[]> free_tmp(tmp);
  if (! rxn.create_reaction(notused, rxcfcro,  mqs, tmp))
  {
    cerr << "Cannot create reaction\n";
    return 0;
  }

  return 1;
}

static int
leftSideSignature(RXN_File & rxn,std::vector<SignatureListItem> &signatures, int &countOfChangedAtoms)
{
  int molCount;

  molCount = rxn.number_reagents();

  for (int reagentIndex = 0 ; reagentIndex < molCount; ++reagentIndex)
  {
    if (!outputAllParts && reagentIndex > 0)
      break;

    ISIS_RXN_FILE_Molecule *thisMol(NULL);
   
    thisMol = &rxn.reagent(reagentIndex);
      
    chemical_standardisation.process(*thisMol);
      
    const int matoms = thisMol->natoms();

    int xx = rxn.highest_atom_map_number();
    if (xx < matoms)
      xx = matoms;

    int * changed = new_int(matoms); std::unique_ptr<int[]> free_changed(changed);
    
  // We need to create a reaction in order for the reagent and product locator arrays to be filled

    ensure_locator_arrays_filled(rxn, xx);

    Changing_Atom_Conditions cac;
    cac.set_is_changing_if_different_neighbours(1);
    cac.set_ignore_lost_atom_if_isolated(1);
    cac.set_include_changing_bonds_in_changing_atom_count(1);
    cac.set_consider_aromatic_bonds(0);


    countOfChangedAtoms +=  rxn.identify_atoms_changing_reagent(reagentIndex, atom_typing_specification, changed, cac);   

    
    if (! rxn.at_least_some_mapped_atoms_common_btw_reagents_and_products())
    {
      reactions_with_no_reagent_atoms_in_products++;
      return 0;
    }

    if (rxn.contains_duplicate_atom_map_numbers())
    {
      reactions_containing_duplicate_atom_map_numbers++;
      return 0;
    }

    if (thisMol->number_fragments() > 1)
    {
      cerr << "Caution, first reagent has multiple fragments\n";
    }
   
      
    atype_t *atype = new atype_t[matoms]; std::unique_ptr<atype_t[]> free_atype(atype);

    if (! atom_typing_specification.assign_atom_types(*thisMol, atype))
    {
      cerr << "Cannot assign atom types '" << rxn.name() << "'\n";
      return 0;
    }

  // Now we can do what we came here for

    for (int i = 0; i < matoms; ++i)
    {
      if (changed[i])
      {
        thisMol->set_isotope(i, 1);
        changed[i] = 1;  // the atoms that changed because of aromatic changes too. These will be "2" in the changed array, we need them to be 1
      }
    }
    
    for (int i = 0; i < nr; ++i)
    {
      IWString thisSig = get_unique_smiles_core(*thisMol, radius[i], atype, changed);
      signatures[i].addReactantSignature(thisSig); 
    }
  }
  
  return 1;
}

static int
rightSideSignature(RXN_File & rxn,std::vector<SignatureListItem>  &signatures, int &countOfChangedAtoms)
{
  countOfChangedAtoms = 0;
  
  int molCount;

  molCount = rxn.number_products();
    
  for (int productIndex = 0 ; productIndex < molCount; ++productIndex)
  {
    if (!outputAllParts && productIndex > 0)
      break;
    
    ISIS_RXN_FILE_Molecule *thisMol(NULL);
 
    thisMol = &rxn.product(productIndex);

    chemical_standardisation.process(*thisMol);
      
    const int matoms = thisMol->natoms();

    int xx = rxn.highest_atom_map_number();
    if (xx < matoms)
      xx = matoms;

    int * changed = new_int(matoms); std::unique_ptr<int[]> free_changed(changed);
    
  // We need to create a reaction in order for the reagent and product locator arrays to be filled

    ensure_locator_arrays_filled(rxn, xx);

    Changing_Atom_Conditions cac;
    cac.set_is_changing_if_different_neighbours(1);
    cac.set_ignore_lost_atom_if_isolated(1);
    cac.set_include_changing_bonds_in_changing_atom_count(1);
    cac.set_consider_aromatic_bonds(0);


    countOfChangedAtoms +=  rxn.identify_atoms_changing_product(productIndex, atom_typing_specification, changed, cac);
 
    if (! rxn.at_least_some_mapped_atoms_common_btw_reagents_and_products())
    {
      reactions_with_no_reagent_atoms_in_products++;
      return 1;
    }

    if (rxn.contains_duplicate_atom_map_numbers())
    {
      reactions_containing_duplicate_atom_map_numbers++;
      return 1;
    }

    if (thisMol->number_fragments() > 1)
    {
      cerr << "Caution, first product has multiple fragments\n";
    }
   
   
    atype_t *atype = new atype_t[matoms]; std::unique_ptr<atype_t[]> free_atype(atype);

    if (! atom_typing_specification.assign_atom_types(*thisMol, atype))
    {
      cerr << "Cannot assign atom types '" << rxn.name() << "'\n";
      return 0;
    }

  // Now we can do what we came here for

    for (int i = 0; i < matoms; ++i)
    {
      if (changed[i])
      {
        thisMol->set_isotope(i, 1);
        changed[i] = 1;  // the atoms that changed because of aromatic changes too. These will be "2" in the changed array, we need them to be 1
      }
    }
    
    for (int i = 0; i < nr; ++i)
    {
      IWString thisSig = get_unique_smiles_core(*thisMol, radius[i], atype, changed);
      signatures[i].addProductSignature(thisSig); 
    }
  }
  
  return 1;
}


static int
rxn_signature(RXN_File & rxn,
                IWString_and_File_Descriptor & output)
{
  const int initial_nr = rxn.number_reagents();

  if (0 == initial_nr)
  {
    cerr << "Skipping reaction with no reagents " << rxn.name() << endl;
 
    return 1;
  }
  
  std::vector<SignatureListItem> signatures;
  
  for (int i = 0; i < nr; ++i)
  {
    signatures.push_back(SignatureListItem());
  }
  int countOfChangedAtoms = 0 , countOfChangedProductAtoms = 0;
  
  if (outputLeftSignature)
    if (!leftSideSignature(rxn,signatures, countOfChangedAtoms))
      return 0;
 
 
  if (outputRightSignature)
  {
    if (!rightSideSignature(rxn,signatures, countOfChangedProductAtoms))
      return 0;
    if (!outputLeftSignature  || countOfChangedProductAtoms > countOfChangedAtoms)
      countOfChangedAtoms = countOfChangedProductAtoms;
  }

     
  if (stream_for_changing_atoms.is_open())
    write_changing_atom_count(rxn, countOfChangedAtoms, stream_for_changing_atoms);
      
  acc_changing_atoms[countOfChangedAtoms]++;
  
  Reaction_Smiles_Options opts;
  opts.set_write_reaction_name(0);
  opts.set_write_agent(1);
  
  IWString s;
  s = rxn.name();

  if (s.ends_with(".rxn"))
    s.chop(4);
  
  rxn.write_rxn_smiles(opts, output);      
  output << " " << countOfChangedAtoms<< " " << s;

  output.write_if_buffer_holds_more_than(4096);
  
  for (int i = 0; i < nr; ++i)
  {
    output << ' ' << signatures[i].toString();
  }


  output << '\n';
 
  output.write_if_buffer_holds_more_than(4096);

  return output.good();
}

static int
echo_bad_data(iwstring_data_source & input,
              const off_t initial_offset,
              IWString_and_File_Descriptor & output)
{
  const auto current_offset = input.tellg();
  input.seekg(initial_offset);
//cerr << "Bad reaction begin " << initial_offset << " now " << current_offset << endl;
  input.echo(output, (current_offset - initial_offset));
  input.seekg(current_offset, 0);
  assert (input.tellg() == current_offset);
  output.write_if_buffer_holds_more_than(4096);

  return input.tellg() == current_offset;
}

static int
next_reaction_smiles(iwstring_data_source & input,
                     RXN_File & rxn)
{
  const_IWSubstring buffer;

  if (! input.next_record(buffer))
    return 0;

//cerr << "SMILES INPUT '" << buffer << "'\n";

  if (! rxn.build_from_reaction_smiles(buffer))
  {
    cerr << "Cannot interpret " << buffer << endl;
    return 0;
  }

  return 1;
}


static int rxn_signature(const char * fname, IWString_and_File_Descriptor & output);


static int
rxn_signature_file_containing_reaction_files(const char * fname,
                                               IWString_and_File_Descriptor & output)
{
  IWString tmp(fname);
  assert (tmp.starts_with("F:"));
  tmp.remove_leading_chars(2);

  iwstring_data_source input(tmp.null_terminated_chars());

  if (! input.good())
  {
    cerr << "Cannot open file of reactions '" << tmp << "'\n";
    return 0;
  }

  IWString buffer;

  while (input.next_record(buffer))
  {
    IWString thisFile(buffer);
    
    if (thisFile.length() != 0 && ! rxn_signature(thisFile.chars(), output))
      return 0;
  }

  return 1;
}

static int
rxn_signature (const char * fname,
                 IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);
  if (! input.good())
  {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  while (1)
  {
    RXN_File rxn;
    rxn.set_do_automatic_atom_mapping(0);

    const auto initial_offset = input.tellg();

    if (! next_reaction_smiles(input, rxn))
    {
      if (input.eof())
        return 1;

      return 0;
    }

    if (rxn.contains_isotopic_reagent_atoms())
    {
      if (verbose)
        cerr << rxn.name() << " contains isotopic atoms\n";
      reactions_discarded_for_isotopic_atoms++;
      continue;
    }

    reactions_read++;

    if (verbose > 1)
      cerr << "Processing " << rxn.name() << "\n";

    if (! rxn_signature(rxn, output))
    {
      cerr << "rxn_signature:fatal error processing '" << rxn.name() << "' now at line " << input.lines_read() << endl;

      if (! ignore_bad_reactions)
        return 0;

      if (stream_for_bad_reactions.is_open())
        echo_bad_data(input, initial_offset, stream_for_bad_reactions);

      bad_reactions_ignored++;
    }
  }

  return 1;

}

static int
get_radii (const Command_Line & cl,
           const char flag,
           resizable_array<int> & radii)
{
  const_IWSubstring r;
  for (int i = 0; cl.value(flag, r, i); ++i)
  {
    int j = 0;
    const_IWSubstring token;
    while (r.nextword(token, j, ','))
    {
      int x;
      if (! token.numeric_value(x) || x < 0)
      {
        cerr << "Invalid radius '" << token << "'\n";
        return 0;
      }

//    cerr << "interpreted as " << x << endl;

      radii.add_if_not_already_present(x);
    }
  }

  radii.iwqsort_lambda([](const int r1, const int r2) {    // not really necessary to sort them
                           if (r1 < r2)
                             return -1;
                           else if (r1 > r2)
                             return 1;
                           else
                             return 0;         // will never happen
                            });

  return radii.number_elements();
}


static int rxn_signature (int argc, char ** argv) 
{ 
  Command_Line cl(argc, argv, "vr:C:F:P:O:a");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('a'))
    outputAllParts= true;
    

  set_global_aromaticity_type(Daylight);

  if (cl.option_present('P'))
  {
    const_IWSubstring p = cl.string_value('P');

    if (! atom_typing_specification.build(p))
    {
      cerr << "Invalid atom typing specification '" << p << "'\n";
      return 1;
    }
  }
  else    // too many difficulties if the atom typing ends up unset
  {
    atom_typing_specification.build("UST:AZUCORS");
  }


//  if (! atom_typing_specification.build("UST:AZUCORS"))
//  {
//    cerr << "INvalid atom typing specification 'UST:AZUCORS'\n";
//    return 1;
//  }

  set_iwreaction_display_no_atoms_in_query_message(0);
  
  resizable_array<int> tmpradii;
  nr = get_radii(cl, 'r', tmpradii);

  if (0 == nr)
  {
    cerr << "Must specify radius from changing atoms via the -r option\n";
    usage(1);
  }

  radius = new_int(nr);std::unique_ptr<int[]> free_radius(radius);

  for (int i = 0; i < nr; ++i)
  {
    radius[i] = tmpradii[i];
  }
  

  if (verbose)
  {
    cerr << "radii";
    for (int i = 0; i < nr; ++i)
    {
      cerr << " " << radius[i];
    }
    cerr << endl;
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.option_present('C'))
  {
    const char * c = cl.option_value('C');

    if (! stream_for_changing_atoms.open(c))
    {
      cerr << "Cannot open stream for changing atom counts '" << c << "'\n";
      return 1;
    }

    if (verbose)
      cerr << "Counts of changing atoms written to '" << c << "'\n";

    stream_for_changing_atoms << "ID RGNT_Change\n";
  }
  
  if (cl.option_present('O'))
  {
    const char * val = cl.option_value('O');

    char thisChar = 'z';   // not an allowed option
    if (strlen(val) == 1)
      thisChar = val[0];
      
    switch (thisChar)
    {
        case 'l':
            outputLeftSignature = true;
            outputRightSignature = false;
            if (verbose)
                cerr << "Reactant Signatures will be generated\n";
            break;         
        case 'r':
            outputLeftSignature = false;
            outputRightSignature = true;
            if (verbose)
                cerr << "Product Signatures will be generated\n";
            break;         
        case 'b':
            outputLeftSignature = true;
            outputRightSignature = true;
            if (verbose)
                cerr << "Combined Reactant and Product Signatures will be generated\n";
           break;
                        
        default:
          cerr << "option -O requires 'r', 'p', or 'b' (reactants, products or both)\n";
          usage(1);

    }
  
  
  }
  
  
  if (cl.option_present('F'))
  {
    IWString f = cl.string_value('F');

    ignore_bad_reactions = 1;

    if ('.' == f)
    {
      if (verbose)
        cerr << "Will ignore failed reactions\n";
    }
    else
    {
      if (! stream_for_bad_reactions.open(f.null_terminated_chars()))
      {
        cerr << "Cannot open stream for bad reactions '" << f << "'\n";
        return 1;
      }

      if (verbose)
        cerr << "Will write failed reactions to '" << f << "'\n";
    }
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    const char *fname = cl[i];
    assert(NULL != fname);

    const_IWSubstring tmp(fname);
    if (tmp.starts_with("F:"))
    {      
      if (! rxn_signature_file_containing_reaction_files(fname, output))
      {
        rc = i + 1;
        break;
      } 
    }
    else
    {
      if (! rxn_signature(cl[i], output))
      {
        rc = i + 1;
        break;
      }
    }
  }

  output.flush();

  if (verbose)
  {
    cerr << "Read " << reactions_read << " reactions\n";
    Accumulator_Int<int> acc;

    int most_common = 0;
    int most_common_count = 0;
    for (int i = 0; i < acc_changing_atoms.number_elements(); ++i)
    {
      const int c = acc_changing_atoms[i];
      if (0 == c)
        continue;

      acc.extra(i, c);
      cerr << c << " reactions had " << i << " changing atoms\n";

      if (c > most_common_count)
      {
        most_common_count = c;
        most_common = i;
      }
    }

    if (acc.n() > 0)
    {
      cerr << "Changing atom counts btw " << acc.minval() << " and " << acc.maxval() << " ave " << static_cast<float>(acc.average()) << endl;
      cerr << "Purity " << static_cast<float>(most_common_count) / static_cast<float>(acc.n()) << endl;
    }
  }

  if (bad_reactions_ignored)
    cerr << "SKIPPED " << bad_reactions_ignored << " bad reactions\n";

  if (reactions_discarded_for_isotopic_atoms)
    cerr << "SKIPPED " << reactions_discarded_for_isotopic_atoms << " reactions containing isotopic atoms\n";

  if (reactions_with_no_reagent_atoms_in_products)
    cerr << "SKIPPED " << reactions_with_no_reagent_atoms_in_products << " reactions_with_no_reagent_atoms_in_products\n";
  
  if (reactions_containing_duplicate_atom_map_numbers)
    cerr << "SKIPPED " << reactions_containing_duplicate_atom_map_numbers << " reactions_containing_duplicate_atom_map_numbers\n";

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = rxn_signature(argc, argv);

  return 0;
}
