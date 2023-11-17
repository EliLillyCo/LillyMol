#include <stdlib.h>

#include "iwstring_data_source.h"
#include "iw_auto_array.h"

#include "molecule.h"
#include "rwmolecule.h"
#include "chiral_centre.h"

#include "inchi_api.h"

int inchi_to_inchi_key(const char * inchi,IWString & key)
{
  char ikey[100];
  int ret=GetStdINCHIKeyFromStdINCHI(inchi,ikey);
  if ( ret == INCHIKEY_OK)
  {
    key=ikey;
    return 1;
  }
  else
  {
    cerr << "invalid InChI '" << inchi << "'\n";
    return 0;
  }
}

int
Molecule::build_from_inchi(const const_IWSubstring & inchi)
{
  if (! inchi.starts_with("InChI=1"))
  {
    cerr << "Molecule::build_from_inchi:invalid InChI '" << inchi << "'\n";
    return 0;
  }

  const_IWSubstring myinchi(inchi);   // so we can split off the molecule name

  if (inchi.nwords() > 1)
  {
    _molecule_name = inchi;
    _molecule_name.remove_leading_words(1);
    myinchi.truncate_at_first(' ');
  }

  IWString tmp(myinchi);   // in case we need a null terminated version

// Now ask INCHI API to build an INCHI molecule from the INCHI string

  inchi_InputINCHI inpInchi;
  inchi_OutputStruct outStructure;
  memset( &inpInchi, 0, sizeof(inpInchi) );
  memset( &outStructure, 0, sizeof(outStructure) );
  
  char *szInChI = new char[tmp.nchars()+1]; iw_auto_array<char> free_szInChI(szInChI);
  tmp.copy_to_char_array(szInChI);

  inpInchi.szInChI = szInChI;
  inpInchi.szOptions =NULL;
  int ret = GetStructFromStdINCHI(&inpInchi,&outStructure);

  if ( ret != inchi_Ret_OKAY && ret != inchi_Ret_WARNING)
  { 
    cerr << "Molecule::build_from_inchi:cannot parse InChI '" << inchi << "'\n";
    if(outStructure.szMessage)
    {
      cerr << "Reason: '" << outStructure.szMessage << "'\n";
    }
    return 0;
  }
    
  // transfer all information from the INCHI molecule to THIS.
  for(int i=0;i<outStructure.num_atoms;++i)
  {
    inchi_Atom &at=outStructure.atom[i];
    
    const Element * e = get_element_from_symbol_no_case_conversion(at.elname);
    if (NULL == e)
    {
      cerr << "Molecule::build_from_inchi:cannot create element from '" << at.elname << "'\n";
      return 0;
    }
    Atom * a=new Atom(e);
    int iso=at.isotopic_mass;
    if(iso)
    {
      a->set_isotope(iso-ISOTOPIC_SHIFT_FLAG+ e->normal_isotope());
    }
    if(at.charge)
    {
      a->set_formal_charge(at.charge);
    }
    if(at.radical)
    {
      a->set_implicit_hydrogens(0,1);
      a->set_implicit_hydrogens_known(1);
    }
    //a->setxyz(at.x,at.y,at.z);
    add(a);
  }
  //add bonds after all atoms are added
  for(int i=0;i<outStructure.num_atoms;++i)
  {
    inchi_Atom &at=outStructure.atom[i];
    for(int j=0;j<at.num_bonds;++j)
    {
       //bond_stereo is ignored..
       if(at.bond_stereo[j] !=0)
         cerr<< "bond_stereo: " <<i <<" "<<at.neighbor[j]<<" type="<<(int)at.bond_type[j]<<" stereo="<<(int)at.bond_stereo[j]<< "\n";

      // bond are set twice in inchi_Atoms
      if(_things[i]->is_bonded_to(at.neighbor[j]))
	continue;

      bond_type_t btype=INVALID_BOND_TYPE;
      switch(at.bond_type[j])
      {
	case INCHI_BOND_TYPE_SINGLE: btype=SINGLE_BOND;
				     break;
	case INCHI_BOND_TYPE_DOUBLE: btype=DOUBLE_BOND;
				     break;
	case INCHI_BOND_TYPE_TRIPLE: btype=TRIPLE_BOND;
				     break;
	default:
				     break;
      }
      add_bond(i,at.neighbor[j],btype,1);
    }
  }
  //add H D T  after all atoms and bonds are added
  for(int i=0;i<outStructure.num_atoms;++i)
  {
    inchi_Atom &at=outStructure.atom[i];
    //cerr<<"Hs: "<<(int)at.num_iso_H[0]<<" "<<(int)at.num_iso_H[1]<<" "<<(int)at.num_iso_H[2]<<" "<<(int)at.num_iso_H[3]<<"\n";
    for(int j=0;j<4;++j)
    {
      for(int k=0;k<at.num_iso_H[j];++k)
      {
        const Element * e = get_element_from_symbol_no_case_conversion("H");
        Atom * a=new Atom(e);
	if(j)
          a->set_isotope(j);
	add(a);
	add_bond(i,_number_elements -1,SINGLE_BOND,1);
      }
    }
  }

  for(int i=0;i<outStructure.num_stereo0D;++i)
  {
    inchi_Stereo0D &st=outStructure.stereo0D[i];
    if(INCHI_StereoType_Tetrahedral == st.type )
    {
      //cerr<<"Chiral "<<(int)st.type<<" "<<(int)st.parity<<" "<<st.central_atom<<" "<<st.neighbor[0]<<" "<<st.neighbor[1]<<" "<<st.neighbor[2]<<" "<<st.neighbor[3]<<"\n";

      if( INCHI_PARITY_EVEN != st.parity && INCHI_PARITY_ODD != st.parity)
	continue;

      Chiral_Centre * c = new Chiral_Centre(st.central_atom);
       
      c->set_chirality_known(1);
      
      if(st.central_atom != st.neighbor[0])
        c->set_top_front(st.neighbor[0]);
      else if (_things[st.central_atom]->ncon() != 4)
	c->set_top_front(CHIRAL_CONNECTION_IS_LONE_PAIR);
      else // the last one is H (already make explicit)
	c->set_top_front(other(st.central_atom,3));

      c->set_top_back(st.neighbor[1]);

      if(INCHI_PARITY_EVEN == st.parity)
      {
        c->set_left_down(st.neighbor[2]);
        c->set_right_down(st.neighbor[3]);
      }
      else if(INCHI_PARITY_ODD == st.parity)
      {
        c->set_right_down(st.neighbor[2]);
	c->set_left_down(st.neighbor[3]);
      }
      _chiral_centres.add(c);
    }
    else if (INCHI_StereoType_DoubleBond == st.type)
    {
      //int begin = st.neighbor[1];
      //int end = st.neighbor[2];
      //int begin_neib = st.neighbor[0];
      //int end_neib = st.neighbor[4];
      cerr << "Molecule::build_from_inchi: ignore INCHI_StereoType_DoubleBond \n";
    }
    else
    {
      cerr << "Molecule::build_from_inchi: ignore unsupported stereo type \n";
    }
    
  }
   
  FreeStructFromStdINCHI(&outStructure);
  return 1;
}

int
Molecule::read_molecule_inchi_ds(iwstring_data_source & input)
{
  resize(0);   // clear out anything already present

  const_IWSubstring buffer;

  EXTRA_STRING_RECORD (input, buffer, "read mol inchi");

  return build_from_inchi(buffer);
}

//#define DEBUG_MAKE_INCHI_STRING

int
Molecule::InChI(IWString & inchi_string)      // maybe this could be a const method, depends on what inchi needs...
{
// compute_aromaticity_if_needed();    does inchi require aromaticity and/or ring perception? Method is non const

// Instantiate an INCHI molecule, transfer all information from THIS to the inchi molecule
// then get the inchi string from the inchi molecule and put into S

  inchi_Input inchi_inp;
  memset( &inchi_inp, 0, sizeof(inchi_inp) );
  inchi_Atom * at =new inchi_Atom [_number_elements]; unique_ptr<inchi_Atom> free_a(at);
  memset(at,0,_number_elements * sizeof(inchi_Atom));
  inchi_inp.atom=at;
  inchi_inp.num_atoms= _number_elements;

  for (int i = 0; i < _number_elements; i++)
  {
    Atom * a = _things[i];
    int ch=a->formal_charge();
    a->atomic_symbol().copy_to_char_array(at[i].elname);
    int iso=a->isotope();
    if(iso)
      iso = ISOTOPIC_SHIFT_FLAG + iso - a->element()->normal_isotope();
    at[i].isotopic_mass = iso;
    at[i].charge = ch;
    at[i].radical = 0;
    at[i].x = a->x();
    at[i].y = a->y();
    at[i].z = a->z();
    at[i].num_iso_H[0]= a->implicit_hydrogens(); 
    int icon = a->ncon();
    at[i].num_bonds = icon;
    for (int j = 0; j < icon; j++)
    {
      const Bond * b = a->item (j);

      S_CHAR bond_type=0,s1=0,s2=0;
      if(b->is_single_bond())
      {
        bond_type=INCHI_BOND_TYPE_SINGLE;
      }
      else if(b->is_double_bond())
      {
        bond_type = INCHI_BOND_TYPE_DOUBLE;
      }
      else if(b->is_triple_bond ())
      {
        bond_type = INCHI_BOND_TYPE_TRIPLE;
      }
      else
      {
        cerr << "Molecule::InChI:unrecogised bond type " << b->btype() << endl;
        bond_type = INCHI_BOND_TYPE_SINGLE;
      }
      at[i].bond_type[j]=bond_type;
      at[i].neighbor[j]=static_cast<AT_NUM>(b->other(i));

      if(b->is_wedge_up())
      {
        s1 = INCHI_BOND_STEREO_SINGLE_1UP;
	s2 = INCHI_BOND_STEREO_SINGLE_2UP;
      }
      else if (b->is_wedge_down())
      { 
        s1 = INCHI_BOND_STEREO_SINGLE_1DOWN;
	s2 = INCHI_BOND_STEREO_SINGLE_2DOWN;
      }
      else if ( b->is_wedge_either())
      {
        s1 = INCHI_BOND_STEREO_SINGLE_1EITHER;
	s2 = INCHI_BOND_STEREO_SINGLE_2EITHER;
      }
      else if ( b->is_cis_trans_either_double_bond())
      {
        s1 = INCHI_BOND_STEREO_DOUBLE_EITHER;
        s2 = INCHI_BOND_STEREO_DOUBLE_EITHER;
      }
      else
      {
        s1 = INCHI_BOND_STEREO_NONE;
        s2 = INCHI_BOND_STEREO_NONE;
      }

      //if(s1 || s2 )
	//cerr<< "Molecule::InChI bond stereo "<< i <<" "<< (int)s1<<" " <<(int)s2<<"\n";

      if(b->a1() == i)
	at[i].bond_stereo[j]=s1;
      else
        at[i].bond_stereo[j]=s2;
    }
  }

  int ncc = _chiral_centres.number_elements();
  inchi_Stereo0D  * st =new inchi_Stereo0D [ncc]; unique_ptr<inchi_Stereo0D> free_s(st);
  memset(st,0,ncc * sizeof(inchi_Stereo0D));
  inchi_inp.stereo0D=st;
  inchi_inp.num_stereo0D=ncc;
  for(int i=0;i<ncc;++i)
  {
    const Chiral_Centre *c =_chiral_centres[i];
    inchi_Stereo0D &s=st[i];
    s.type=INCHI_StereoType_Tetrahedral;
    if(c->chirality_known())
      s.parity=INCHI_PARITY_EVEN;
    else
      s.parity=INCHI_PARITY_UNKNOWN;
    s.central_atom=c->a();

#ifdef DEBUG_MAKE_INCHI_STRING
    cerr << "Central atom is " << c->a() << endl;
#endif

    int chiral_centre_iterator = 0;

    int nas=c->number_atoms_specified();

    //When implicit_H or lone pair, neighbor[0] is set to central_atom
    if(nas < 4) s.neighbor[0]= c->a();

    for (int j=0 ; j < nas; j++)
    {
      atom_number_t a = c->next_atom(chiral_centre_iterator);
#ifdef DEBUG_MAKE_INCHI_STRING
      cerr << " bonded to " << a << endl;
#endif
      s.neighbor[4-nas+j] = a;
    }
   
  }

  inchi_inp.szOptions = NULL;

  inchi_Output inchi_out;
  memset(&inchi_out, 0, sizeof(inchi_out) );

  int ret=GetStdINCHI(&inchi_inp,&inchi_out);
  if(ret == inchi_Ret_OKAY || ret == inchi_Ret_WARNING)
  {
    inchi_string=inchi_out.szInChI;

#ifdef DEBUG_MAKE_INCHI_STRING
    cerr << "Successful INCHI Production\n";
    cerr<<inchi_string<<endl;
    cerr<<inchi_out.szAuxInfo<<endl;
#endif
  }
  else
  {
    cerr << "Error from GetStdINCHI " << ret << ", processing '" << _molecule_name << "'\n";
    cerr << inchi_out.szMessage <<endl;
    inchi_string.resize(0);
    return 0;
  }

  FreeStdINCHI(&inchi_out);

  return 1;
}

int
Molecule::write_molecule_inchi (ostream & output)
{
  IWString tmp;

  InChI (tmp);   // generate string version of InChI

  output << tmp << ' ' << _molecule_name << '\n';

  return 1;
}
