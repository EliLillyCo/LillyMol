#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "Foundational/iwaray/iwaray.h"

#include "iwreaction.h"

template class resizable_array_p<Inter_Particle_Bond>;
template class resizable_array_base <Inter_Particle_Bond *>;

template class resizable_array_p<Sidechain_Reaction_Site>;
template class resizable_array_base<Sidechain_Reaction_Site *>;

template class resizable_array_p<Molecule_and_Embedding>;
template class resizable_array<Molecule_and_Embedding *>;
template class resizable_array_base<Molecule_and_Embedding *>;

template class resizable_array<const Atom *>;
template class resizable_array_base<const Atom *>;

template class resizable_array_p<Reaction_Change_Element>;
template class resizable_array_base<Reaction_Change_Element *>;

template class resizable_array_p<Reaction_Formal_Charge>;
template class resizable_array_base<Reaction_Formal_Charge *>;

template class resizable_array_p<Reaction_Change_Formal_Charge>;
template class resizable_array_base<Reaction_Change_Formal_Charge *>;

template class resizable_array_p<Reaction_Place_Isotope>;
template class resizable_array_base<Reaction_Place_Isotope *>;

template class resizable_array_p<Reaction_Increment_Isotope>;
template class resizable_array_base<Reaction_Increment_Isotope *>;

template class resizable_array_p<Reaction_Invert_Isotope>;
template class resizable_array_base<Reaction_Invert_Isotope *>;

template class resizable_array_p<IWReaction>;
template class resizable_array_base<IWReaction *>;

template class resizable_array_p<No_Reaction>;
template class resizable_array_base<No_Reaction *>;

template class resizable_array_p<Reaction_Stereo_Centre>;
template class resizable_array_base<Reaction_Stereo_Centre *>;

template class resizable_array_p<Reaction_Wedge_Bond>;
template class resizable_array_base<Reaction_Wedge_Bond *>;

template class resizable_array_p<Reaction_Dihedral_Angle>;
template class resizable_array_base<Reaction_Dihedral_Angle *>;

template class resizable_array_p<Reaction_Bond_Length>;
template class resizable_array_base<Reaction_Bond_Length *>;

template class resizable_array_p<Reaction_Bond_Angle>;
template class resizable_array_base<Reaction_Bond_Angle *>;

template class resizable_array_p<Reaction_3D_Replace>;
template class resizable_array_base<Reaction_3D_Replace *>;

template class resizable_array_p<Replace_Atom>;
template class resizable_array_base<Replace_Atom *>;
template class resizable_array<Replace_Atom *>;

#ifdef DONTNEEDTHISANYMORE

static void
unused ()
{
  resizable_array<const Atom *> ra;
  ra.resize (0);
  ra.add (nullptr);
  ra[0];

  resizable_array_p<Sidechain_Reaction_Site> srs;
  resizable_array_base<Sidechain_Reaction_Site *> srsb;

  srs.add (nullptr);
  srsb[0];

  resizable_array_p<Reaction_Formal_Charge> rfc;
  resizable_array_base<Reaction_Formal_Charge *> rfcb;

  rfc.add (nullptr);

  rfc[0];

  resizable_array_p<Reaction_Place_Isotope> rpi;
  resizable_array_base<Reaction_Place_Isotope *> rpib;

  rpi.add (nullptr);

  rpi[0];

  resizable_array_p<Reaction_Change_Element> rce;
  resizable_array_base<Reaction_Change_Element *> rceb;

  rce.add (nullptr);
  rce[0];

  resizable_array_p<Reaction_Stereo_Centre> rsc;
  resizable_array_base<Reaction_Stereo_Centre *> rscb;

  rsc.add (nullptr);
  rsc[0];
  rscb.last_item ();

  resizable_array_p<Inter_Particle_Bond> ipb;
  resizable_array_base<Inter_Particle_Bond *> ipbb;

  ipb.add (nullptr);
  ipb[0];

  resizable_array_p<Molecule_and_Embedding> mae;
  resizable_array<Molecule_and_Embedding *> maep;
  resizable_array_base<Molecule_and_Embedding *> maeb;

  mae.add (nullptr);
  mae[0];

  resizable_array_p<No_Reaction> nr;
  resizable_array_base<No_Reaction *> nrb;

  nr.add (nullptr);
  nr[0];

  resizable_array_p<Substitute_Atom> sa;
  resizable_array_base<Substitute_Atom *> sab;

  sab.add (nullptr);
  (void) sab[0];

  return;
}

#endif
