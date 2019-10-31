# Copyright 2018 Eli Lilly and Company 
# 
# 
# Licensed under the Apache License, Version 2.0 (the "License"); 
# you may not use this file except in compliance with the License. 
# You may obtain a copy of the License at  
# 
#     http://www.apache.org/licenses/LICENSE-2.0  
# 
# Unless required by applicable law or agreed to in writing, software 
# distributed under the License is distributed on an "AS IS" BASIS, 
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
# See the License for the specific language governing permissions and 
# limitations under the License. 
########################################################################
IWPROGRAMMES ?= $(PWD)/../..
BUILD_DIR ?= .
UNAME ?= Linux-gcc-6.2.0
include $(IWPROGRAMMES)/makefile.public.$(UNAME)

FLAGS_FOR_TEMPLATE_INSTANTIATORS = $(CXXFLAGS) -fno-implicit-templates $(PROJ) $(DEBUG) $(PROF) $(COMPILER) $(CPU) $(OPT) $(INC) $(TEMPLATES)
NEW_INC_DIR = $(INCDIR)/mol 
########################## Template for libary #########################
PROJECT_NAME = iwtemplate

LIBRARY_OBJECTS_LIST = \
	_iwaray_set_of_atoms.o\
	_resizable_array_const_set_of_atoms.o\
	_resizable_array_p_ring.o\
	_resizable_array_ring.o\
	_resizable_array_p_element.o\
	_resizable_array_element.o\
	_resizable_array_atom.o\
	_resizable_array_bond.o\
	_resizable_array_p_bond.o\
	_resizable_array_subs.o\
	_iwaray_ssq.o\
	_iwaray_csubs.o\
	_iwaray_sym.o\
	_iwaray_p_rs.o\
	_iwaray_p_rmele.o\
	_resizable_array_molecule.o\
	_resizable_array_lors.o\
	_resizable_array_target.o\
	_resizable_array_molecule_to_match.o\
	_resizable_array_p_etr.o\
	_resizable_array_p_pthmsg.o\
	_resizable_array_p_beep.o\
	_resizable_array_p_ctb.o\
	_resizable_array_p_ema.o\
	_resizable_array_p_cs.o\
	_resizable_array_atom_list.o\
	_iwaray_ofsmtype.o\
	_iwaray_coords.o\
	_iwaray_path_scoring.o\
	_vector_float.o\
	_vector_double.o\
	_mdl_molecule.o\
	_resizable_array_possible_lactim.o\
	_collection_template_uint.o\
	_collection_template_float.o


SHARED_INC_FILES = molecule.h atom.h element.h bond.h bond_list.h \
	chiral_centre.h collection_template.h set_of_atoms.h coordinates.h \
	space_vector.h aromatic.h istream_and_type.h rwmolecule.h \
	ostream_and_type.h dihedral.h path.h iwstandard.h ematch.h \
	etrans.h iwmtypes.h charge_assigner.h iw_vdw.h iwrcb.h \
	iwrnm.h moleculej.h smiles.h atom_typing.h donor_acceptor.h \
	mpr.h temp_detach_atoms.h qry_wstats.h qry_wcharge.h substructure.h \
	target.h rmele.h atom_alias.h output.h rotbond_common.h \
	molecule_to_query.h mdl_file_data.h mdl_atom_record.h allowed_elements.h \
	numass.h is_actually_chiral.h path_scoring.h mdl.h mdl_molecule.h \
	moleculed.h toggle_kekule_form.h reaction_match_conditions.h iwreaction.h \
	rxn_file.h rwsubstructure.h charge_calculation.h misc2.h tmpsssr.h \
	pearlman.h marvin.h symmetry.h
	
LIBRARY_OBJECTS = $(addprefix $(BUILD_DIR)/, $(LIBRARY_OBJECTS_LIST))

LIBRARY_STEM = lib$(PROJECT_NAME)

LIBRARY_NAME = $(BUILD_DIR)/$(LIBRARY_STEM).a

# Rule for making library obj
$(LIBRARY_OBJECTS):
	$(CXX) -c $(FLAGS_FOR_TEMPLATE_INSTANTIATORS) $(basename $(notdir $@)).cc -o $@

# Rule for making library
$(LIBRARY_NAME): $(LIBRARY_OBJECTS)
	$(AR) r $@ $(LIBRARY_OBJECTS)

########################## Template for libary #########################
PROJECT_NAME_1 = iwmolecule

MOLECULE_IO_OBJECTS = mdl.o mdl_v30.o pdb.o mmod.o smi.o msi.o cif.o\
	tripos.o charmm.o mrk.o ostream_and_type.o output.o marvin.o 

SMILES_OBJECTS = smiles.o smiles_support.o

REACTION_OBJECTS = iwreaction.o iwreaction2.o _iwaray_rxn.o reaction_iterator.o rxn_dihedral_angle.o rxn_bond_length.o\
	rxn_bond_angle.o rxn_matched_sors.o reaction_match_conditions.o rxn_3d_replace.o reaction_from_smirks.o\
	rxnfile.o rxnfile2.o 

SUBSTRUCTURE_OBJECTS = csubstructure.o iwsubstructure.o rwsubstructure.o \
	substructure_a.o substructure_spec.o\
	target.o\
	substructure_results.o\
	substructure_env.o\
	substructure_nmab.o\
	ss_bonds.o\
	ss_ring.o\
	ss_ring_sys.o\
	ss_ring_base.o\
	ss_atom_env.o\
	element_hits_needed.o\
	tokenise_atomic_smarts.o\
	substructure_chiral.o\
	atom_alias.o\
	molecule_to_query.o\
	mdl_file_data.o\
	mdl_atom_record.o\
	parse_smarts_tmp.o\
	isis_link_atom.o\
	mdl_molecule.o\
	ISIS_Atom_List.o\
	charge_calculation.o\
	_resizable_array_mdl_atom.o\
	_resizable_array_scc.o\
	_substructure.o

C1 = molecule.o\
	path_scoring.o\
	moleculeb.o\
	bond_list.o\
	atom.o\
	bond.o\
	element.o\
	iwrnm.o\
	pearlman.o\
	iwrcb.o\
	moleculer.o\
	aromatic.o\
	temp_detach_atoms.o\
	moleculeh.o\
	path.o\
	ring_bond_iterator.o\
	frag.o\
	moleculed.o\
	molecule_smarts.o\
	$(SMILES_OBJECTS)\
	$(SUBSTRUCTURE_OBJECTS)\
	$(REACTION_OBJECTS)\
	set_of_atoms.o\
	rwmolecule.o\
	symmetry.o\
	$(MOLECULE_IO_OBJECTS)\
	unique.o\
	symm_class_can_rank.o\
	dihedral.o\
	coordinates.o\
	coordinates_double.o\
	cis_trans_bond.o\
	ematch.o\
	_istream_and_type.o\
	chiral_centre.o\
	toggle_kekule_form.o\
	atom_typing.o\
	donor_acceptor.o\
	charge_assigner.o\
	iwqry_wstats.o\
	_qry_wstats.o\
	_qry_wcharge.o\
	allowed_elements.o\
	is_actually_chiral.o\
	numass.o\
	inchi_dummy.o


LIBRARY_OBJECTS_LIST_1 = \
	standardise.o \
	rmele.o etrans.o \
	van_der_waals.o \
	mpr.o \
	$(C1) \
	misc2.o \
	careful_frag.o \
	rotbond_common.o
	

LIBRARY_OBJECTS_1 = $(addprefix $(BUILD_DIR)/, $(LIBRARY_OBJECTS_LIST_1))

LIBRARY_STEM_1 = lib$(PROJECT_NAME_1)

LIBRARY_NAME_1 = $(BUILD_DIR)/$(LIBRARY_STEM_1).a

# Rule for making library obj
$(LIBRARY_OBJECTS_1):
	$(CXX) -c $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES) $(basename $(notdir $@)).cc -o $@

# Fortran code in the library: override the general cxx compiling rule
$(BUILD_DIR)/u3b.o: u3b.f
	$(FC) -c -o $@ $< $(FFLAGS)
	
# Rule for making library
$(LIBRARY_NAME_1): $(LIBRARY_OBJECTS_1) $(BUILD_DIR)/u3b.o
	$(AR) r $@ $(LIBRARY_OBJECTS_1) $(BUILD_DIR)/u3b.o

###################### Template for build command ######################
.PHONY: copy_include
.PHONY: library
.PHONY: copy_library
.PHONY: exe
.PHONY: copy_exe
.PHONY: clean
.PHONY: veryclean
.PHONY: make_include_dir

make_include_dir:
	echo "Making ${NEW_INC_DIR} for Molecule_Lib"
	if [ ! -d $(NEW_INC_DIR) ] ; then mkdir $(NEW_INC_DIR) ; fi

exe copy_exe:
	@echo "Nothing to do for $(PROJECT_NAME) target $@"

copy_include: make_include_dir
	$(CP) $(SHARED_INC_FILES) $(NEW_INC_DIR)

# Rule for making static library
library: $(LIBRARY_NAME) $(LIBRARY_NAME_1)

# Rule for installing static library
copy_library: copy_include library
	$(CP) $(LIBRARY_NAME) $(LIBRARY_NAME_1) $(LIBDIR)
	
clean:
	@$(RM) $(BUILD_DIR)/*.o || true

veryclean: clean
	@$(RM) $(BUILD_DIR)/*.a || true

