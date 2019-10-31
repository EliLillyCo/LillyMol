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

########################## Shared variables ############################
ML = $(INCDIR)/mol

LIBS = $(LIB) -liwmolecule -liwtemplate -liwrandom -liwstring_data_source -liwcrex -laccumulator\
	-liwbits -liwmisc -lcmdline -liwstring -liwaray -lxmlParser $(LIBZ)
GFPLIB = -lgfp -liwmolecule -lsparsefp -liwtdt
MOLECULE_H = $(ML)/molecule.h $(ML)/element.h $(ML)/atom.h $(ML)/bond.h $(ML)/path.h $(ML)/iwrcb.h $(ML)/coordinates.h $(ML)/set_of_atoms.h

$(BUILD_DIR)/%.o : %.cc
	$(CXX) -c -o $@ $< $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES)

################# Template for public executables ######################
COMMON_NAMES_EXE = $(BUILD_DIR)/common_names
COMMON_NAMES_OBJECTS = $(BUILD_DIR)/common_names.o
$(BUILD_DIR)/common_names.o: common_names.cc
$(BUILD_DIR)/common_names: $(COMMON_NAMES_OBJECTS) 
	$(LD) -o $@ $(COMMON_NAMES_OBJECTS) $(LIBS) $(SYSTEM_LIBS)
EXECUTABLES += $(COMMON_NAMES_EXE)

################# Template for public executables ######################
UNIQUE_MOLECULE_EXE = $(BUILD_DIR)/unique_molecules
UNIQUE_MOLECULE_OBJECTS = $(BUILD_DIR)/unique_molecules.o
$(BUILD_DIR)/unique_molecules.o: unique_molecules.cc smiles_hash.h 
$(BUILD_DIR)/unique_molecules: $(UNIQUE_MOLECULE_OBJECTS)
	$(LD) -o $@ $(UNIQUE_MOLECULE_OBJECTS) $(LIBS) -liwtdt $(SYSTEM_LIBS) 
EXECUTABLES += $(UNIQUE_MOLECULE_EXE)

################# Template for public executables ######################
FILECONV_EXE = $(BUILD_DIR)/fileconv
FILECONV_OBJECTS = $(BUILD_DIR)/fileconv.o $(BUILD_DIR)/fix_structures.o \
	$(BUILD_DIR)/do_remove_duplicate_fragments.o $(BUILD_DIR)/known_fragment_data.o
$(BUILD_DIR)/fileconv.o: fileconv.cc fix_structures.h known_fragment_data.h
$(BUILD_DIR)/fileconv: $(FILECONV_OBJECTS)
	$(LD) $(LDFLAGS) -o $@ $(FILECONV_OBJECTS) $(LIBS) $(SYSTEM_LIBS)
EXECUTABLES += $(FILECONV_EXE)

################# Template for public executables ######################
TRXN_EXE = $(BUILD_DIR)/trxn
TRXN_OBJECTS = $(BUILD_DIR)/trxn.o
$(BUILD_DIR)/trxn.o: trxn.cc
$(BUILD_DIR)/trxn: $(TRXN_OBJECTS)
	$(LD) -o $@ $(TRXN_OBJECTS) -liwtdt $(LIBS) $(SYSTEM_LIBS)
EXECUTABLES += $(TRXN_EXE)

################# Template for public executables ######################
TSUBSTRUCTURE_EXE = $(BUILD_DIR)/tsubstructure
# Spliter from internal version
#TSUBSTRUCTURE_OBJECTS = $(BUILD_DIR)/tsubstructure.o
TSUBSTRUCTURE_OBJECTS = $(BUILD_DIR)/tsubstructure.o $(BUILD_DIR)/INT_tsubstructure_fp.o
$(BUILD_DIR)/tsubstructure.o: tsubstructure.cc 
$(BUILD_DIR)/tsubstructure: $(TSUBSTRUCTURE_OBJECTS)
	$(LD) $(STATIC) -o $@ $(TSUBSTRUCTURE_OBJECTS) -lsparsefp $(LIBS) $(SYSTEM_LIBS) -liwtdt
EXECUTABLES += $(TSUBSTRUCTURE_EXE)


################# Template for public executables ######################
RXN_SIGNATURE_EXE = $(BUILD_DIR)/rxn_signature
RXN_SIGNATURE_OBJECTS = $(BUILD_DIR)/rxn_signature.o
$(BUILD_DIR)/rxn_signature.o: rxn_signature.cc 
$(BUILD_DIR)/rxn_signature: $(RXN_SIGNATURE_OBJECTS)
	$(LD) -o $@ $(RXN_SIGNATURE_OBJECTS) $(LIBS) -lsparsefp -liwstring -liwmisc $(SYSTEM_LIBS)
EXECUTABLES += $(RXN_SIGNATURE_EXE)

################# Template for public executables ######################
RXN_STANDARDIZE_EXE = $(BUILD_DIR)/rxn_standardize
RXN_STANDARDIZE_OBJECTS = $(BUILD_DIR)/rxn_standardize.o 
$(BUILD_DIR)/rxn_standardize.o: rxn_standardize.cc 
$(BUILD_DIR)/rxn_standardize: $(RXN_STANDARDIZE_OBJECTS)
	$(LD) -o $@ $(RXN_STANDARDIZE_OBJECTS) $(LIBTEMPLATE) -lsparsefp -liwstring $(LIBS) $(SYSTEM_LIBS) $(OTHER_LIBS) $(IWMALLOC)
EXECUTABLES += $(RXN_STANDARDIZE_EXE)

################# Template for public executables ######################
RETROSYNTHESIS_EXE = $(BUILD_DIR)/retrosynthesis
RETROSYNTHESIS_OBJECTS = $(BUILD_DIR)/retrosynthesis.o
$(BUILD_DIR)/retrosynthesis.o: retrosynthesis.cc 
$(BUILD_DIR)/retrosynthesis: $(RETROSYNTHESIS_OBJECTS)
	$(LD) -o $@ $(RETROSYNTHESIS_OBJECTS) $(LIBS) $(SYSTEM_LIBS)
EXECUTABLES += $(RETROSYNTHESIS_EXE)

################# Template for public executables ######################
IWDEMERIT_EXE = $(BUILD_DIR)/iwdemerit
IWDEMERIT_OBJECTS = $(BUILD_DIR)/iwdemerit.o $(BUILD_DIR)/demerit.o $(BUILD_DIR)/substructure_demerits.o
$(BUILD_DIR)/iwdemerit: $(IWDEMERIT_OBJECTS)
	$(LD) -o $@ $(IWDEMERIT_OBJECTS) $(LIBS) $(SYSTEM_LIBS)
EXECUTABLES += $(IWDEMERIT_EXE)

################# Template for public executables ######################
PREFERRED_SMILES_EXE = $(BUILD_DIR)/preferred_smiles
PREFERRED_SMILES_OBJECTS = $(BUILD_DIR)/preferred_smiles.o 
$(BUILD_DIR)/preferred_smiles: $(PREFERRED_SMILES_OBJECTS)
	$(LD) -o $@  $(PREFERRED_SMILES_OBJECTS) $(LIBS) $(SYSTEM_LIBS) 
EXECUTABLES += $(PREFERRED_SMILES_EXE)

### Release 2
################# Template for public executables ######################
ACTIVITY_CONSISTENCY_EXE = $(BUILD_DIR)/activity_consistency
ACTIVITY_CONSISTENCY_OBJECTS = $(BUILD_DIR)/activity_consistency.o
$(BUILD_DIR)/activity_consistency.o: activity_consistency.cc
$(BUILD_DIR)/activity_consistency: $(ACTIVITY_CONSISTENCY_OBJECTS)
	$(LD) -o $@ $(ACTIVITY_CONSISTENCY_OBJECTS) $(LIBS) $(SYSTEM_LIBS)
EXECUTABLES += $(ACTIVITY_CONSISTENCY_EXE)

################# Template for public executables ######################
RANDOM_SMILES_EXE = $(BUILD_DIR)/random_smiles
RANDOM_SMILES_OBJECTS = $(BUILD_DIR)/random_smiles.o 
$(BUILD_DIR)/random_smiles.o: random_smiles.cc 
$(BUILD_DIR)/random_smiles: $(RANDOM_SMILES_OBJECTS)
	$(LD) -o $@ $(RANDOM_SMILES_OBJECTS) $(LIBS) $(SYSTEM_LIBS)
EXECUTABLES += $(RANDOM_SMILES_EXE)

################# Template for public executables ######################
ROTATABLE_BONDS_EXE = $(BUILD_DIR)/rotatable_bonds
ROTATABLE_BONDS_OBJECTS = $(BUILD_DIR)/rotatable_bonds.o
$(BUILD_DIR)/rotatable_bonds: $(ROTATABLE_BONDS_OBJECTS)
	$(LD) -o $@ $(ROTATABLE_BONDS_OBJECTS) $(LIBS) $(SYSTEM_LIBS)
EXECUTABLES += $(ROTATABLE_BONDS_EXE)


################# Template for public executables ######################
MSORT_EXE = $(BUILD_DIR)/msort
MSORT_OBJECTS = $(BUILD_DIR)/msort.o
$(BUILD_DIR)/msort.o: msort.cc
$(BUILD_DIR)/msort: $(MSORT_OBJECTS)
	$(LD) -o $@ $(MSORT_OBJECTS) $(LIBS) $(SYSTEM_LIBS) 
EXECUTABLES += $(MSORT_EXE)	

################# Template for public executables ######################
RXN_SUBSTRUCTURE_SEARCH_EXE = $(BUILD_DIR)/rxn_substructure_search
RXN_SUBSTRUCTURE_SEARCH_OBJECTS = $(BUILD_DIR)/rxn_substructure_search.o
$(BUILD_DIR)/rxn_substructure_search: $(RXN_SUBSTRUCTURE_SEARCH_OBJECTS)
	$(LD) -o $@ $(RXN_SUBSTRUCTURE_SEARCH_OBJECTS) $(LIBS) $(SYSTEM_LIBS) 
EXECUTABLES += $(RXN_SUBSTRUCTURE_SEARCH_EXE)

################# Template for public executables ######################
SMILES_MUTATION_EXE = $(BUILD_DIR)/smiles_mutation
SMILES_MUTATION_OBJECTS = $(BUILD_DIR)/smiles_mutation.o  $(BUILD_DIR)/count_atoms_in_smiles.o $(BUILD_DIR)/random_reactions.o
$(BUILD_DIR)/smiles_mutation.o: smiles_mutation.cc random_reactions.h
$(BUILD_DIR)/smiles_mutation: $(SMILES_MUTATION_OBJECTS)
	$(LD) -o $@ $(SMILES_MUTATION_OBJECTS) $(LIBS) $(SYSTEM_LIBS)
EXECUTABLES += $(SMILES_MUTATION_EXE)

################# Template for public executables ######################
TP1_EXE = $(BUILD_DIR)/tp_first_pass
TP1_OBJECTS = $(BUILD_DIR)/tp_first_pass.o 
$(BUILD_DIR)/tp_first_pass.o: tp_first_pass.cc
$(BUILD_DIR)/tp_first_pass: $(TP1_OBJECTS) 
	$(LD) -o $@ $(TP1_OBJECTS) $(LIBS) $(SYSTEM_LIBS)
EXECUTABLES += $(TP1_EXE)

################# Template for public executables ######################
MOL2QRY_EXE = $(BUILD_DIR)/mol2qry
MOL2QRY_OBJECTS = $(BUILD_DIR)/mol2qry.o 
$(BUILD_DIR)/mol2qry.o:mol2qry.cc 
$(BUILD_DIR)/mol2qry: $(MOL2QRY_OBJECTS)
	$(LD) -o $@ $(MOL2QRY_OBJECTS) $(LIBS) $(SYSTEM_LIBS) 
EXECUTABLES += $(MOL2QRY_EXE)

################# Template for public executables ######################
MOLECULAR_SCAFFOLD_EXE = $(BUILD_DIR)/molecular_scaffold
MOLECULAR_SCAFFOLD_OBJECTS = $(BUILD_DIR)/molecular_scaffold.o  
$(BUILD_DIR)/molecular_scaffold: $(MOLECULAR_SCAFFOLD_OBJECTS)
	$(LD) -o $@ $(MOLECULAR_SCAFFOLD_OBJECTS) $(LIBS) $(SYSTEM_LIBS) 
EXECUTABLES += $(MOLECULAR_SCAFFOLD_EXE)

################# Template for public executables ######################
MOLECULE_SUBSET_EXE = $(BUILD_DIR)/molecule_subset
MOLECULE_SUBSET_OBJECTS = $(BUILD_DIR)/molecule_subset.o 
$(BUILD_DIR)/molecule_subset: $(MOLECULE_SUBSET_OBJECTS)
	$(LD) -o $@ $(MOLECULE_SUBSET_OBJECTS) $(LIBS) $(SYSTEM_LIBS)
EXECUTABLES += $(MOLECULE_SUBSET_EXE)

################# Template for public executables ######################
RGROUP_EXE = $(BUILD_DIR)/rgroup
RGROUP_OBJECTS = $(BUILD_DIR)/rgroup.o 
$(BUILD_DIR)/rgroup.o: rgroup.cc
$(BUILD_DIR)/rgroup: $(RGROUP_OBJECTS) $(TEMPLATE_LIBRARY)
	$(LD) -o $@ $(RGROUP_OBJECTS) $(LIBS) $(SYSTEM_LIBS) 
EXECUTABLES += $(RGROUP_EXE)

################# Template for public executables ######################
RING_EXTRACTION_EXE = $(BUILD_DIR)/ring_extraction
RING_EXTRACTION_OBJECTS = $(BUILD_DIR)/ring_extraction.o $(BUILD_DIR)/ring_ext_rep.o $(BUILD_DIR)/reaction_duplicate.o
$(BUILD_DIR)/ring_extraction.o: ring_extraction.cc ring_ext_rep.h
$(BUILD_DIR)/ring_ext_rep.o: ring_ext_rep.cc ring_ext_rep.h
$(BUILD_DIR)/ring_extraction: $(RING_EXTRACTION_OBJECTS) 
	$(LD) -o $@ $(RING_EXTRACTION_OBJECTS) $(LIBS) $(SYSTEM_LIBS) 
EXECUTABLES += $(RING_EXTRACTION_EXE)

################# Template for public executables ######################
RING_TRIMMING_EXE = $(BUILD_DIR)/ring_trimming
RING_TRIMMING_OBJECTS = $(BUILD_DIR)/ring_trimming.o  
$(BUILD_DIR)/ring_trimming: $(RING_TRIMMING_OBJECTS)
	$(LD) -o $@ $(RING_TRIMMING_OBJECTS)  $(LIBS) -liwmisc $(SYSTEM_LIBS) 
EXECUTABLES += $(RING_TRIMMING_EXE)

################# Template for public executables ######################
SP3_FILTER_EXE = $(BUILD_DIR)/sp3_filter
SP3_FILTER_OBJECTS = $(BUILD_DIR)/sp3_filter.o
$(BUILD_DIR)/sp3_filter: $(SP3_FILTER_OBJECTS)
	$(LD) -o $@ $(SP3_FILTER_OBJECTS) $(LIBS) $(SYSTEM_LIBS) 
EXECUTABLES += $(SP3_FILTER_EXE)


################# Template for public executables ######################
TAUTOMER_GENERATION_EXE = $(BUILD_DIR)/tautomer_generation
TAUTOMER_GENERATION_OBJECTS = $(BUILD_DIR)/tautomer_generation.o 
$(BUILD_DIR)/tautomer_generation.o: tautomer_generation.cc 
$(BUILD_DIR)/tautomer_generation: $(TAUTOMER_GENERATION_OBJECTS)
	$(LD) -o $@ $(TAUTOMER_GENERATION_OBJECTS) $(LIBS) $(SYSTEM_LIBS) 
EXECUTABLES += $(TAUTOMER_GENERATION_EXE)

### Release 4
################# Template for public executables ###################### 
DBF_EXE = $(BUILD_DIR)/dbf
DBF_OBJECTS  = $(BUILD_DIR)/dbf.o
$(BUILD_DIR)/dbf: $(DBF_OBJECTS)
	$(LD) -o $@ $(DBF_OBJECTS) -lsparsefp $(LIBS) -liwhistogram 
EXECUTABLES += $(DBF_EXE)

################# Template for public executables ######################
HYDROPHOBIC_SECTIONS_EXE = $(BUILD_DIR)/hydrophobic_sections
HYDROPHOBIC_SECTIONS_OBJECTS = $(BUILD_DIR)/hydrophobic_sections.o
$(BUILD_DIR)/hydrophobic_sections.o:hydrophobic_sections.cc 
$(BUILD_DIR)/hydrophobic_sections: $(HYDROPHOBIC_SECTIONS_OBJECTS)
	$(LD) -o $@ $(HYDROPHOBIC_SECTIONS_OBJECTS) $(LIBS) $(SYSTEM_LIBS)
EXECUTABLES += $(HYDROPHOBIC_SECTIONS_EXE)

################# Template for public executables ###################### 
IWECFP_EXE = $(BUILD_DIR)/iwecfp
IWECFP_OBJECTS = $(BUILD_DIR)/iwecfp.o
$(BUILD_DIR)/iwecfp: $(IWECFP_OBJECTS) 
	$(LD) -o $@ $(IWECFP_OBJECTS) -lsparsefp $(LIBS) -liwstring
EXECUTABLES += $(IWECFP_EXE)

################# Template for public executables ######################
MAKE_THESE_MOLECULES_EXE = $(BUILD_DIR)/make_these_molecules
MAKE_THESE_MOLECULES_OBJECTS = $(BUILD_DIR)/make_these_molecules.o  
$(BUILD_DIR)/make_these_molecules.o: make_these_molecules.cc 
$(BUILD_DIR)/make_these_molecules: $(MAKE_THESE_MOLECULES_OBJECTS) 
	$(LD) -o $@ $(MAKE_THESE_MOLECULES_OBJECTS) $(LIBS) -liwmisc $(SYSTEM_LIBS)
EXECUTABLES += $(MAKE_THESE_MOLECULES_EXE)

################# Template for public executables ###################### 
MOLECULES_FROM_REAGENTS_EXE = $(BUILD_DIR)/molecules_from_reagents
MOLECULES_FROM_REAGENTS_OBJECTS = $(BUILD_DIR)/molecules_from_reagents.o 
$(BUILD_DIR)/molecules_from_reagents: $(MOLECULES_FROM_REAGENTS_OBJECTS)
	$(LD) -o $@ $(MOLECULES_FROM_REAGENTS_OBJECTS) $(LIBS) $(OPENMP)
EXECUTABLES += $(MOLECULES_FROM_REAGENTS_EXE)

################# Template for public executables ######################
REMOVE_AND_LABEL_EXE = $(BUILD_DIR)/remove_and_label
REMOVE_AND_LABEL_OBJECTS = $(BUILD_DIR)/remove_and_label.o 
$(BUILD_DIR)/remove_and_label.o: remove_and_label.cc 
$(BUILD_DIR)/remove_and_label: $(REMOVE_AND_LABEL_OBJECTS)
	$(LD) -o $@ $(REMOVE_AND_LABEL_OBJECTS) $(LIBS) $(SYSTEM_LIBS)
EXECUTABLES += $(REMOVE_AND_LABEL_EXE)

################# Template for public executables ######################
REMOVE_MATCHED_ATOMS_EXE = $(BUILD_DIR)/remove_matched_atoms
REMOVE_MATCHED_ATOMS_OBJECTS = $(BUILD_DIR)/remove_matched_atoms.o  $(BUILD_DIR)/spatially_common_matched_atoms.o 
$(BUILD_DIR)/remove_matched_atoms: $(REMOVE_MATCHED_ATOMS_OBJECTS)
	$(LD) -o $@ $(REMOVE_MATCHED_ATOMS_OBJECTS) $(LIBS) $(SYSTEM_LIBS) 
EXECUTABLES += $(REMOVE_MATCHED_ATOMS_EXE)

################# Template for public executables ###################### 
RING_FINGERPRINT_EXE = $(BUILD_DIR)/ring_fingerprint
RING_FINGERPRINT_OBJECTS = $(BUILD_DIR)/ring_fingerprint.o 
$(BUILD_DIR)/ring_fingerprint: $(RING_FINGERPRINT_OBJECTS)
	$(LD) -o $@ $(RING_FINGERPRINT_OBJECTS) -lsparsefp $(LIBS)
EXECUTABLES += $(RING_FINGERPRINT_EXE)

################# Template for public executables ######################
RING_REPLACEMENT_EXE = $(BUILD_DIR)/ring_replacement
RING_REPLACEMENT_OBJECTS = $(BUILD_DIR)/ring_replacement.o $(BUILD_DIR)/ring_ext_rep.o $(BUILD_DIR)/reaction_duplicate.o
$(BUILD_DIR)/ring_replacement.o: ring_replacement.cc
$(BUILD_DIR)/ring_replacement: $(RING_REPLACEMENT_OBJECTS) 
	$(LD) -o $@ $(RING_REPLACEMENT_OBJECTS)  $(LIBS) $(SYSTEM_LIBS) 
EXECUTABLES += $(RING_REPLACEMENT_EXE)

################# Template for public executables ######################
SUBSTITUTIONS_EXE = $(BUILD_DIR)/substitutions
SUBSTITUTIONS_OBJECTS = $(BUILD_DIR)/substitutions.o 
$(BUILD_DIR)/substitutions.o: substitutions.cc
$(BUILD_DIR)/substitutions: $(SUBSTITUTIONS_OBJECTS) 
	$(LD) -o $@ $(SUBSTITUTIONS_OBJECTS) $(LIBS) $(SYSTEM_LIBS)
EXECUTABLES += $(SUBSTITUTIONS_EXE)

################# Template for public executables ###################### 
TEMPERATURE_EXE = $(BUILD_DIR)/temperature
TEMPERATURE_OBJECTS = $(BUILD_DIR)/temperature.o 
$(BUILD_DIR)/temperature: $(TEMPERATURE_OBJECTS) 
	$(LD) -o $@ $(TEMPERATURE_OBJECTS) -lsparsefp $(LIBS)  
EXECUTABLES += $(TEMPERATURE_EXE)

################# Template for public executables ###################### 
TNASS_EXE = $(BUILD_DIR)/tnass 
TNASS_OBJECTS = $(BUILD_DIR)/tnass.o $(BUILD_DIR)/nass.o $(BUILD_DIR)/_resizable_array_nas.o
$(BUILD_DIR)/tnass: $(TNASS_OBJECTS) 
	$(LD) -o $@ $(TNASS_OBJECTS) -liwtdt -lsparsefp $(LIBS)
EXECUTABLES += $(TNASS_EXE)

#second path

################# Template for public executables ###################### 
DICER_EXE = $(BUILD_DIR)/dicer
DICER_OBJECTS = $(BUILD_DIR)/dicer.o  $(BUILD_DIR)/swap_kekule_forms.o $(BUILD_DIR)/count_atoms_in_smiles.o
$(BUILD_DIR)/dicer: $(DICER_OBJECTS) 
	$(LD) -o $@ $(DICER_OBJECTS) -lsparsefp $(LIBS) 
EXECUTABLES += $(DICER_EXE)

################# Template for public executables ###################### 
IWDESCR_EXE = $(BUILD_DIR)/iwdescr 
IWDESCR_OBJECTS = $(BUILD_DIR)/iwdescr.o $(BUILD_DIR)/qry_wcharge.o $(BUILD_DIR)/path_around_ring.o $(BUILD_DIR)/nvrtspsa.o
$(BUILD_DIR)/iwdescr: $(IWDESCR_OBJECTS)
	$(LD) -o $@ $(IWDESCR_OBJECTS) -lsparsefp $(LIBS) -laccumulator -liwstring 
EXECUTABLES += $(IWDESCR_EXE)

################# Template for public executables ###################### 
IWFP_EXE = $(BUILD_DIR)/iwfp         
IWFP_OBJECTS = $(BUILD_DIR)/iwfp.o $(BUILD_DIR)/iwmfingerprint.o
$(BUILD_DIR)/iwfp: $(IWFP_OBJECTS)
	$(LD) -o $@ $(IWFP_OBJECTS) -lsparsefp $(LIBS) 
EXECUTABLES += $(IWFP_EXE)

################## Template for public executables ###################### 
JWSADB_EXE = $(BUILD_DIR)/jwsadb    
JWSADB_OBJECTS = $(BUILD_DIR)/jw_sa_db_descriptor.o $(BUILD_DIR)/read_control_file.o $(BUILD_DIR)/comma_related_procedures.o $(BUILD_DIR)/jw_vector_array_operation.o $(BUILD_DIR)/surface_area_molvol.o $(BUILD_DIR)/smolvol.o
$(BUILD_DIR)/smolvol.o: smolvol.f
	$(FC) -c -o $@ $< $(FFLAGS)
$(BUILD_DIR)/jwsadb: $(JWSADB_OBJECTS)
	$(FORT_LD) -o $@ $(JWSADB_OBJECTS) $(FORTRAN_LIBRARY) $(LIBS) -lsparsefp -liwstring -liwmisc
EXECUTABLES += $(JWSADB_EXE)

################# Template for public executables ###################### 
MACCSKEYS_FN5_EXE = $(BUILD_DIR)/maccskeys    
MACCSKEYS_FN5_OBJECTS = $(BUILD_DIR)/maccskeys_fn5_main.o $(BUILD_DIR)/maccskeys_fn5.o
$(BUILD_DIR)/maccskeys: $(MACCSKEYS_FN5_OBJECTS)  
	$(LD) -o $@ $(MACCSKEYS_FN5_OBJECTS) $(LIBS)  -lsparsefp -liwstring -liwmisc
EXECUTABLES += $(MACCSKEYS_FN5_EXE)

################# Template for public executables ###################### 
MAV2_EXE = $(BUILD_DIR)/molecular_abstraction 
MAV2_OBJECTS = $(BUILD_DIR)/mav2.o $(BUILD_DIR)/molecular_abstraction_functions.o $(BUILD_DIR)/molecular_abstraction_specifications.o $(BUILD_DIR)/iwmfingerprint.o
$(BUILD_DIR)/molecular_abstraction: $(MAV2_OBJECTS)
	$(LD) -o $@ $(MAV2_OBJECTS) -lsparsefp  $(LIBS) 
EXECUTABLES += $(MAV2_EXE)

################# Template for public executables ###################### 
MOLECULAR_TRANSFORMATION_EXE = $(BUILD_DIR)/molecular_transformations
MOLECULAR_TRANSFORMATION_OBJECTS = $(BUILD_DIR)/molecular_transformations.o $(BUILD_DIR)/iwmfingerprint.o
$(BUILD_DIR)/molecular_transformations: $(MOLECULAR_TRANSFORMATION_OBJECTS)
	$(LD) -o $@ $(MOLECULAR_TRANSFORMATION_OBJECTS) $(LIBS) 
EXECUTABLES += $(MOLECULAR_TRANSFORMATION_EXE)

################# Template for public executables ###################### 
PERMUTATIONS_EXE = $(BUILD_DIR)/random_molecular_permutations 
PERMUTATIONS_OBJECTS = $(BUILD_DIR)/random_molecular_permutations.o $(BUILD_DIR)/swap_kekule_forms.o \
	$(BUILD_DIR)/random_reactions.o $(BUILD_DIR)/set_of_target_molecules.o $(BUILD_DIR)/linear_path_fingerprint.o $(BUILD_DIR)/maccskeys_fn5.o 
#~ 	$(MOLECULE_TOOLS_BIN_DIR)/random_reactions.o $(BUILD_DIR)/set_of_target_molecules.o $(BUILD_DIR)/linear_path_fingerprint.o $(BUILD_DIR)/maccskeys_fn5.o 
$(BUILD_DIR)/random_molecular_permutations: $(PERMUTATIONS_OBJECTS)
	$(LD) -o $@ $(PERMUTATIONS_OBJECTS) $(GFPLIB) $(LIBS) 
EXECUTABLES += $(PERMUTATIONS_EXE)	

################# Template for public executables ###################### 
RING_SUBSTITUTION_EXE = $(BUILD_DIR)/ring_substitution      
RING_SUBSTITUTION_OBJECTS = $(BUILD_DIR)/ring_substitution.o $(BUILD_DIR)/path_around_ring.o
$(BUILD_DIR)/ring_substitution: $(RING_SUBSTITUTION_OBJECTS) 
	$(LD) -o $@ $(RING_SUBSTITUTION_OBJECTS) -lsparsefp $(LIBS) 
EXECUTABLES += $(RING_SUBSTITUTION_EXE)	

################# Template for public executables ###################### 
RULE_OF_FIVE_EXE = $(BUILD_DIR)/rule_of_five 
RULE_OF_FIVE_OBJECTS = $(BUILD_DIR)/rule_of_five.o $(BUILD_DIR)/nvrtspsa.o
$(BUILD_DIR)/rule_of_five: $(RULE_OF_FIVE_OBJECTS)
	$(LD) -o $@ $(RULE_OF_FIVE_OBJECTS) -liwtdt $(LIBS) -liwstring 
EXECUTABLES += $(RULE_OF_FIVE_EXE)

################# Template for public executables ###################### 
RXN_FINGERPRINT_EXE = $(BUILD_DIR)/rxn_fingerprint 
RXN_FINGERPRINT_OBJECTS = $(BUILD_DIR)/rxn_fingerprint.o $(BUILD_DIR)/rxnfile3.o \
$(BUILD_DIR)/linear_path_fingerprint.o $(BUILD_DIR)/extended_connectivity_fp.o $(BUILD_DIR)/iwmfingerprint.o 
$(BUILD_DIR)/rxn_fingerprint: $(RXN_FINGERPRINT_OBJECTS) 
	$(LD) -o $@ $(RXN_FINGERPRINT_OBJECTS) -lsparsefp $(LIBS) 
EXECUTABLES += $(RXN_FINGERPRINT_EXE)

###################### Template for build command ######################
.PHONY: copy_include
.PHONY: library
.PHONY: copy_library
.PHONY: exe
.PHONY: copy_exe
.PHONY: clean
.PHONY: veryclean

copy_include library copy_library:
	@echo "Nothing to do for Molecule_Tools target $@"
	
# Rule for making executables
exe: $(EXECUTABLES)

# Rule for install executables
copy_exe: exe
	$(CP) $(EXECUTABLES) $(EXEDIR)
	
clean:
	@$(RM) $(BUILD_DIR)/*.o || true

veryclean: clean
	@$(RM) $(EXECUTABLES) || true
	
