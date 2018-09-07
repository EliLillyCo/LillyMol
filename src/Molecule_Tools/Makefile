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
TSUBSTRUCTURE_OBJECTS = $(BUILD_DIR)/tsubstructure.o
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
	
