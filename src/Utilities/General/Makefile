
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
IWPROGRAMMES ?= $(PWD)/../../..
BUILD_DIR ?= .
UNAME ?= Linux-gcc-6.2.0
include $(IWPROGRAMMES)/makefile.public.$(UNAME)

########################## Shared variables ############################
LIBS =  $(LIB) -lcmdline -liwrandom -laccumulator\
	-liwstring_data_source -liwcrex -liwstring \
	-liwbits $(LIBZ) -liwaray

# If we are using the new command line object
NCLIBS =  $(LIB) -lcmdline_v2 -liwmisc -liwrandom -laccumulator\
	-liwstring_data_source -liwcrex -liwstring \
	-liwbits -liwaray\
	$(LIBZ) $(USRLIB) 

$(BUILD_DIR)/%.o : %.cc
	$(CXX) -c -o $@ $< $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES)
########################## Template for libary #########################
# Change the project name
PROJECT_NAME = iwnormalisation

LIBRARY_OBJECTS_LIST = normalisation.o

SHARED_INC_FILES =  normalisation.h
	
LIBRARY_OBJECTS = $(addprefix $(BUILD_DIR)/, $(LIBRARY_OBJECTS_LIST))

LIBRARY_STEM = lib$(PROJECT_NAME)

LIBRARY_NAME = $(BUILD_DIR)/$(LIBRARY_STEM).a

# Rule for making library obj
$(LIBRARY_OBJECTS):
	$(CXX) -c $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES) $(basename $(notdir $@)).cc -o $@

# Rule for making library
$(LIBRARY_NAME): $(LIBRARY_OBJECTS)
	$(AR) r $@ $(LIBRARY_OBJECTS)

################# Template for public executables ######################
UNIQUE_ROWS_EXE = $(BUILD_DIR)/unique_rows
UNIQUE_ROWS_OBJECTS = $(BUILD_DIR)/unique_rows.o $(BUILD_DIR)/iwtokeniser.o
$(BUILD_DIR)/unique_rows: $(UNIQUE_ROWS_OBJECTS)
	$(LD) -o $@ $(LDFLAGS) $(PROF) $(UNIQUE_ROWS_OBJECTS) $(NCLIBS) $(SYSTEM_LIBS) $(IWMALLOC)
EXECUTABLES += $(UNIQUE_ROWS_EXE)

################# Template for public executables ######################
IWCUT_EXE = $(BUILD_DIR)/iwcut
IWCUT_OBJECTS = $(BUILD_DIR)/iwcut.o
$(BUILD_DIR)/iwcut: $(IWCUT_OBJECTS)
	$(LD) -o $@ $(IWCUT_OBJECTS) $(LIBS) -lcmdline -liwstring_data_source -liwmisc -liwstring -liwcrex -liwaray $(SYSTEM_LIBS) 
EXECUTABLES += $(IWCUT_EXE)

################# Template for public executables ######################
FSQ_EXE = $(BUILD_DIR)/fetch_smiles_quick
FSQ_OBJECTS = $(BUILD_DIR)/fetch_smiles_quick.o
$(BUILD_DIR)/fetch_smiles_quick: $(FSQ_OBJECTS)
	$(LD) -o $@ $(LDFLAGS) $(PROF) $(FSQ_OBJECTS) $(LIBS) -liwmisc $(SYSTEM_LIBS)
EXECUTABLES += $(FSQ_EXE)

################# Template for public executables ######################
CONCAT_FILES_EXE = $(BUILD_DIR)/concat_files
CONCAT_FILES_OBJECTS = $(BUILD_DIR)/concat_files.o $(BUILD_DIR)/iwtokeniser.o
$(BUILD_DIR)/concat_files: $(CONCAT_FILES_OBJECTS)
	$(LD) -o $@ $(CONCAT_FILES_OBJECTS) $(LIBS) -liwstring -lcmdline -liwstring_data_source -liwcrex -liwstring -liwaray -liwmisc $(SYSTEM_LIBS) 
EXECUTABLES += $(CONCAT_FILES_EXE)

################# Template for public executables ######################
DESCRIPTORS_TO_FINGERPRINT_EXE = $(BUILD_DIR)/descriptors_to_fingerprint
DESCRIPTORS_TO_FINGERPRINT_OBJECTS = $(BUILD_DIR)/descriptors_to_fingerprint.o
$(BUILD_DIR)/descriptors_to_fingerprint: $(DESCRIPTORS_TO_FINGERPRINT_OBJECTS)
	$(LD) -o $@ $(PROF) $(DESCRIPTORS_TO_FINGERPRINT_OBJECTS) $(LIBS) -lsparsefp -liwmisc -liwbits $(SYSTEM_LIBS)	
EXECUTABLES += $(DESCRIPTORS_TO_FINGERPRINT_EXE)

###################### Template for build command ######################
.PHONY: copy_include
.PHONY: library
.PHONY: copy_library
.PHONY: exe
.PHONY: copy_exe
.PHONY: clean
.PHONY: veryclean

copy_include:
	$(CP) $(SHARED_INC_FILES) $(INCDIR)

# Rule for making static library
library: $(LIBRARY_NAME)

# Rule for installing static library
copy_library: copy_include library
	$(CP) $(LIBRARY_NAME) $(LIBDIR)
	
# Rule for making executables
exe: $(EXECUTABLES)

# Rule for install executables
copy_exe: exe
	$(CP) $(EXECUTABLES) $(EXEDIR)
	
clean:
	@$(RM) $(BUILD_DIR)/*.o || true

veryclean: clean
	@$(RM) $(BUILD_DIR)/*.a || true
	@$(RM) $(EXECUTABLES) || true
