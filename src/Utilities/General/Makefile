
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

################# Template for public executables ######################
AVERAGE_EXE = $(BUILD_DIR)/average
AVERAGE_OBJECTS = $(BUILD_DIR)/average.o $(BUILD_DIR)/iwtokeniser.o
$(BUILD_DIR)/average: $(AVERAGE_OBJECTS)
	$(LD) -o $@ $(LDFLAGS) $(PROF) $(AVERAGE_OBJECTS) -liwmisc $(LIBS) $(SYSTEM_LIBS) $(IWMALLOC)
EXECUTABLES += $(AVERAGE_EXE)

################# Template for public executables ######################
DFILEFILTER_EXE = $(BUILD_DIR)/dfilefilter
DFILEFILTER_OBJECTS = $(BUILD_DIR)/dfilefilter.o
$(BUILD_DIR)/dfilefilter: $(DFILEFILTER_OBJECTS)
	$(LD) -o $@ $(DFILEFILTER_OBJECTS) $(LIBS) -liwmisc
EXECUTABLES += $(DFILEFILTER_EXE)

################# Template for public executables ######################
DISTRIBUTION_EXE = $(BUILD_DIR)/distribution
DISTRIBUTION_OBJECTS = $(BUILD_DIR)/distribution.o
$(BUILD_DIR)/distribution: $(DISTRIBUTION_OBJECTS)
	$(LD) -o $@ $(DISTRIBUTION_OBJECTS) $(LIBS) -liwmisc  $(SYSTEM_LIBS)
EXECUTABLES += $(DISTRIBUTION_EXE)

################# Template for public executables ######################
FETCH_SDF_QUICK_EXE = $(BUILD_DIR)/fetch_sdf_quick
FETCH_SDF_QUICK_OBJECTS = $(BUILD_DIR)/fetch_sdf_quick.o
$(BUILD_DIR)/fetch_sdf_quick: $(FETCH_SDF_QUICK_OBJECTS)
	$(LD) -o $@ $(LDFLAGS) $(PROF) $(FETCH_SDF_QUICK_OBJECTS) $(LIBS) -liwmisc $(SYSTEM_LIBS)
EXECUTABLES += $(FETCH_SDF_QUICK_EXE)

################# Template for public executables ######################
IWSPLIT_EXE = $(BUILD_DIR)/iwsplit
IWSPLIT_OBJECTS = $(BUILD_DIR)/iwsplit.o
$(BUILD_DIR)/iwsplit: $(IWSPLIT_OBJECTS)
	$(LD) -o $@ $(PROF) $(IWSPLIT_OBJECTS) $(NCLIBS) -liwmisc -liwcrex $(SYSTEM_LIBS) $(LIBZ)
EXECUTABLES += $(IWSPLIT_EXE)

################# Template for public executables ######################
NORMALISE_EXE = $(BUILD_DIR)/normalise
NORMALISE_OBJECTS = $(BUILD_DIR)/normalise.o
$(BUILD_DIR)/normalise: $(NORMALISE_OBJECTS)
	$(LD) -o $@ $(LDFLAGS) $(NORMALISE_OBJECTS) $(LIBS) -liwmisc $(LIBTEMPLATE) $(SYSTEM_LIBS) 
EXECUTABLES += $(NORMALISE_EXE)

################# Template for public executables ######################
NENV_EXE = $(BUILD_DIR)/notenoughvariance
NENV_OBJECTS = $(BUILD_DIR)/notenoughvariance.o
$(BUILD_DIR)/notenoughvariance: $(NENV_OBJECTS)
	$(LD) -o $@ $(LDFLAGS) $(NENV_OBJECTS) $(LIBS) $(LIBTEMPLATE) $(SYSTEM_LIBS)
EXECUTABLES += $(NENV_EXE)

################# Template for public executables ######################
RANDOM_RECORDS_EXE = $(BUILD_DIR)/random_records
RANDOM_RECORDS_OBJECTS = $(BUILD_DIR)/random_records.o
$(BUILD_DIR)/random_records: $(RANDOM_RECORDS_OBJECTS)
	$(LD) -o $@ $(RANDOM_RECORDS_OBJECTS) $(LIBS) -laccumulator -lcmdline -liwstring_data_source -liwstring -liwcrex -liwaray $(SYSTEM_LIBS) 
EXECUTABLES += $(RANDOM_RECORDS_EXE)

################# Template for public executables ######################
TCOUNT_EXE = $(BUILD_DIR)/tcount
TCOUNT_OBJECTS = $(BUILD_DIR)/tcount.o $(BUILD_DIR)/tokens_in_quoted_string.o
$(BUILD_DIR)/tcount: $(TCOUNT_OBJECTS)
	$(LD) -o $@ $(LDFLAGS) $(PROF) $(TCOUNT_OBJECTS) $(LIBS) $(LIBTEMPLATE) $(SYSTEM_LIBS)
EXECUTABLES += $(TCOUNT_EXE)

################# Template for public executables ######################
WHATSMISSING_EXE = $(BUILD_DIR)/whatsmissing
WHATSMISSING_OBJECTS = $(BUILD_DIR)/whatsmissing.o
$(BUILD_DIR)/whatsmissing: $(WHATSMISSING_OBJECTS)
	$(LD) -o $@ $(LDFLAGS) $(PROF) $(WHATSMISSING_OBJECTS) $(LIBS) -liwmisc $(SYSTEM_LIBS) $(IWMALLOC)
EXECUTABLES += $(WHATSMISSING_EXE)

################# Template for public executables ######################
REARRANGE_COLUMNS_EXE = $(BUILD_DIR)/rearrange_columns
REARRANGE_COLUMNS_OBJECTS = $(BUILD_DIR)/rearrange_columns.o $(BUILD_DIR)/valid_descriptor_name.o
$(BUILD_DIR)/rearrange_columns: $(REARRANGE_COLUMNS_OBJECTS)
	$(LD) -o $@ $(REARRANGE_COLUMNS_OBJECTS) $(LIBS) -lcmdline -liwstring_data_source -liwstring -liwcrex -liwaray
EXECUTABLES += $(REARRANGE_COLUMNS_EXE)

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
