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
LIBS = $(LIB) -lcmdline -liwrandom \
	-liwstring_data_source -liwcrex -liwstring -liwaray

BITS_IN_COMMON_OPTIONS = $(COMPILER) $(PROF) $(CPU) -O3 -unroll-loop=2 $(INC)

########################## Template for libary #########################
# Change the project name
PROJECT_NAME = iwbits

LIBRARY_OBJECTS_LIST = iwbits.o du_bin2ascii.o fastprint.o hexbits.o _resizable_array_iwbits.o fixed_size_counted_fingerprint.o

LIBRARY_C_OBJECTS_LIST = bits_in_common.o
 
SHARED_INC_FILES = iwbits.h dy_fingerprint.h

LIBRARY_OBJECTS = $(addprefix $(BUILD_DIR)/, $(LIBRARY_OBJECTS_LIST))

LIBRARY_C_OBJECTS = $(addprefix $(BUILD_DIR)/, $(LIBRARY_C_OBJECTS_LIST))

LIBRARY_STEM = lib$(PROJECT_NAME)

LIBRARY_NAME = $(BUILD_DIR)/$(LIBRARY_STEM).a

# Rule for making library obj
$(LIBRARY_OBJECTS):
	$(CXX) -c $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES) $(basename $(notdir $@)).cc -o $@

$(LIBRARY_C_OBJECTS):
	$(CC) -c $(BITS_IN_COMMON_OPTIONS) $(basename $(notdir $@)).c -o $@ $<

# Rule for making library
$(LIBRARY_NAME): $(LIBRARY_OBJECTS) $(LIBRARY_C_OBJECTS)
	$(AR) r $@ $(LIBRARY_OBJECTS) $(LIBRARY_C_OBJECTS)

################# Template for public executables ######################
TBITS_EXE = $(BUILD_DIR)/test_bits
TBITS_OBJECTS = $(BUILD_DIR)/tbits.o $(LIBRARY_OBJECTS) $(LIBRARY_C_OBJECTS)
$(BUILD_DIR)/tbits.o: tbits.cc
	$(CXX) -c -o $@ $< $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES)
$(BUILD_DIR)/test_bits: $(TBITS_OBJECTS)
	$(LD) -o $@ $(PROF) $(TBITS_OBJECTS) $(LIBS) $(LIBTEMPLATE) $(SYSTEM_LIBS) 
#~ EXECUTABLES += $(TBITS_EXE)

################# Template for public executables ######################
TBIC_EXE = $(BUILD_DIR)/test_bic
TBIC_OBJECTS = $(BUILD_DIR)/tbic.o $(LIBRARY_OBJECTS) $(LIBRARY_C_OBJECTS)
$(BUILD_DIR)/tbic.o: tbic.cc
	$(CXX) -c -o $@ $< $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES)
$(BUILD_DIR)/test_bic: $(TBIC_OBJECTS)
	$(LD) -o $@ $(PROF) $(TBIC_OBJECTS) $(LIBS) $(LIBTEMPLATE) $(SYSTEM_LIBS) 
#~ EXECUTABLES += $(TBIC_EXE)

################# Template for public executables ######################
TEST_EQ_EXE = $(BUILD_DIR)/test_eq
TEST_EQ_OBJECTS = $(BUILD_DIR)/testeq.o $(BUILD_DIR)/iwbits.o $(BUILD_DIR)/du_bin2ascii.o $(BUILD_DIR)/bits_in_common.o
$(BUILD_DIR)/testeq.o: testeq.cc
	$(CXX) -c -o $@ $< $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES)
$(BUILD_DIR)/test_eq: $(TEST_EQ_OBJECTS)
	$(LD) -o $@ $(PROF) $(TEST_EQ_OBJECTS) $(LIBS) $(LIBTEMPLATE) $(SYSTEM_LIBS) 
#~ EXECUTABLES += $(TEST_EQ_EXE)

################# Template for public executables ######################
TDAYLIGHT_EXE = $(BUILD_DIR)/test_daylight
TDAYLIGHT_OBJECTS = $(BUILD_DIR)/tdaylight.o $(LIBRARY_OBJECTS) $(LIBRARY_C_OBJECTS)
$(BUILD_DIR)/tdaylight.o: tdaylight.cc
	$(CXX) -c -o $@ $< $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES)
$(BUILD_DIR)/test_daylight: $(TDAYLIGHT_OBJECTS)
	$(LD) -o $@ $(PROF) $(TDAYLIGHT_OBJECTS) $(LIBS) $(LIBTEMPLATE) $(SYSTEM_LIBS) 
#~ EXECUTABLES += $(TDAYLIGHT_EXE)

################# Template for public executables ######################
TEST_DAYLIGHT_EXE = $(BUILD_DIR)/test_reading_daylight_represenation
TEST_DAYLIGHT_OBJECTS = $(BUILD_DIR)/test_reading_daylight_represenation.o $(LIBRARY_OBJECTS) $(LIBRARY_C_OBJECTS)
$(BUILD_DIR)/test_reading_daylight_represenation.o: test_reading_daylight_represenation.cc
	$(CXX) -c -o $@ $< $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES)
$(BUILD_DIR)/test_reading_daylight_represenation: $(TEST_DAYLIGHT_OBJECTS)
	$(LD) -o $@ $(PROF) $(TEST_DAYLIGHT_OBJECTS) $(LIBS) -laccumulator -liwcrex $(LIBTEMPLATE) $(SYSTEM_LIBS) $(LIBZ)
#~ EXECUTABLES += $(TEST_DAYLIGHT_EXE)

################# Template for public executables ######################
DETERMINE_FIRST_BIT_EXE = $(BUILD_DIR)/test_determine_first_bit
DETERMINE_FIRST_BIT_OBJECTS = $(BUILD_DIR)/determine_first_bit.o $(LIBRARY_OBJECTS) $(LIBRARY_C_OBJECTS)
$(BUILD_DIR)/determine_first_bit.o: determine_first_bit.cc
	$(CXX) -c -o $@ $< $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES)
$(BUILD_DIR)/test_determine_first_bit: $(DETERMINE_FIRST_BIT_OBJECTS)
	$(LD) -o $@ $(PROF) $(DETERMINE_FIRST_BIT_OBJECTS) $(LIBS) $(LIBTEMPLATE) $(SYSTEM_LIBS) 
#~ EXECUTABLES += $(DETERMINE_FIRST_BIT_EXE)

################# Template for public executables ######################
MKPRECOMPBIT_EXE = $(BUILD_DIR)/test_mkprecompbit
MKPRECOMPBIT_OBJECTS = $(BUILD_DIR)/mkprecompbit.o $(LIBRARY_OBJECTS) $(LIBRARY_C_OBJECTS)
$(BUILD_DIR)/mkprecompbit.o: mkprecompbit.cc
	$(CXX) -c -o $@ $< $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES)
$(BUILD_DIR)/test_mkprecompbit: $(MKPRECOMPBIT_OBJECTS)
	$(LD) -o $@ $(PROF) $(MKPRECOMPBIT_OBJECTS) $(LIBS) $(LIBTEMPLATE) $(SYSTEM_LIBS) 
#~ EXECUTABLES += $(MKPRECOMPBIT_EXE)

################# Template for public executables ######################
TOMDOMAN_FP_TO_BITS_EXE = $(BUILD_DIR)/test_tomdoman_fingerprints_to_bits
TOMDOMAN_FP_TO_BITS_OBJECTS = $(BUILD_DIR)/tomdoman_fingerprints_to_bits.o $(LIBRARY_OBJECTS) $(LIBRARY_C_OBJECTS)
$(BUILD_DIR)/tomdoman_fingerprints_to_bits.o: mkprecompbit.cc
	$(CXX) -c -o $@ $< $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES)
$(BUILD_DIR)/test_tomdoman_fingerprints_to_bits: $(TOMDOMAN_FP_TO_BITS_OBJECTS)
	$(LD) -o $@ $(PROF) $(TOMDOMAN_FP_TO_BITS_OBJECTS) $(LIBS) $(LIBTEMPLATE) $(SYSTEM_LIBS) 
#~ EXECUTABLES += $(TOMDOMAN_FP_TO_BITS_EXE)

###################### Template for build command ######################
.PHONY: copy_include
.PHONY: library
.PHONY: copy_library
.PHONY: exe
.PHONY: copy_exe
.PHONY: clean
.PHONY: veryclean

exe copy_exe:
	@echo "Nothing to do for $(PROJECT_NAME) target $@"
	
copy_include:
	$(CP) $(SHARED_INC_FILES) $(INCDIR)

# Rule for making static library
library: $(LIBRARY_NAME)

# Rule for installing static library
copy_library: copy_include library
	$(CP) $(LIBRARY_NAME) $(LIBDIR)
	
# Rule for making executables
#exe: $(EXECUTABLES)

# Rule for install executables
#copy_exe:
#	$(CP) $(EXECUTABLES) $(EXEDIR)
	
clean:
	@$(RM) $(BUILD_DIR)/*.o || true

veryclean: clean
	@$(RM) $(BUILD_DIR)/*.a || true
	@$(RM) $(EXECUTABLES) || true
