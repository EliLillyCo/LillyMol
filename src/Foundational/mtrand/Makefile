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

########################## Template for libary #########################
# Change the project name
PROJECT_NAME = iwrandom

LIBRARY_OBJECTS_LIST = mtrand.o\
	iwrandom.o\
	normal_distribution.o\
	_random_number_between_int.o\
	_random_number_between_float.o\
	_random_number_between_double.o\
	_random_number_between_uint.o\
	_random_number_between_long.o\
	_random_number_between_unsigned_long.o\
	_random_number_between_long_long.o
 
SHARED_INC_FILES = iwrandom.h mtrand.h 

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
TIWRANDOM_EXE = $(BUILD_DIR)/test_iwrandom
TIWRANDOM_OBJECTS = $(BUILD_DIR)/tiwrandom.o $(LIBRARY_OBJECTS)
$(BUILD_DIR)/tiwrandom.o: tiwrandom.cc
	$(CXX) -c -o $@ $< $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES)
$(BUILD_DIR)/test_iwrandom: $(TIWRANDOM_OBJECTS)
	$(LD) -o $@ $(TIWRANDOM_OBJECTS) $(LIB) -lcmdline -liwstring -liwaray
#~ EXECUTABLES += $(TIWRANDOM_EXE)

################# Template for public executables ######################
TNORMAL_EXE = $(BUILD_DIR)/test_normal
TNORMAL_OBJECTS = $(BUILD_DIR)/tnormal.o $(LIBRARY_OBJECTS)
$(BUILD_DIR)/tnormal.o: tnormal.cc
	$(CXX) -c -o $@ $< $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES)
$(BUILD_DIR)/test_normal: $(TNORMAL_OBJECTS)
	$(LD) -o $@ $(TNORMAL_OBJECTS) $(LIB) -liwhistogram -laccumulator -lcmdline -liwstring -liwaray
#~ EXECUTABLES += $(TNORMAL_EXE)

################# Template for public executables ######################
TIWRANDOM_MTRAND_EXE = $(BUILD_DIR)/test_iwrandom_mtrand
TIWRANDOM_MTRAND_OBJECTS = $(BUILD_DIR)/tiwrandom_mtrand.o $(LIBRARY_OBJECTS)
$(BUILD_DIR)/tiwrandom_mtrand.o: tiwrandom_mtrand.cc
	$(CXX) -c -o $@ $< $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES)
$(BUILD_DIR)/test_iwrandom_mtrand: $(TIWRANDOM_MTRAND_OBJECTS)
	$(LD) -o $@ $(TIWRANDOM_MTRAND_OBJECTS) $(LIB) -lcmdline -liwstring -liwaray
#~ EXECUTABLES += $(TIWRANDOM_MTRAND_EXE)

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
#	@$(RM) $(EXECUTABLES) || true

