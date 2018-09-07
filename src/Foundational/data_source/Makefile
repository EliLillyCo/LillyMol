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
LIBS = $(LIB) -L$(IWPUB)/lib -liwcrex -liwrandom -lcmdline -liwmisc -liwstring $(LIBZ) -liwaray $(SYSTEM_LIBS)

########################## Template for libary #########################
# Change the project name
PROJECT_NAME = iwstring_data_source

LIBRARY_OBJECTS_LIST = iwstring_data_source.o string_data_source.o iwmmap.o _storage_reader_ram.o

SHARED_INC_FILES = iwmmap.h iwstring_data_source.h string_data_source.h

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
TIWDS_EXE = $(BUILD_DIR)/test_iwds
TIWDS_OBJECTS = $(BUILD_DIR)/tiwds.o $(BUILD_DIR)/iwstring_data_source.o 
$(BUILD_DIR)/tiwds.o: tiwds.cc iwstring_data_source.h
	$(CXX) -c -o $@ $< $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES)
$(BUILD_DIR)/test_iwds: $(TIWDS_OBJECTS)
	$(LD) -o $@ $(TIWDS_OBJECTS) $(LIBS) 
#~ EXECUTABLES += $(TIWDS_EXE)

################# Template for public executables ######################
TIWDS_STDIN_EXE = $(BUILD_DIR)/test_iwds_stdin
TIWDS_STDIN_OBJECTS = $(BUILD_DIR)/tiwds_stdin.o $(BUILD_DIR)/iwstring_data_source.o
$(BUILD_DIR)/tiwds_stdin.o: tiwds_stdin.cc iwstring_data_source.h
	$(CXX) -c -o $@ $< $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES)
$(BUILD_DIR)/test_iwds_stdin: $(TIWDS_STDIN_OBJECTS)
	$(LD) -o $@ $(TIWDS_STDIN_OBJECTS) $(LIBS) 
#~ EXECUTABLES += $(TIWDS_STDIN_EXE)

################# Template for public executables ######################
TEST_IWMMAP_EXE = $(BUILD_DIR)/test_iwmmap
TEST_IWMMAP_OBJECTS = $(BUILD_DIR)/test_iwmmap.o $(BUILD_DIR)/iwmmap.o $(BUILD_DIR)/_storage_reader_ram.o
$(BUILD_DIR)/test_iwmmap.o: test_iwmmap.cc iwmmap.h
	$(CXX) -c -o $@ $< $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES)
$(BUILD_DIR)/test_iwmmap: $(TEST_IWMMAP_OBJECTS)
	$(LD) -o $@ $(TEST_IWMMAP_OBJECTS) $(LIBS) 
#~ EXECUTABLES += $(TEST_IWMMAP_EXE)
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

# Rule for copying executables
#copy_exe:
#	$(CP) $(EXECUTABLES) $(EXEDIR)
	
clean:
	@$(RM) $(BUILD_DIR)/*.o || true

veryclean: clean
	@$(RM) $(BUILD_DIR)/*.a || true
#	@$(RM) $(EXECUTABLES) || true
