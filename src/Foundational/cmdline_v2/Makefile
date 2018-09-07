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
PROJECT_NAME = cmdline_v2

LIBRARY_OBJECTS_LIST = cmdline.o

SHARED_INC_FILES = cmdline_v2.h 

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
TEST_CMDLINE_EXE = $(BUILD_DIR)/test_iwcmdline
TEST_CMDLINE_OBJECTS = $(BUILD_DIR)/test_iwcmdline.o $(LIBRARY_OBJECTS)
$(BUILD_DIR)/test_iwcmdline.o: test_iwcmdline.cc cmdline_v2.h
	$(CXX) -c -o $@ $< $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES)
$(BUILD_DIR)/test_iwcmdline: $(TEST_CMDLINE_OBJECTS)
	$(LD) -o $@ $(TEST_CMDLINE_OBJECTS) $(LIB) -liwstring_data_source -liwcrex -liwmisc -liwstring -liwaray $(LIBZ)
#~ EXECUTABLES += $(TEST_CMDLINE_EXE)

################# Template for public executables ######################
TEST_IWCMDLINE_INTERACTIVE_EXE = $(BUILD_DIR)/test_iwcmdline_interactive
TEST_IWCMDLINE_INTERACTIVE_OBJECTS = $(BUILD_DIR)/test_iwcmdline_interactive.o $(LIBRARY_OBJECTS)
$(BUILD_DIR)/test_iwcmdline_interactive.o: test_iwcmdline_interactive.cc cmdline_v2.h
	$(CXX) -c -o $@ $< $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES)
$(BUILD_DIR)/test_iwcmdline_interactive: $(TEST_IWCMDLINE_INTERACTIVE_OBJECTS)
	$(LD) -o $@ $(TEST_IWCMDLINE_INTERACTIVE_OBJECTS) $(LIB) -liwstring_data_source -liwcrex -liwmisc -liwstring -liwaray $(LIBZ)
#~ EXECUTABLES += $(TEST_IWCMDLINE_INTERACTIVE_EXE)

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

# Rule for copying static library
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
	
