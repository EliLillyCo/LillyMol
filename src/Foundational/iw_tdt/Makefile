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
LIBS = $(LIB) -liwstring_data_source $(USRLIB) -lcmdline -liwstring -liwmisc -liwcrex -liwaray $(LIBZ)
NCLIBS = $(LIB) -liwstring_data_source $(USRLIB) -lcmdline_v2 -liwstring -liwmisc -liwcrex -liwaray $(LIBZ)

########################## Template for libary #########################
# Change the project name
PROJECT_NAME = iwtdt

LIBRARY_OBJECTS_LIST = iw_tdt.o _iwaray_p_tdt.o _iw_tdt_iwstring.o \
	_set_dataitem_value.o _add_dataitem.o _add_dataitem_float.o \
	_iw_tdt_int.o _iw_tdt_double.o _iw_tdt_float.o \
	iw_tdt_filter.o iw_tdt_string.o _resizable_array_filter.o

SHARED_INC_FILES = iw_tdt.h  iw_tdt_filter.h

LIBRARY_OBJECTS = $(addprefix $(BUILD_DIR)/, $(LIBRARY_OBJECTS_LIST))

LIBRARY_STEM = lib$(PROJECT_NAME)

LIBRARY_NAME = $(BUILD_DIR)/$(LIBRARY_STEM).a

# Rule for making library obj
$(LIBRARY_OBJECTS):
	$(CXX) -c $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES) $(basename $(notdir $@)).cc -o $@

# Rule for making library
$(LIBRARY_NAME): $(LIBRARY_OBJECTS)
	$(AR) r $@ $(LIBRARY_OBJECTS)

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
exe:
	@echo "Nothing to do for iw_tdt target $@"

# Rule for install executables
copy_exe:
	@echo "Nothing to do for iw_tdt target $@"
	
clean:
	@$(RM) $(BUILD_DIR)/*.o || true

veryclean: clean
	@$(RM) $(BUILD_DIR)/*.a || true

