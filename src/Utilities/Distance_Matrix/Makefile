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
LIBS = $(LIBZ) $(LIB)

########################## Template for libary #########################
# Change the project name
PROJECT_NAME = iwdistance_matrix

LIBRARY_OBJECTS_LIST = \
	_iwdistance_matrix_byte.o\
	_iwdistance_matrix_int.o\
	_iwdistance_matrix_short.o\
	_iwdistance_matrix_float.o\
	_iwdistance_matrix_double.o\
	_iwdistance_matrix_float_byte.o\
	_iwdistance_matrix_double_byte.o\
    _masquerading_as_byte_float.o\
    _masquerading_as_byte_double.o\
    _dm_to_nn_conditions_float.o\
    _dm_to_nn_conditions_int.o\
    _dm_to_nn_conditions_double.o\
    _dm_to_nn_conditions_uchar.o\
	iwdmsupport.o

SHARED_INC_FILES =  iwdmsupport.h IWDistanceMatrixBase.h Masquerading_as_Byte.h DM_to_NN_Conditions_.h
	
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


replace_exe:
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
#copy_exe: exe
#	$(CP) $(EXECUTABLES) $(EXEDIR)
	
clean:
	@$(RM) $(BUILD_DIR)/*.o || true

veryclean: clean
	@$(RM) $(BUILD_DIR)/*.a || true

