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
PROJECT_NAME = accumulator

LIBRARY_OBJECTS_LIST = _accumulator_int.o\
	_accumulator_unsigned_int.o\
	_accumulator_float.o\
	_accumulator_long.o\
	_accumulator_double.o\
	_accumulator_w_missing_float.o\
	_resizable_array_p_accumulator_float.o\
	_accumulator_int_v2.o\
	_accumulator_long_v2.o\
	_accumulator_float_v2.o\
	_accumulator_double_v2.o\
	KahanSum.o\

SHARED_INC_FILES = accumulator.h accumulator_v2.h kahan_sum.h

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


exe copy_exe:
	@echo "Nothing to do for $(PROJECT_NAME) target $@"

copy_include:
	$(CP) $(SHARED_INC_FILES) $(INCDIR)

# Rule for making static library
library: $(LIBRARY_NAME)

# Rule for copying static library
copy_library: copy_include library
	$(CP) $(LIBRARY_NAME) $(LIBDIR)
		
clean:
	@$(RM) $(BUILD_DIR)/*.o || true

veryclean: clean
	@$(RM) $(BUILD_DIR)/*.a || true
