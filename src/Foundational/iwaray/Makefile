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
PROJECT_NAME = iwaray

LIBRARY_OBJECTS_LIST = _resizable_array_int.o\
	_resizable_array_uint.o\
	_resizable_array_long.o\
	_resizable_array_long_long.o\
	_resizable_array_unsigned_long.o\
	_resizable_array_p_int.o\
	_resizable_array_float.o\
	_resizable_array_char.o\
	_resizable_array_char2.o\
	_resizable_array_uchar.o\
	_resizable_array_ushort.o\
	_resizable_array_double.o\
	_extending_resizable_array_int.o\
	_extending_resizable_array_ushort.o\
	_extending_resizable_array_long.o\
	_extending_resizable_array_long_long.o\
	_extending_resizable_array_uchar.o\
	_extending_resizable_array_float.o\
	_extending_resizable_array_double.o\
	_iwaray_int.o \
	_iwaray_float.o

SHARED_INC_FILES = iwaray.h

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
TEST_IWARAY_EXE = $(BUILD_DIR)/test_iwaray
TEST_IWARAY_OBJECTS = $(BUILD_DIR)/test_iwaray.o
$(BUILD_DIR)/test_iwaray.o: test_iwaray.cc iwaray.h iwaray_.h
	$(CC) -c $(FLAGS) $< -o $@
$(BUILD_DIR)/test_iwaray: $(TEST_IWARAY_OBJECTS) iwaray.h iwaray_.h $(TEMPLATE_LIBRARY)
	$(LD) -o $@ $(TEST_IWARAY_OBJECTS) $(LIB) -liwaray  $(SYSTEM_LIBS) -liwrandom -lm
#~ EXECUTABLES += $(TEST_IWARAY_EXE)

################# Template for public executables ######################
TEST_EXTENTING_EXE = $(BUILD_DIR)/test_extending 
TEST_EXTENDING_OBJECTS = $(BUILD_DIR)/test_extending.o
$(BUILD_DIR)/test_extending.o: test_extending.cc iwaray.h iwaray_.h extending_resizable_array_.h
	$(CXX) -c $(FLAGS) $< -o $@
$(BUILD_DIR)/test_extending: $(TEST_EXTENDING_OBJECTS)
	$(LD) -o $@ $(TEST_EXTENDING_OBJECTS) $(LIB) -liwaray  $(SYSTEM_LIBS) -liwrandom
#~ EXECUTABLES += $(TEST_EXTENTING_EXE)

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
