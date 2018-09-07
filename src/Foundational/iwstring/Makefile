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
LIBS = $(LIB) 

########################## Template for libary #########################
# Change the project name
PROJECT_NAME = iwstring

LIBRARY_OBJECTS_LIST = iwstring.o IWString_class.o iwgetline.o iwwrite.o \
	iwwrite_block.o _iwaray_p_iwstring.o _iwaray_iwstring.o \
	_iwaray_p_iwsubstring.o _iwarayq_iwstring.o \
	iw_stl_hash_map.o _iw_stl_hash_map_int.o _iw_stl_hash_map_uint.o \
	_iw_stl_hash_map_float.o _iw_stl_hash_map_iwstring.o \
	_iw_stl_hash_map_long.o _iw_stl_hash_map_long_long.o \
	_iw_stl_vector_string.o string_change.o remove_suffix.o iwzlib.o\
	_split_dv_int.o _split_dv_float.o _split_dv_double.o iwstring_strstream.o\
	string_relationals.o _iwarayq_substring.o hexprint.o \
	iwstring_and_file_descriptor.o iwstrncasecmp.o _append_int_form.o \
	append_number.o unhtml.o

SHARED_INC_FILES = iwstring.h iwhash.h iw_stl_hash_map.h iw_stl_hash_set.h \
	iw_stl_hash_multimap.h iw_stl_hash_multiset.h iw_stl_vector.h iwzlib.h 

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
TSTRING_EXE = $(BUILD_DIR)/test_iwstring
TSTRING_OBJECTS = $(BUILD_DIR)/tiwstring.o $(BUILD_DIR)/iwstring.o 
$(BUILD_DIR)/tiwstring.o: tiwstring.cc iwstring.h
	$(CXX) -c -o $@ $< $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES)
$(BUILD_DIR)/test_iwstring: $(TSTRING_OBJECTS)
	$(LD) -o $@ $(TSTRING_OBJECTS) $(LIB) -liwaray $(SYSTEM_LIBS)
#~ EXECUTABLES += $(TSTRING_EXE)

################# Template for public executables ######################
TSTL_EXE = $(BUILD_DIR)/test_stl
TSTL_OBJECTS = $(BUILD_DIR)/tstl.o $(BUILD_DIR)/iw_stl_hash_map.o
$(BUILD_DIR)/tstl.o: tstl.cc iwstring.h iw_stl_hash_map.h
	$(CXX) -c -o $@ $< $(FLAGS_FOR_TEMPLATE_INSTANTIATORS)
$(BUILD_DIR)/test_stl: $(TSTL_OBJECTS) library
	$(LD) -o $@ $(TSTL_OBJECTS) $(LIBS) -liwstring_data_source -liwcrex -lcmdline -liwstring -liwaray $(LIBZ)
#~ EXECUTABLES += $(TSTL_EXE)

################# Template for public executables ######################
TEST_IWZLIB_EXE = $(BUILD_DIR)/test_iwzlib
TEST_IWZLIB_OBJECTS = $(BUILD_DIR)/test_iwzlib.o $(BUILD_DIR)/iwzlib.o
$(BUILD_DIR)/test_iwzlib.o: test_iwzlib.cc iwzlib.h
	$(CXX) -c -o $@ $< $(FLAGS_FOR_TEMPLATE_INSTANTIATORS)
$(BUILD_DIR)/test_iwzlib: $(TEST_IWZLIB_OBJECTS)
	$(LD) -o $@ $(TEST_IWZLIB_OBJECTS) $(LIBS) $(LIBZ) -lcmdline -liwstring -liwaray 	
#~ EXECUTABLES += $(TEST_IWZLIB_EXE)

################# Template for public executables ######################
TEST_WRITE_SPEED_EXE = $(BUILD_DIR)/test_write_speed
TEST_WRITE_SPEED_OBJECTS = $(BUILD_DIR)/test_write_speed.o $(BUILD_DIR)/IWString_class.o $(BUILD_DIR)/iwwrite.o
$(BUILD_DIR)/test_write_speed.o: test_write_speed.cc
	$(CXX) -c -o $@ $< $(FLAGS_FOR_TEMPLATE_INSTANTIATORS)
$(BUILD_DIR)/test_write_speed: $(TEST_WRITE_SPEED_OBJECTS)
	$(LD) -o $@ $(TEST_WRITE_SPEED_OBJECTS) $(LIBS) -liwstring_data_source -liwcrex -lcmdline -liwstring -liwaray $(LIBZ) -lm
#~ EXECUTABLES += $(TEST_WRITE_SPEED_EXE)

################# Template for public executables ######################
WRITE_TEST_EXE = $(BUILD_DIR)/test_write
WRITE_TEST_OBJECTS = $(BUILD_DIR)/write_test.o $(BUILD_DIR)/IWString_class.o $(BUILD_DIR)/iwstring_and_file_descriptor.o
$(BUILD_DIR)/write_test.o: write_test.cc
	$(CXX) -c -o $@ $< $(FLAGS_FOR_TEMPLATE_INSTANTIATORS)
$(BUILD_DIR)/test_write: $(WRITE_TEST_OBJECTS)
	$(LD) -o $@ $(WRITE_TEST_OBJECTS) $(LIBS) -liwstring_data_source -liwcrex -lcmdline -liwstring -liwaray $(LIBZ) -lm
#~ EXECUTABLES += $(WRITE_TEST_EXE)

################# Template for public executables ######################
TSCLASS_EXE = $(BUILD_DIR)/test_sclass
TSCLASS_OBJECTS = $(BUILD_DIR)/tsclass.o $(BUILD_DIR)/should_match.o $(LIBRARY_OBJECTS)
$(BUILD_DIR)/tsclass.o: tsclass.cc
	$(CXX) -c $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES) $(basename $(notdir $@)).cc -o $@
$(BUILD_DIR)/should_match.o: should_match.cc
	$(CXX) -c $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES) $(basename $(notdir $@)).cc -o $@
$(BUILD_DIR)/test_sclass: $(TSCLASS_OBJECTS)
	$(LD) -o $@ $(TSCLASS_OBJECTS) $(LIBS) -lcmdline -liwstring -liwaray $(SYSTEM_LIBS) $(LIBZ)
#~ EXECUTABLES += $(TSCLASS_EXE)

################# Template for public executables ######################
TEST_HASH_FUNCTION_EXE = $(BUILD_DIR)/test_hash_function
TEST_HASH_FUNCTION_OBJECTS = $(BUILD_DIR)/test_hash_function.o $(LIBRARY_OBJECTS)
$(BUILD_DIR)/test_hash_function.o:
	$(CXX) -c $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES) $(basename $(notdir $@)).cc -o $@
$(BUILD_DIR)/test_hash_function: $(TEST_HASH_FUNCTION_OBJECTS)
	$(LD) -o $@ $(TEST_HASH_FUNCTION_OBJECTS) $(LIBS) -lcmdline -liwstring -liwaray -liwstring_data_source -liwmisc -liwcrex $(SYSTEM_LIBS) $(LIBZ)
#~ EXECUTABLES += $(TEST_HASH_FUNCTION_EXE)

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
#copy_exe: exe
#	$(CP) $(EXECUTABLES) $(EXEDIR)
	
clean:
	@$(RM) $(BUILD_DIR)/*.o || true

veryclean: clean
	@$(RM) $(BUILD_DIR)/*.a || true
#	@$(RM) $(EXECUTABLES) || true
