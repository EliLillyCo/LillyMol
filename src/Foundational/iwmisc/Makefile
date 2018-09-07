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
PROJECT_NAME = iwmisc

LIBRARY_C_OBJECTS_LIST = md5c.o

LIBRARY_OBJECTS_LIST = new_int.o\
	new_unsigned_int.o\
	new_float.o\
	new_double.o\
	write_space_suppressed_string.o\
	write_first_token_of_string.o\
	append_first_token_of_name.o\
	iwdigits.o\
	file_find_dividing_points.o\
	fraction_as_string.o\
	rhswab.o\
	endian.o\
	logical_expression.o\
	hton_stuff.o\
	msi_object.o\
	int_comparator.o\
	normalisation.o\
	iw_tabular_data.o\
	dash_s.o\
	dash_f.o\
	dash_d.o\
	dash_x.o\
	iwuuencode.o iwuudecode.o\
	fetch_via_regexp.o\
	_set_vector_int.o\
	_set_vector_unsigned_int.o\
	_set_vector_float.o\
	_set_vector_double.o\
	_set_vector_uchar.o\
	_count_occurrences_int.o\
	_count_occurrences_float.o\
	_copy_vector_int.o\
	_copy_vector_uint.o\
	_copy_vector_float.o\
	_copy_vector_double.o\
	_sum_vec_double.o\
	_sum_vec_float.o\
	_sum_vec_unsigned_int.o\
	_sum_vec_int.o\
	_iwswap_int.o\
	_locate_item_int.o\
	_iwfactorial_int.o\
	_iwfactorial_uint.o\
	_iwfactorial_uint64.o\
	_iwfactorial_double.o\
	_count_non_zero_occurrences_int.o\
	_stl_stuff.o\
	_stl_stuff_unsigned_int.o\
	_stl_fill_int.o\
	_stl_fill_unsigned_int.o\
	_numeric_data_from_file_float.o\
	_numeric_data_from_file_int.o\
	_iw_max_of_array_int.o\
	_iw_max_of_array_double.o\
	_iw_write_array_int.o\
	_iw_write_array_float.o\
	_iw_write_array_double.o\
	_resizable_array_msi.o\
	_report_progress_uint.o\
	_iwarchive_int.o\
	_iwarchive_uint.o\
	_iwarchive_float.o\
	_iwarchive_double.o\
	_set_or_unset_int.o\
	_set_or_unset_float.o\
	_set_or_unset_long.o\
	_set_or_unset_long_long.o\
	_set_or_unset_streampos.o\
	_set_or_unset_double.o\
	_set_or_unset_iwstring.o\
	_minmaxspc_double.o\
	_minmaxspc_float.o\
	_minmaxspc_int.o\
	_resizable_array_mmspc_int.o\
	_iwminmax_int.o\
	_iwmin_int.o\
	_iwmin_float.o\
	_iwmin_double.o\
	_iwmax_int.o\
	_iwmax_float.o\
	_iwmax_double.o\
	_iwminmax_float.o\
	_iwminmax_double.o\
	_iwminid_int_int.o\
	_iwmaxid_int_int.o\
	_iwminid_double_int.o\
	_iwmaxid_double_int.o\
	_iwminid_double_iwstring.o\
	_iwmaxid_double_iwstring.o\
	_iwmaxid_float_iwstring.o

SHARED_INC_FILES = iwconfig.h misc.h iw_time.h primes.h sparse_fp_creator.h \
	numeric_data_from_file.h iw_auto_array.h msi_object.h logical_expression.h \
	iwdigits.h report_progress.h iwarchive.h set_or_unset.h normalisation.h \
	iwminmax.h minmaxspc.h fetch_via_regexp.h iw_tabular_data.h md5.h md5_global.h

LIBRARY_OBJECTS = $(addprefix $(BUILD_DIR)/, $(LIBRARY_OBJECTS_LIST))
LIBRARY_C_OBJECTS = $(addprefix $(BUILD_DIR)/, $(LIBRARY_C_OBJECTS_LIST))


LIBRARY_STEM = lib$(PROJECT_NAME)

LIBRARY_NAME = $(BUILD_DIR)/$(LIBRARY_STEM).a

# Rule for making library obj
$(LIBRARY_OBJECTS):
	$(CXX) -c $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES) $(basename $(notdir $@)).cc -o $@

$(LIBRARY_C_OBJECTS):
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $(basename $(notdir $@)).c -o $@ $<
	 
# Rule for making library
$(LIBRARY_NAME): $(LIBRARY_OBJECTS) $(LIBRARY_C_OBJECTS)
	$(AR) r $@ $(LIBRARY_OBJECTS) $(LIBRARY_C_OBJECTS)
	
########################## Template for libary #########################
# Change the project name
PROJECT_NAME_2 = sparsefp

LIBRARY_OBJECTS_LIST_2 = sparse_fp_creator.o sparse_fp_creator_support.o _unordered_map_to_md5.o

LIBRARY_OBJECTS_2 = $(addprefix $(BUILD_DIR)/, $(LIBRARY_OBJECTS_LIST_2))

LIBRARY_STEM_2 = lib$(PROJECT_NAME_2)

LIBRARY_NAME_2 = $(BUILD_DIR)/$(LIBRARY_STEM_2).a

# Rule for making library obj
$(LIBRARY_OBJECTS_2):
	$(CXX) -c $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES) $(basename $(notdir $@)).cc -o $@

# Rule for making library
$(LIBRARY_NAME_2): $(LIBRARY_OBJECTS_2)
	$(AR) r $@ $(LIBRARY_OBJECTS_2)

################# Template for public executables ######################
TEST_SWAB_EXE = $(BUILD_DIR)/test_swab
TEST_SWAB_OBJECTS = $(BUILD_DIR)/test_swab.o $(BUILD_DIR)/rhswab.o
$(BUILD_DIR)/test_swab.o: test_swab.cc
	$(CXX) -c -o $@ $< $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES)
$(BUILD_DIR)/test_swab: $(TEST_SWAB_OBJECTS)
	$(LD) -o $@ $(TEST_SWAB_OBJECTS) $(LIB) -liwrandom -lcmdline -liwstring -liwaray
#~ EXECUTABLES += $(TEST_SWAB_EXE)

################# Template for public executables ######################
TEST_IWUUENCODE_EXE = $(BUILD_DIR)/test_iwuuencode
TEST_IWUUENCODE_OBJECTS = $(BUILD_DIR)/test_iwuuencode.o $(BUILD_DIR)/iwuuencode.o
$(BUILD_DIR)/test_iwuuencode.o: test_iwuuencode.cc misc.h
	$(CXX) -c -o $@ $< $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES)
$(BUILD_DIR)/test_iwuuencode: $(TEST_IWUUENCODE_OBJECTS)
	$(LD) -o $@ $(TEST_IWUUENCODE_OBJECTS) $(LIB) -liwstring_data_source -liwcrex -lcmdline -liwstring -liwaray
#~ EXECUTABLES += $(TEST_IWUUENCODE_EXE)

################# Template for public executables ######################
TEST_IWUUDECODE_EXE = $(BUILD_DIR)/test_iwuudecode
TEST_IWUUDECODE_OBJECTS = $(BUILD_DIR)/test_iwuudecode.o $(BUILD_DIR)/iwuudecode.o
$(BUILD_DIR)/test_iwuudecode.o: test_iwuudecode.cc misc.h
	$(CXX) -c -o $@ $< $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES)
$(BUILD_DIR)/test_iwuudecode: $(TEST_IWUUDECODE_OBJECTS)
	$(LD) -o $@ $(TEST_IWUUDECODE_OBJECTS) $(LIB) -liwstring_data_source -liwcrex -lcmdline -liwstring -liwaray $(LIBZ)
#~ EXECUTABLES += $(TEST_IWUUDECODE_EXE)

################# Template for public executables ######################
TEST_IWDIGITS_EXE = $(BUILD_DIR)/test_iwdigits
TEST_IWDIGITS_OBJECTS = $(BUILD_DIR)/test_iwdigits.o $(BUILD_DIR)/iwdigits.o $(BUILD_DIR)/fraction_as_string.o $(BUILD_DIR)/dash_d.o $(BUILD_DIR)/dash_s.o $(BUILD_DIR)/dash_x.o
$(BUILD_DIR)/test_iwdigits.o: test_iwdigits.cc iwdigits.h
	$(CXX) -c -o $@ $< $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES)
$(BUILD_DIR)/test_iwdigits: $(TEST_IWDIGITS_OBJECTS)
	$(LD) -o $@ $(TEST_IWDIGITS_OBJECTS) $(LIB) -lcmdline_v2 -liwmisc -liwstring -liwaray $(LIBZ)
#~ EXECUTABLES += $(TEST_IWDIGITS_EXE)

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
library: $(LIBRARY_NAME) $(LIBRARY_NAME_2)

# Rule for copying static library
copy_library: copy_include library
	$(CP) $(LIBRARY_NAME) $(LIBDIR)
	$(CP) $(LIBRARY_NAME_2) $(LIBDIR)
	
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


