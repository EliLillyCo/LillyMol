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
LIBS = $(LIB) -lcmdline -liwrandom -liwstring -liwaray

SYSTEM_LIBS = -lm

########################## Template for libary #########################
# Change the project name
PROJECT_NAME = iwqsort

LIBRARY_OBJECTS_LIST = _iwqsort_int.o

SHARED_INC_FILES = iwqsort.h iwqsort_mbrfn.h iwqsort_tbb.h iwqsort_ptr_fo_.h	

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
SORT_FILE_EXE = $(BUILD_DIR)/test_sort_file
SORT_FILE_OBJECTS = $(BUILD_DIR)/sort_file.o 
$(BUILD_DIR)/sort_file.o: sort_file.cc iwqsort.h iwqsort_.h
	$(CXX) -c -o $@ $< $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES)
$(BUILD_DIR)/test_sort_file: $(SORT_FILE_OBJECTS)
	$(LD) -o $@ $(PROF) $(SORT_FILE_OBJECTS) $(LIBS) -liwstring_data_source -liwcrex -liwstring -lcmdline $(LIBTEMPLATE) $(SYSTEM_LIBS) $(LIBZ)
#~ EXECUTABLES += $(SORT_FILE_EXE)

################# Template for public executables ######################
TEST_IWQSORT_EXE = $(BUILD_DIR)/test_iwqsort 
TEST_IWQSORT_OBJECTS = $(BUILD_DIR)/test_iwqsort.o $(BUILD_DIR)/_foo_int.o $(BUILD_DIR)/_foo_float.o
$(BUILD_DIR)/test_iwqsort.o: test_iwqsort.cc iwqsort.h iwqsort_.h foo.h foo_.h
	$(CXX) -c -o $@ $< $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES)
$(BUILD_DIR)/_foo_int.o: _foo_int.cc foo.h foo_.h
	$(CXX) -c -o $@ $< $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES)
$(BUILD_DIR)/_foo_float.o: _foo_float.cc foo.h foo_.h
	$(CXX) -c -o $@ $< $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES)
$(BUILD_DIR)/test_iwqsort: $(TEST_IWQSORT_OBJECTS)
	$(LD) -o $@ $(PROF) $(TEST_IWQSORT_OBJECTS) $(LIBS) $(LIBTEMPLATE) $(SYSTEM_LIBS)
#~ EXECUTABLES += $(TEST_IWQSORT_EXE)

################# Template for public executables ######################
TEST_QSORT_EXE = $(BUILD_DIR)/test_qsort 
TEST_QSORT_OBJECTS = $(BUILD_DIR)/test_qsort.o $(BUILD_DIR)/_foo_int.o $(BUILD_DIR)/_foo_float.o
# See above for _foo_float.o and _foo_int.o build rule
$(BUILD_DIR)/test_qsort.o: test_qsort.cc foo.h foo_.h
	$(CXX) -c -o $@ $< $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES)
$(BUILD_DIR)/test_qsort: $(TEST_QSORT_OBJECTS)
	$(LD) -o $@ $(PROF) $(TEST_QSORT_OBJECTS) $(LIBS) $(LIBTEMPLATE) $(SYSTEM_LIBS)
#~ EXECUTABLES += $(TEST_QSORT_EXE)

################# Template for public executables ######################
TEST_IWQSORT_PTR_EXE = $(BUILD_DIR)/test_iwqsort_ptr 
TEST_IWQSORT_PTR_OBJECTS = $(BUILD_DIR)/test_iwqsort_ptr.o $(BUILD_DIR)/_foo_int.o $(BUILD_DIR)/_foo_float.o
# See above for rules for _foo_init and foo_float 
$(BUILD_DIR)/test_iwqsort_ptr.o: test_iwqsort_ptr.cc iwqsort_ptr.h iwqsort_ptr_.h foo.h foo_.h
	$(CXX) -c -o $@ $< $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES)
$(BUILD_DIR)/test_iwqsort_ptr: $(TEST_IWQSORT_PTR_OBJECTS)
	$(LD) -o $@ $(PROF) $(TEST_IWQSORT_PTR_OBJECTS) $(LIBS) $(LIBTEMPLATE) $(SYSTEM_LIBS)
#~ EXECUTABLES += $(TEST_IWQSORT_PTR_EXE)

################# Template for public executables ######################
TEST_IWQSORT_MFN_EXE = $(BUILD_DIR)/test_iwqsort_mfn 
TEST_IWQSORT_MFN_OBJECTS = $(BUILD_DIR)/test_iwqsort_mfn.o $(BUILD_DIR)/_foo_int.o $(BUILD_DIR)/_foo_float.o
# See above for rules for _foo_init and foo_float 
$(BUILD_DIR)/test_iwqsort_mfn.o: test_iwqsort_mfn.cc iwqsort_mfn.h iwqsort_mfn_.h foo.h foo_.h
	$(CXX) -c -o $@ $< $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES)
$(BUILD_DIR)/test_iwqsort_mfn: $(TEST_IWQSORT_MFN_OBJECTS)
	$(LD) -o $@ $(PROF) $(TEST_IWQSORT_MFN_OBJECTS) $(LIBS) $(LIBTEMPLATE) $(SYSTEM_LIBS) $(LIBZ)
#~ EXECUTABLES += $(TEST_IWQSORT_MFN_EXE)

################# Template for public executables ######################
TEST_IWQSORT_PMFN_EXE = $(BUILD_DIR)/test_iwqsort_pmfn
TEST_IWQSORT_PMFN_OBJECTS = $(BUILD_DIR)/test_iwqsort_pmfn.o $(BUILD_DIR)/_foo_int.o $(BUILD_DIR)/_foo_float.o
# See above for rules for _foo_init and foo_float 
$(BUILD_DIR)/test_iwqsort_pmfn.o: test_iwqsort_pmfn.cc iwqsort_pmfn.h iwqsort_pmfn_.h foo.h foo_.h
	$(CXX) -c -o $@ $< $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES)
$(BUILD_DIR)/test_iwqsort_pmfn: $(TEST_IWQSORT_PMFN_OBJECTS)
	$(LD) -o $@ $(PROF) $(TEST_IWQSORT_PMFN_OBJECTS) $(LIBS) $(LIBTEMPLATE) $(SYSTEM_LIBS)
#~ EXECUTABLES += $(TEST_IWQSORT_PMFN_EXE)

################# Template for public executables ######################
TEST_IWQSORT_FO_EXE = $(BUILD_DIR)/test_iwqsort_fo 
TEST_IWQSORT_FO_OBJECTS = $(BUILD_DIR)/test_iwqsort_fo.o
# See above for rules for _foo_init and foo_float 
$(BUILD_DIR)/test_iwqsort_fo.o: test_iwqsort_fo.cc iwqsort.h
	$(CXX) -c -o $@ $< $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES)
$(BUILD_DIR)/test_iwqsort_fo: $(TEST_IWQSORT_FO_OBJECTS)
	$(LD) -o $@ $(PROF) $(TEST_IWQSORT_FO_OBJECTS) $(LIBS) $(LIBTEMPLATE) $(SYSTEM_LIBS)
#~ EXECUTABLES += $(TEST_IWQSORT_FO_EXE)

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
#	@$(RM) $(EXECUTABLES) || more
