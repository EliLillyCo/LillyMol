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
LIBS = $(LIB) -liwhistogram -lcmdline -liwrandom -laccumulator -liwmisc\
	-liwstring_data_source -liwcrex -liwstring -liwaray -liwbits

SYSTEM_LIBS = -lm $(LIBZ)

########################## Template for libary #########################
# Change the project name
PROJECT_NAME = iwhistogram

LIBRARY_OBJECTS_LIST = iwhistogram.o resizable_histogram.o fuzzy_histogram.o histogram_fp.o overlapping_buckets.o \
	 _most_recent_int.o _most_recent_float.o _most_recent_double.o \
	 _running_average_int.o _running_average_float.o _running_average_double.o

SHARED_INC_FILES = iwhistogram.h 

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
TDIST_EXE = $(BUILD_DIR)/test_dist
TDIST_OBJECTS = $(BUILD_DIR)/tdist.o $(BUILD_DIR)/iwhistogram.o $(BUILD_DIR)/resizable_histogram.o $(BUILD_DIR)/fuzzy_histogram.o
$(BUILD_DIR)/tdist.o: tdist.cc
	$(CXX) -c -o $@ $< $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES)
$(BUILD_DIR)/test_dist: $(TDIST_OBJECTS)
	$(LD) -o $@ $(PROF) $(TDIST_OBJECTS) $(LIBS) $(SYSTEM_LIBS)
#~ EXECUTABLES += $(TDIST_EXE)

################# Template for public executables ######################
TMR_EXE = $(BUILD_DIR)/test_mr
TMR_OBJECTS = $(BUILD_DIR)/tmr.o $(BUILD_DIR)/_most_recent_int.o $(BUILD_DIR)/_most_recent_float.o $(BUILD_DIR)/_most_recent_double.o
$(BUILD_DIR)/tmr.o: tmr.cc
	$(CXX) -c -o $@ $< $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES)
$(BUILD_DIR)/test_mr: $(TMR_OBJECTS)
	$(LD) -o $@ $(PROF) $(TMR_OBJECTS) $(LIBS) $(SYSTEM_LIBS)
#~ EXECUTABLES += $(TMR_EXE)

################# Template for public executables ######################
TRA_EXE = $(BUILD_DIR)/test_ra
TRA_OBJECTS = $(BUILD_DIR)/tra.o $(BUILD_DIR)/_most_recent_int.o $(BUILD_DIR)/_most_recent_float.o $(BUILD_DIR)/_most_recent_double.o
$(BUILD_DIR)/tra.o: tra.cc
	$(CXX) -c -o $@ $< $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES)
$(BUILD_DIR)/test_ra: $(TRA_OBJECTS)
	$(LD) -o $@ $(PROF) $(TRA_OBJECTS) $(LIBS) $(SYSTEM_LIBS)
#~ EXECUTABLES += $(TRA_EXE)

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

# Rule for install executables
#copy_exe:
#	$(CP) $(EXECUTABLES) $(EXEDIR)
	
clean:
	@$(RM) $(BUILD_DIR)/*.o || true

veryclean: clean
	@$(RM) $(BUILD_DIR)/*.a || true
	@$(RM) $(EXECUTABLES) || true
