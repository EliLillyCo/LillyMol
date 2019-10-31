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
LIBS = $(LIB) -lgfp -liwmolecule -liwstring_data_source -liwtemplate -liwstring -lsparsefp -laccumulator -lcmdline  -liwhistogram -liwbits -liwaray -liwtdt -liwmisc -liwrandom -liwcrex  -lxmlParser  $(LIBZ) $(USRLIB) $(SYSTEM_LIBS)

LIBTEMPLATE = -liwtemplate
 
$(BUILD_DIR)/%.o : %.cc
	$(CXX) -c -o $@ $< $(CXX_FLAGS_WITH_IMPLICIT_TEMPLATES)
	
########################## Template for libary #########################
PROJECT_NAME = gfp

LIBRARY_OBJECTS_LIST = \
	dyfp.o gfp.o \
	sparsefp.o \
	various_distance_metrics.o \
	sparse_collection.o \
    fixed_size_counted_fingerprint_uchar.o \
    _fixed_size_counted_fingerprint_uchar.o \
    fixed_size_counted_fingerprint_uint.o \
    _fixed_size_counted_fingerprint_uint.o \
    bic_table.o \
	multi_conformer.o \
	tversky.o \
	gfp_contiguous.o \
	count_tdts.o \
	gfp_standard.o \
	smiles_id_dist.o neighbour_list.o fp_and_smiles.o \
	_resizable_array_sid.o _neighbour_list_fp_and_smiles.o \
	spread_v2.o _resizable_array_spread_v2_object.o

SHARED_INC_FILES = gfp.h tversky.h sparsefp.h dyfp.h leader.h spread_v2.h \
	multi_conformer.h gfp_standard.h fixed_size_counted_fingerprint.h \
	sparse_collection.h \
	smiles_id_dist.h neighbour_list.h
	
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
ADD_DESCRIPTORS_EXE = $(BUILD_DIR)/gfp_add_descriptors
ADD_DESCRIPTORS_OBJECTS = $(BUILD_DIR)/gfp_add_descriptors.o
$(BUILD_DIR)/gfp_add_descriptors: $(ADD_DESCRIPTORS_OBJECTS) 
	$(LD) $(PROF) -o $@ $(OPENMP) $(ADD_DESCRIPTORS_OBJECTS) $(LIBS) 
EXECUTABLES += $(ADD_DESCRIPTORS_EXE)

################## Template for public executables ######################
DESCRIPTORS_TO_01_FINGERPRINT_EXE = $(BUILD_DIR)/descriptor_file_to_01_fingerprints
DESCRIPTORS_TO_01_FINGERPRINT_OBJECTS = $(BUILD_DIR)/descriptor_file_to_01_fingerprints.o
$(BUILD_DIR)/descriptor_file_to_01_fingerprints: $(DESCRIPTORS_TO_01_FINGERPRINT_OBJECTS)
	$(LD) -o $@ $(PROF) $(DESCRIPTORS_TO_01_FINGERPRINT_OBJECTS) $(LIBS) -lsparsefp -liwmisc -liwbits $(SYSTEM_LIBS)	
EXECUTABLES += $(DESCRIPTORS_TO_01_FINGERPRINT_EXE)

################# Template for public executables ######################
DISTANCE_MATRIX_EXE = $(BUILD_DIR)/gfp_distance_matrix
DISTANCE_MATRIX_OBJECTS = $(BUILD_DIR)/gfp_distance_matrix.o $(BUILD_DIR)/build_pool.o
$(BUILD_DIR)/gfp_distance_matrix.o: gfp_distance_matrix.cc
$(BUILD_DIR)/gfp_distance_matrix: $(DISTANCE_MATRIX_OBJECTS)
	$(LD) $(PROF) -o $@ $(DISTANCE_MATRIX_OBJECTS) $(LIBS) -liwstring -liwaray
EXECUTABLES += $(DISTANCE_MATRIX_EXE)	

################## Template for public executables ######################
FINGERPRINTS_TO_DESCRIPTORS_EXE = $(BUILD_DIR)/gfp_to_descriptors
FINGERPRINTS_TO_DESCRIPTORS_OBJECTS = $(BUILD_DIR)/gfp_to_descriptors.o 
$(BUILD_DIR)/gfp_to_descriptors.o: gfp_to_descriptors.cc 
$(BUILD_DIR)/gfp_to_descriptors: $(FINGERPRINTS_TO_DESCRIPTORS_OBJECTS)
	$(LD) $(PROF) -o $@ $(FINGERPRINTS_TO_DESCRIPTORS_OBJECTS) $(LIBS)
EXECUTABLES += $(FINGERPRINTS_TO_DESCRIPTORS_EXE)

################# Template for public executables ######################
GFP_PAIRWISE_DISTANCES_EXE = $(BUILD_DIR)/gfp_pairwise_distances
GFP_PAIRWISE_DISTANCES_OBJECTS = $(BUILD_DIR)/gfp_pairwise_distances.o 
$(BUILD_DIR)/gfp_pairwise_distances.o: gfp_pairwise_distances.cc
$(BUILD_DIR)/gfp_pairwise_distances: $(GFP_PAIRWISE_DISTANCES_OBJECTS)
	$(LD) $(PROF) -o $@ $(GFP_PAIRWISE_DISTANCES_OBJECTS) $(LIBS) -liwtdt
EXECUTABLES += $(GFP_PAIRWISE_DISTANCES_EXE)

################# Template for public executables ######################
DISTANCE_FILTER_EXE = $(BUILD_DIR)/gfp_distance_filter
DISTANCE_FILTER_OBJECTS = $(BUILD_DIR)/gfp_distance_filter.o
$(BUILD_DIR)/gfp_distance_filter.o: gfp_distance_filter.cc
$(BUILD_DIR)/gfp_distance_filter: $(DISTANCE_FILTER_OBJECTS)
	$(LD) $(PROF) -o $@ $(DISTANCE_FILTER_OBJECTS) $(LIBS) -liwbits
EXECUTABLES += $(DISTANCE_FILTER_EXE)

################# Template for public executables ######################
SPARSE_TO_FIXED_EXE = $(BUILD_DIR)/gfp_sparse_to_fixed
SPARSE_TO_FIXED_OBJECTS = $(BUILD_DIR)/gfp_sparse_to_fixed.o 
$(BUILD_DIR)/gfp_sparse_to_fixed.o: gfp_sparse_to_fixed.cc
$(BUILD_DIR)/gfp_sparse_to_fixed: $(SPARSE_TO_FIXED_OBJECTS)
	$(LD) $(PROF) -o $@ $(SPARSE_TO_FIXED_OBJECTS) $(LIBS) -liwbits -liwtdt -laccumulator $(LIBTEMPLATE) $(SYSTEM_LIBS)
EXECUTABLES += $(SPARSE_TO_FIXED_EXE)

################# Template for public executables ######################
SINGLE_LINKAGE_EXE = $(BUILD_DIR)/gfp_single_linkage
SINGLE_LINKAGE_OBJECTS = $(BUILD_DIR)/gfp_single_linkage.o $(BUILD_DIR)/distance_cache.o 
$(BUILD_DIR)/gfp_single_linkage.o: gfp_single_linkage.cc 
$(BUILD_DIR)/gfp_single_linkage: $(SINGLE_LINKAGE_OBJECTS)
	$(LD) $(PROF) -o $@ $(SINGLE_LINKAGE_OBJECTS) -liwrandom $(LIBS)
EXECUTABLES += $(SINGLE_LINKAGE_EXE)


################# Template for public executables ######################
NEARNEIGHBOURS_SINGLE_FILE_EXE = $(BUILD_DIR)/gfp_nearneighbours_single_file
NEARNEIGHBOURS_SINGLE_FILE_OBJECTS = $(BUILD_DIR)/gfp_nearneighbours_single_file.o 
$(BUILD_DIR)/gfp_nearneighbours_single_file.o: gfp_nearneighbours_single_file.cc 
$(BUILD_DIR)/gfp_nearneighbours_single_file: $(NEARNEIGHBOURS_SINGLE_FILE_OBJECTS)
	$(LD) $(PROF) -o $@ $(NEARNEIGHBOURS_SINGLE_FILE_OBJECTS) $(LIBS) $(IWMALLOC)
EXECUTABLES += $(NEARNEIGHBOURS_SINGLE_FILE_EXE)	


################# Template for public executables ######################
NEARNEIGHBOURS_EXE = $(BUILD_DIR)/gfp_nearneighbours
NEARNEIGHBOURS_OBJECTS = $(BUILD_DIR)/gfp_nearneighbours.o 
$(BUILD_DIR)/gfp_nearneighbours.o: gfp_nearneighbours.cc
$(BUILD_DIR)/gfp_nearneighbours: $(NEARNEIGHBOURS_OBJECTS)
	$(LD) $(PROF) -o $@ $(NEARNEIGHBOURS_OBJECTS) $(LIBS)
EXECUTABLES += $(NEARNEIGHBOURS_EXE)

################# Template for public executables ######################
LEADER_V2_EXE = $(BUILD_DIR)/gfp_leader_v2
LEADER_V2_OBJECTS = $(BUILD_DIR)/gfp_leader_v2.o $(BUILD_DIR)/_rsa_leader.o
$(BUILD_DIR)/gfp_leader_v2.o: gfp_leader_v2.cc 
$(BUILD_DIR)/gfp_leader_v2: $(LEADER_V2_OBJECTS)
	$(LD) $(PROF) -o $@ $(LEADER_V2_OBJECTS) $(LIBS) -liwrandom 
EXECUTABLES += $(LEADER_V2_EXE)	

################# Template for public executables ######################
LARGE_NEARNEIGHBOURS_EXE = $(BUILD_DIR)/gfp_lnearneighbours
LARGE_NEARNEIGHBOURS_OBJECTS = $(BUILD_DIR)/gfp_lnearneighbours.o $(BUILD_DIR)/fb_bits_and_weights.o 
$(BUILD_DIR)/gfp_lnearneighbours: $(LARGE_NEARNEIGHBOURS_OBJECTS)
	$(LD) $(PROF) -o $@ $(LARGE_NEARNEIGHBOURS_OBJECTS) $(LIBS)
EXECUTABLES += $(LARGE_NEARNEIGHBOURS_EXE)	

################# Template for public executables ######################
PROFILE_ACTIVITY_EXE = $(BUILD_DIR)/gfp_profile_activity_by_bits
PROFILE_ACTIVITY_OBJECTS = $(BUILD_DIR)/gfp_profile_activity_by_bits.o $(BUILD_DIR)/bit_and_weight.o 
$(BUILD_DIR)/gfp_profile_activity_by_bits.o: gfp_profile_activity_by_bits.cc nbc.h
$(BUILD_DIR)/gfp_profile_activity_by_bits: $(PROFILE_ACTIVITY_OBJECTS)
	$(LD) $(PROF) -o $@ $(PROFILE_ACTIVITY_OBJECTS) $(LIBS) -liwbits   -liwcrex -laccumulator $(LIBTEMPLATE) $(SYSTEM_LIBS)
EXECUTABLES += $(PROFILE_ACTIVITY_EXE)

################# Template for public executables ######################
SPREAD_V2_EXE = $(BUILD_DIR)/gfp_spread_v2
SPREAD_V2_OBJECTS = $(BUILD_DIR)/gfp_spread_v2.o \
	$(BUILD_DIR)/_resizable_array_spread_v2_object.o 
$(BUILD_DIR)/gfp_spread_v2.o: gfp_spread_v2.cc 
$(BUILD_DIR)/gfp_spread_v2: $(SPREAD_V2_OBJECTS)
	$(LD) $(PROF) -o $@ $(SPREAD_V2_OBJECTS) -liwrandom $(LIBS) -liwstring -liwaray -lm
EXECUTABLES += $(SPREAD_V2_EXE)

################# Template for public executables ######################
PLOTNN_EXE = $(BUILD_DIR)/nplotnn
PLOTNN_OBJECTS = $(BUILD_DIR)/nplotnn.o $(BUILD_DIR)/distance_scaling.o
$(BUILD_DIR)/nplotnn.o: nplotnn.cc distance_scaling.h
$(BUILD_DIR)/nplotnn: $(PLOTNN_OBJECTS)
	$(LD) $(PROF) -o $@ $(PLOTNN_OBJECTS) $(LIBS)
EXECUTABLES += $(PLOTNN_EXE)

################# Template for public executables ######################
SPREAD_BUCKETS_V2_EXE = $(BUILD_DIR)/gfp_spread_buckets_v2
SPREAD_BUCKETS_V2_OBJECTS = $(BUILD_DIR)/gfp_spread_buckets_v2.o 
$(BUILD_DIR)/gfp_spread_buckets_v2.o: gfp_spread_buckets_v2.cc
$(BUILD_DIR)/gfp_spread_buckets_v2: $(SPREAD_BUCKETS_V2_OBJECTS)
	$(LD) $(PROF) -o $@ $(SPREAD_BUCKETS_V2_OBJECTS) -liwrandom -liwmisc $(LIBS) -liwbits
EXECUTABLES += $(SPREAD_BUCKETS_V2_EXE)

################# Template for public executables ######################
TDT_SORT_EXE = $(BUILD_DIR)/tdt_sort
TDT_SORT_OBJECTS = $(BUILD_DIR)/tdt_sort.o
$(BUILD_DIR)/tdt_sort: $(TDT_SORT_OBJECTS)
	$(LD) $(PROF) -o $@ $(TDT_SORT_OBJECTS) -lcmdline_v2 -liwrandom -liwmisc $(LIBS)
EXECUTABLES += $(TDT_SORT_EXE)

################# Template for public executables ######################
TDT_JOIN_EXE = $(BUILD_DIR)/tdt_join
TDT_JOIN_OBJECTS = $(BUILD_DIR)/tdt_join.o
$(BUILD_DIR)/tdt_join: $(TDT_JOIN_OBJECTS)
	$(LD) -o $@ $(TDT_JOIN_OBJECTS) $(LIBS) -liwmisc $(SYSTEM_LIBS)
EXECUTABLES += $(TDT_JOIN_EXE)

################## Template for public executables ######################
#DISTANCE_MATRIX_DMC_EXE = $(BUILD_DIR)/gfp_distance_matrix_iwdm
#DISTANCE_MATRIX_DMC_OBJECTS = $(BUILD_DIR)/gfp_distance_matrix_iwdm.o 
#$(BUILD_DIR)/gfp_distance_matrix_iwdm.o: gfp_distance_matrix_iwdm.cc
#$(BUILD_DIR)/gfp_distance_matrix_iwdm: $(DISTANCE_MATRIX_DMC_OBJECTS)
#	$(LD) $(PROF) -o $@ $(DISTANCE_MATRIX_DMC_OBJECTS) -liwdistance_matrix $(LIBS)
#EXECUTABLES += $(DISTANCE_MATRIX_DMC_EXE)	

################# Template for public executables ######################
LEADER_STD_EXE = $(BUILD_DIR)/gfp_leader_standard
LEADER_STD_OBJECTS = $(BUILD_DIR)/gfp_leader_standard.o 
$(BUILD_DIR)/gfp_leader_standard.o: gfp_leader_standard.cc 
$(BUILD_DIR)/gfp_leader_standard: $(LEADER_STD_OBJECTS)
	$(LD) $(PROF) -o $@ $(OPENMP) $(LEADER_STD_OBJECTS) $(LIBS) $(SYSTEM_LIBS) -liwmisc
EXECUTABLES += $(LEADER_STD_EXE)	

################# Template for public executables ######################
LARGE_NEARNEIGHBOURS_STANDARD_EXE = $(BUILD_DIR)/gfp_lnearneighbours_standard
LARGE_NEARNEIGHBOURS_STANDARD_OBJECTS = $(BUILD_DIR)/gfp_lnearneighbours_standard.o
$(BUILD_DIR)/gfp_lnearneighbours_standard.o: gfp_lnearneighbours_standard.cc
$(BUILD_DIR)/gfp_lnearneighbours_standard: $(LARGE_NEARNEIGHBOURS_STANDARD_OBJECTS)
	$(LD) $(PROF) -o $@ $(LARGE_NEARNEIGHBOURS_STANDARD_OBJECTS) $(LIBS)
EXECUTABLES += $(LARGE_NEARNEIGHBOURS_STANDARD_EXE)	

################## Template for public executables ######################
GFP_NAIVE_BAYESIAN_EXE = $(BUILD_DIR)/gfp_naive_bayesian
GFP_NAIVE_BAYESIAN_OBJECTS = $(BUILD_DIR)/gfp_naive_bayesian.o 
$(BUILD_DIR)/gfp_naive_bayesian.o: gfp_naive_bayesian.cc 
#bit_subset_v2.h
$(BUILD_DIR)/gfp_naive_bayesian: $(GFP_NAIVE_BAYESIAN_OBJECTS)
	$(LD) $(PROF) -o $@ $(GFP_NAIVE_BAYESIAN_OBJECTS) $(LIBS) -liwbits   -liwcrex -laccumulator $(LIBTEMPLATE) $(SYSTEM_LIBS)
EXECUTABLES += $(GFP_NAIVE_BAYESIAN_EXE)

################# Template for public executables ######################
SPREAD_STANDARD_EXE = $(BUILD_DIR)/gfp_spread_standard
SPREAD_STANDARD_OBJECTS = $(BUILD_DIR)/gfp_spread_standard.o 
$(BUILD_DIR)/gfp_spread_standard: $(SPREAD_STANDARD_OBJECTS) 
	$(LD) $(PROF) -o $@ $(OPENMP) $(SPREAD_STANDARD_OBJECTS) -liwrandom -liwbits -lm $(LIBS) 
EXECUTABLES += $(SPREAD_STANDARD_EXE)

################# Template for public executables ######################
FINGERPRINTS_TO_DESCRIPTORS_MULTIPLE_EXE = $(BUILD_DIR)/gfp_to_descriptors_multiple
FINGERPRINTS_TO_DESCRIPTORS_MULTIPLE_OBJECTS = $(BUILD_DIR)/gfp_to_descriptors_multiple.o
$(BUILD_DIR)/gfp_to_descriptors_multiple.o: gfp_to_descriptors_multiple.cc 
$(BUILD_DIR)/gfp_to_descriptors_multiple: $(FINGERPRINTS_TO_DESCRIPTORS_MULTIPLE_OBJECTS)
	$(LD) $(PROF) -o $@ $(FINGERPRINTS_TO_DESCRIPTORS_MULTIPLE_OBJECTS) $(LIBS)
EXECUTABLES += $(FINGERPRINTS_TO_DESCRIPTORS_MULTIPLE_EXE)

################# Template for public executables ######################
LEADER_MULTIPLE_THRESHOLD_EXE = $(BUILD_DIR)/nn_leader_and_jp
LEADER_MULTIPLE_THRESHOLD_OBJECTS = $(BUILD_DIR)/nn_leader_and_jp.o 
$(BUILD_DIR)/nn_leader_and_jp.o: nn_leader_and_jp.cc
$(BUILD_DIR)/nn_leader_and_jp: $(LEADER_MULTIPLE_THRESHOLD_OBJECTS)
	$(LD) $(PROF) -o $@ $(LEADER_MULTIPLE_THRESHOLD_OBJECTS) $(LIBS) -laccumulator -liwmisc
EXECUTABLES += $(LEADER_MULTIPLE_THRESHOLD_EXE)


################# Template for public executables ######################
DISTANCE_MATRIX_DMC_EXE = $(BUILD_DIR)/gfp_distance_matrix_iwdm
DISTANCE_MATRIX_DMC_OBJECTS = $(BUILD_DIR)/gfp_distance_matrix_iwdm.o 
$(BUILD_DIR)/gfp_distance_matrix_iwdm.o: gfp_distance_matrix_iwdm.cc
$(BUILD_DIR)/gfp_distance_matrix_iwdm: $(DISTANCE_MATRIX_DMC_OBJECTS)
	$(LD) $(PROF) -o $@ $(DISTANCE_MATRIX_DMC_OBJECTS) -liwdistance_matrix $(LIBS)
EXECUTABLES += $(DISTANCE_MATRIX_DMC_EXE)

################# Template for public executables ######################
IWSTATS_EXE = $(BUILD_DIR)/iwstats
IWSTATS_OBJECTS = $(BUILD_DIR)/iwstats.o $(BUILD_DIR)/bsquared.o $(BUILD_DIR)/Metric.o $(BUILD_DIR)/Enrichment.o
$(BUILD_DIR)/iwstats.o: iwstats.cc bsquared.h
$(BUILD_DIR)/iwstats: $(IWSTATS_OBJECTS)
	$(LD) $(PROF) -o $@ $(IWSTATS_OBJECTS) -liwmisc $(LIBS) -liwrandom 
EXECUTABLES += $(IWSTATS_EXE)
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
exe: $(EXECUTABLES)

# Rule for install executables
copy_exe: exe
	$(CP) $(EXECUTABLES) $(EXEDIR)
	
clean:
	@$(RM) $(BUILD_DIR)/*.o || true

veryclean: clean
	@$(RM) $(BUILD_DIR)/*.a || true
	@$(RM) $(EXECUTABLES) || true
