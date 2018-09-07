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

IWPROGRAMMES ?= $(PWD)/../..
UNAME ?= Linux-gcc-6.2.0

prefix = $(IWPROGRAMMES)
exec_prefix = $(prefix)

libdir = ${exec_prefix}/lib/$(UNAME)
includedir = ${prefix}/include
bindir = ${exec_prefix}/bin/$(UNAME)

CP=cp -p

SUBDIRS := $(wildcard */.)

.PHONY : $(SUBDIRS)
.PHONY : $(DO_X)
.PHONY: copy_include
.PHONY: library
.PHONY: copy_library
.PHONY: exe
.PHONY: copy_exe
.PHONY: clean
.PHONY: veryclean


DO_X = do-clean\
	do-veryclean\
	do-copy_include\
	do-library\
	do-copy_library\
	do-exe\
	do-copy_exe

# Enter each directory to build
# Build to Build directory
#~ 			echo "Enter sub-makefile" ; \
#~ 			if [ "$${target}" = "veryclean" ] ; \
#~ 			then \
#~ 				echo "Removing $(UNAME)" ; \
#~ 				if [ -d ./$(UNAME) ] ; \
#~ 				then \
#~ 					rm -r ./$(UNAME) ; \
#~ 				fi ; \			
#~ 			else \
#~ 				if [ ! -d ./$(UNAME) ] ; \
#~ 				then \
#~ 					mkdir ./$(UNAME) ; \
#~ 				fi ; \
#~ 			fi ; \

$(DO_X):
	@target=`echo $@ | sed -e 's/^do-//'`; \
	for dir in $(SUBDIRS) ; do \
		cd $${dir} ; \
		if [ -f ./Makefile ]; \
		then \
			if [ "$${target}" != "veryclean" ] ; \
			then \
				if [ ! -d ./$(UNAME) ] ; \
				then \
					mkdir ./$(UNAME) ; \
				fi ; \
			fi; \
			IWPROGRAMMES=$(IWPROGRAMMES) UNAME=$(UNAME) BUILD_DIR=$(UNAME) CXX=$(CXX) $(MAKE) -f Makefile $${target} ; \
			if [ "$${target}" = "veryclean" ] ; \
			then \
				if [ -d ./$(UNAME) ] ; \
				then \
					rm -r ./$(UNAME) ; \
				fi ; \
			fi; \
		fi; \
		cd ..;\
	done;

include: includedir copy_include_files

includedir:
	echo "Making ${includedir}"
	if [ ! -d ${includedir} ] ; then mkdir ${includedir} ; fi
	cp -p include/* $(includedir)

libdir: $(libdir)

bindir: $(bindir)

copy_include: includedir do-copy_include
	
veryclean: do-veryclean

clean: do-clean

exe: do-exe

copy_exe: do-copy_exe

library: do-library

copy_library: do-copy_library
