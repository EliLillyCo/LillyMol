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
IWPROGRAMMES ?= $(PWD)
UNAME ?= Linux-gcc-8.3.0
#Linux-gcc-6.2.0

prefix = $(IWPROGRAMMES)
exec_prefix = $(prefix)

LIB_ROOT_DIR = ${exec_prefix}/lib
INC_ROOT_DIR =${prefix}/include
BIN_ROOT_DIR = ${exec_prefix}/bin
libdir = $(LIB_ROOT_DIR)/$(UNAME)
includedir = $(INC_ROOT_DIR)
bindir = $(BIN_ROOT_DIR)/$(UNAME)

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
$(DO_X):
	@target=`echo $@ | sed -e 's/^do-//'`; \
	if [ "$${target}" != "veryclean" ] ; \
	then \
	  if [ ! -d ${libdir} ] ; \
	  then \
	    mkdir -p ${libdir} ; \
	  fi ; \
	  if [ ! -d ${bindir} ] ; \
	  then \
	    mkdir -p ${bindir} ; \
	  fi ; \
	  if [ ! -d ${includedir} ] ; \
	  then \
	    mkdir -p ${includedir} ; \
	  fi ; \
	fi ; \
	for dir in $(SUBDIRS) ; \
	do \
		cd $${dir};\
		if [ -f ./Makefile ]; \
		then \
			if [ "$${target}" = "veryclean" ] ; \
			then \
				if [ -d ./$(UNAME) ] ; \
				then \
					rm -r ./$(UNAME) ; \
				fi ; \
			else \
				if [ ! -d ./$(UNAME) ] ; \
				then \
					mkdir ./$(UNAME) ; \
				fi ; \
			fi; \
			IWPROGRAMMES=$(IWPROGRAMMES) UNAME=$(UNAME) BUILD_DIR=$(UNAME) CXX=$(CXX) $(MAKE) -j 5 -f Makefile $${target} ; \
		fi; \
		cd ..;\
	done; \
	if [ "$${target}" = "veryclean" ] ; \
	then \
	  if [ -d ${LIB_ROOT_DIR} ] ; \
	  then \
	    rm -r ${LIB_ROOT_DIR} ; \
	  fi ; \
	  if [ -d ${BIN_ROOT_DIR} ] ; \
	  then \
	    rm -r ${BIN_ROOT_DIR} ; \
	  fi ; \
	  if [ -d ${INC_ROOT_DIR} ] ; \
	  then \
	    rm -r ${INC_ROOT_DIR} ; \
	  fi ; \
	fi;

include: includedir copy_include_files

includedir:
	echo "Making ${includedir}"
	if [ ! -d ${includedir} ] ; then mkdir ${includedir} ; fi
	cp -p include/* $(includedir)

libdir: $(libdir)

bindir: $(bindir)

copy_include: do-copy_include

veryclean: do-veryclean

clean: do-clean

exe: do-exe

copy_exe: do-copy_exe

library: do-library

copy_library: do-copy_library


build_docker:
	docker build  -f Dockerfile -t lillymolprivate .
     
test_lillymol:
	docker container exec  lilly_mol  bash -c "cd test/ && ./run_all_test.sh"

start_s3:
	docker-compose up -d
	sleep 10
	echo 'run s3 commannd with: aws s3 --endpoint "http://localhost:4566" <s3 command>'
stop_s3:
	docker-compose down

	
