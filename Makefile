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

.PHONY: default
.PHONY: all

# Determine the operating system
UNAME := $(shell uname)

ifeq ($(UNAME),Darwin)
	BUILD_SCRIPT := build_macos.zsh
else
	BUILD_SCRIPT := build_linux.sh
endif

# Determine REPO_HOME
REPO_HOME := $(CURDIR)

# A default target that will probably work in most cases.
# Note it does not build BerkeleyDB dependent tools or Python bindings.
default:
	@echo "Build platform: $(UNAME)"
	cd src && REPO_HOME=$(REPO_HOME) ./$(BUILD_SCRIPT) 

all:
	@echo "Build platform: $(UNAME)"
	cd src && REPO_HOME=$(REPO_HOME) BUILD_BDB=1 BUILD_PYTHON=1 BUILD_XGBOOST=1 BUILD_VENDOR=1 ./$(BUILD_SCRIPT)

advance:
	@echo "Build platform: $(UNAME)"
	cd src && REPO_HOME=$(REPO_HOME) BUILD_BDB=1 BUILD_PYTHON=1 ./$(BUILD_SCRIPT)

vendor:
	@echo "Build platform: $(UNAME)"
	cd src && REPO_HOME=$(REPO_HOME) BUILD_VENDOR=1 ./$(BUILD_SCRIPT)

build_docker:
	docker build -f Dockerfile -t lillymolprivate .
     
test_lillymol:
	docker container exec lilly_mol bash -c 'cd bin/Linux/ && echo "Number of executables built: `ls -1 | wc -l`" && ls -1'
	docker container exec lilly_mol bash -c "cd test/ && ./run_all_test.sh"
	docker container exec lilly_mol bash -c "cd src/ && ./run_python_unit_tests.sh 2>&1"

start_s3:
	docker-compose up -d
	sleep 30
	@echo 'run s3 commannd with: aws s3 --endpoint "http://localhost:4566" <s3 command>'

stop_s3:
	docker-compose down
