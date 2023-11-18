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

# A default target that will probably work in most cases.
# Note it does not build BerkeleyDB dependent tools or Python bindings.
default:
	bash -c 'if [[ -z "$$(type -p bazelisk)" && -z "$$(type -p bazel)" ]] ; then echo "No bazel/bazelisk, see README.md" && exit 1 ; fi'
	echo "Default build does not build targets 'berkeleydb' and 'python'"
	cd src && ./update_bazel_configs.sh
	cd src && ./build_third_party.sh
	cd src && ./build_from_src.sh

all:
	bash -c 'if [[ -z "$$(type -p bazelisk)" && -z "$$(type -p bazel)" ]] ; then echo "No bazel/bazelisk, see README.md" && exit 1 ; fi'
	cd src && BUILD_BDB=1 BUILD_PYTHON=1 BUILD_VENDOR=1 ./update_bazel_configs.sh
	cd src && BUILD_BDB=1 BUILD_PYTHON=1 BUILD_VENDOR=1 ./build_third_party.sh
	cd src && BUILD_BDB=1 BUILD_PYTHON=1 BUILD_VENDOR=1 ./build_from_src.sh

berkeleydb:
	bash -c 'if [[ -z "$$(type -p bazelisk)" && -z "$$(type -p bazel)" ]] ; then echo "No bazel/bazelisk, see README.md" && exit 1 ; fi'
	cd src && BUILD_BDB=1 BUILD_PYTHON=1 BUILD_VENDOR=1 ./update_bazel_configs.sh
	cd src && BUILD_BDB=1 ./build_third_party.sh
	cd src && BUILD_BDB=1 ./build_from_src.sh

python:
	bash -c 'if [[ -z "$$(type -p bazelisk)" && -z "$$(type -p bazel)" ]] ; then echo "No bazel/bazelisk, see README.md" && exit 1 ; fi'
	cd src && BUILD_BDB=1 BUILD_PYTHON=1 BUILD_VENDOR=1 ./update_bazel_configs.sh
	cd src && BUILD_PYTHON=1 ./build_third_party.sh
	cd src && BUILD_PYTHON=1 ./build_from_src.sh

vendor:
	bash -c 'if [[ -z "$$(type -p bazelisk)" && -z "$$(type -p bazel)" ]] ; then echo "No bazel/bazelisk, see README.md" && exit 1 ; fi'
	cd src && BUILD_BDB=1 BUILD_PYTHON=1 BUILD_VENDOR=1 ./update_bazel_configs.sh
	cd src && BUILD_VENDOR=1 ./build_third_party.sh
	cd src && BUILD_VENDOR=1 ./build_from_src.sh

build_docker:
	docker build -f Dockerfile -t lillymolprivate .
     
test_lillymol:
	docker container exec lilly_mol bash -c 'cd bin/Linux/ && echo "Number of executables built: `ls -1 | wc -l`" && ls -1'
	docker container exec lilly_mol bash -c "cd test/ && ./run_all_test.sh"
	docker container exec lilly_mol bash -c "cd src/ && ./run_python_unit_tests.sh 2>&1"

start_s3:
	docker-compose up -d
	sleep 30
	echo 'run s3 commannd with: aws s3 --endpoint "http://localhost:4566" <s3 command>'
stop_s3:
	docker-compose down
