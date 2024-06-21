#! /bin/bash

if [ -z "$LILLYMOL_HOME" ] || [ -z "$BUILD_DIR" ]
then 
    # undefined BIN_DIR
    echo "System variables LILLYMOL_HOME and BUILD_DIR are required for running the test"
    echo "Please export LILLYMOL_HOME(local path to LillyMol code)"
    echo "Please export BUILD_DIR(the folder name under the bin folder after build)"
    echo "Example: export LILLYMOL_HOME=/home/user/LillyMol"
    echo "Example: export BUILD_DIR=Linux-gcc-7.2.1" 
    exit 1
else
    BIN_DIR=$LILLYMOL_HOME/bin/$BUILD_DIR
fi

shared_data_dir=$LILLYMOL_HOME/data
command=$BIN_DIR/gfp_nearneighbours_single_file
case_id="Case 1"
echo "Testing:  ${command}"

if [ ! -e "${command}" ]
then
  echo "Executable ${command} is not found"
  exit 1
fi

golden=out/out.txt

stdout=out.txt
stderr='stderr'

# Support linux and mac 
if [[ "$OSTYPE" == "linux-gnu" ]]; then
  golden='out/linux/out.txt'
elif [[ "$OSTYPE" == "darwin"* ]]; then
  golden='out/osx/out.txt'
else
  echo "OS is not supported"
  golden='out/linux/out.txt'  # Might work...
fi

diff_tool='../../fileDiff.sh'

# Need to remove neighbours with a distance of 0.2 because that can vary across
# library versions - the -T 0.2 is an exact floating point comparison.
# This does not really solve the problem, but lessens the probability
# of issues. 
# Alternative would be to write a custom diff tool
${command} -p -z -T 0.2 -F FPDSC,w=0.2 -F NCSELW,nc,w=0.8 -V a=0.3 -V b=1.7 -j 3 ${shared_data_dir}/pubchem.gfp 2>${stderr} | sed -e '/ 0\.2$/d' > ${stdout}
${diff_tool} ${stdout} ${golden}

if [ $? -eq 1 ]
then
  echo "$case_id : TEST PASS"
else
  echo "$case_id : TEST FAIL"
  diff -w ${stdout} ${golden} | head -n 20
fi

rm ${stdout} ${stderr}
