########################################################################################
                                       Warning: 
Please pay special attention to the order of arguments and options when
you are running command under Mac
1. All the command options shall be placed right after the command name.
2. All the other input information for the command shall be appended to the end of the command after options
Sample:
fetch_smiles_quick -j -c 1 -C 2 -X notInRecord -Y notInIdentifier record.w structures.sm  

########################################################################################         


Welcome to the Eli Lilly LillyMol implementation.

REQUIREMENT:

This software requires following packages to build
1. GCC >= 6.2.0 (see https://gcc.gnu.org/install/index.html)
Example command to intall gcc: `brew install gcc`

2. zlib >= 1.2.11 (see http://www.zlib.net/)
Example command to install zlib: `brew install zlib`
You need to define the location for zlib.a in makefile.public.*


Note:
Following command can be used to install commandline build tool on Mac OS X
`xcode-select --install`

Following command can be used to install the homebrew tool if it does not exist
`ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"`

BUILD:
1. Start terminal

2. Pull down the code from repo

3. Enter the root directory of the code

4. Update the path for the zlib in the makefile.public.OSX-gcc-8

5. Run makeall_osx.sh (Skip rest of the steps if you run this command)

6. Export compilers
 * `export CC=gcc-8`
 * `export CXX=g++-8`
 * `export FC=gfortran-8`

7. Export the UNAME for makefile
 * `export UNAME=OSX-gcc-8`


8. Alternatively, you can run following commands:
```bash
    make veryclean
    make copy_include
    make library
    make copy_library
    make exe
    make copy_exe
```

EXECUTION:
See Wiki page for sample commands
See the example folder for data used in the sample commands
Example to verify a generated command:
```bash
    cd bin/OSX-gcc-8/
    ./common_names
```
Note: This command shall print out the help menu on screen if it is built successfully
      This approach can be used to verify if each command is built successfully
Example to run sample commands:
```bash
    cd example/common_names/
    ../../bin/OSX-gcc-8/common_names input1.smi input2.smi -S ./output -s 10000 -r 10000 -D + -v
```


DIRECTORY:
src:                             source code
example:                         data for sample commands (see Wiki page)
test:                            test scripts for each command
bin(generated after build):      all generated executables
lib(generated after build):      all generated library files
include(generated after build):  shared include files
contrib:                         legacy tools


TEST:
See the test folder for the test case for each command
Example to run a test:
```
    cd test/common_names/case_1/
    ./test_case_1.sh
```
Note: Test shall print out TEST PASS if it is successful, otherwise it shall
      print out TEST FAIL

LICENSE:
Consult the LICENSE file for details of the license
