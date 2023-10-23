# Constrained single-particle tomography (CSPT)

## 1. First, follow the instructions to build the [ETTK](https://github.com/nextpyp/ettk) and [cistem-1.0](https://github.com/nextpyp/cistem-1.0) packages using the intel compiler

## 2. Build cisTEM library (statically linked)
```
# Go to the cistem-1.0 repo
cd cistem-1.0

# Use cistem container and enable intel compiler, for example:
singularity shell -B /opt/apps cistem_dev_env_latest.sif
source  /opt/apps/rhel8/intel-2020/compilers_and_libraries/linux/bin/compilervars.sh intel64

cd src
./build_cspt_lib.sh

# Exit the container
exit
```
(If the build was successful, you should see the file `cspt_lib.a`)

## 3. Build CSP 
```
# Set up the path to the ettk & cistem builds
export ETTK_PATH=path_to_your_ettk_project
export CISTEM_PATH=path_to_your_cistem-1.0_project

# Use ettk container for building
singularity shell ettk-devel.sif

# Go to csp repo and link header file for cistem library
cd csp/src/
ln -s ${CISTEM_PATH}/src/programs/refine3d_cspt/refine3d_cspt.h .

rm CMakeCache.txt && rm -rf CMakeFiles

ccmake -D ITK_DIR=${ETTK_PATH}/ettk/external/InsightToolkit-4.2.1/build \
       -D VTK_DIR=${ETTK_PATH}/ettk/external/VTK5.10.1/build 
       -D BLITZ_DIR=${ETTK_PATH}/ettk/external/blitz-0.9 
       -D CMAKE_CXX_COMPILER=/opt/apps/rhel8/intel-2020/compilers_and_libraries_2020.4.304/linux/bin/intel64/icpc 
       -D CMAKE_C_COMPILER=/opt/apps/rhel8/intel-2020/compilers_and_libraries_2020.4.304/linux/bin/intel64/icc 
       -D CMAKE_CXX_FLAGS="-no-multibyte-chars" 
       -D CMAKE_C_FLAGS="-no-multibyte-chars"
       -D CISTEM_LIB=${CISTEM_PATH}/src/cspt_lib.a  
       .

# The command will bring you to an interactive window, run the following in order
(1) c
(2) c
(3) g

make csp
```
