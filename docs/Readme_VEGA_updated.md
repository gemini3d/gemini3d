# GEMINI Setup and Build Guide for ERAU HPC system

## Overview
This outlines the steps to build and install GEMINI3D and its dependencies on an ERAU HPC system.

Before building GEMINI, ensure the following dependencies are available:
- **GCC 8.5.0**
- **CMake**
- **OpenMPI**
- **Intel OneAPI Compilers** (to work across nodes)
- **LAPACK**

### 1. Load Necessary Modules
First, load the required modules for the build process. This ensures that your system has access to the necessary compilers and libraries:

```
module load intel-oneapi-compilers-classic
module load intel-oneapi-mkl
module load intel-oneapi-mpi
module load gcc/8.5.0-gcc-8.5.0-cokvw3c
module load cmake
module load openmpi
```
### 2. Clone the GEMINI3D Repository
Choose the location where you want to download the repository. For this example, we'll use ~/Projects:

```
cd ~/Projects

Use git to clone the repository:
```
git clone https://github.com/gemini3d/gemini3d.git

Ensure you have git installed and configured before attempting to clone.

### 3. Build LAPACK Locally
#### Create a directory for LAPACK under .local
```
mkdir ~/.local/lapack

#### Download the LAPACK source code
```
cd ~/.local/lapack
git clone https://github.com/Reference-LAPACK/lapack.git

#### Build and install LAPACK
```
cmake -B build -DCMAKE_INSTALL_PREFIX=$HOME/.local/lapack
cmake --build build -j16 (-j16 specifies the number of cores to use for the build (16 in this case). You can adjust this depending on how many cores your system has.)
cmake --install build

#### Verify Installation
Confirm the installation by checking the contents of the .local/lapack/ directory.

#### Set LAPACK Environment Variable
Set the environment variable to point to the installed LAPACK library.
```
export LAPACK_ROOT=$HOME/.local/lapack/
To make this consistent across sessions, add the line above to your .bashrc

### 4. Build GEMINI
Navigate to the cloned gemini3d directory (if not already cloned, perform step 2)
```
cd ~/Projects/gemini3d

#### Run cmake to configure the build
cmake -B build
This will configure the build system and generate the necessary files in the build directory. This script will take a few minutes to complete; it will note multiple missing packages but will download source code for these to be compiled during the build step. This step can take up to 10 minutes depending on download times. If this completes successfully you can move on to the compile step:

#### Compile
cmake --build build -j16 (adjust j# accordingly based on how many cores your system has)
A full compile will take approximate 5-10 minutes depending on which packages need to be compiled; you will see lots of warnings but these can be safely ignored. Executables are placed in the build directory and can be run from there.

## To avoid manually loading the required modules every time you log in, you can add the module load commands to your .bashrc file. This ensures that the modules are automatically loaded when you start a new session.
```
# .bashrc
# Source global definitions
if [ -f /etc/bashrc ]; then
. /etc/bashrc
fi

# Path definitions
source /apps/spack/share/spack/setup-env.sh
source $(spack location -i lmod)/lmod/lmod/init/bash
#cd /apps/spack/share/spack
#source setup-env.sh
cd $HOME

# modules
module load intel-oneapi-compilers-classic
module load intel-oneapi-mkl
module load intel-oneapi-mpi
module load cmake

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/ihome/home2/zettergm/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/ihome/home2/zettergm/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/ihome/home2/zettergm/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/ihome/home2/zettergm/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

If your system runs .bash_profile on login (as is common in many systems, including Vega), you can include a line in your .bash_profile to source .bashrc. This ensures that the settings and environment variables in .bashrc are applied during login.
# .bash_profile
  
# Get the aliases and functions
if [ -r ~/.bashrc ]; then
        source ~/.bashrc
fi
```
This checks if .bashrc exists and sources it to apply the environment settings defined there.

# Running GEMINI on the ERAU HPC
To start an interactive session
```
qsub -I -l walltime=24:000:00 -l nodes=1:ppn=192

For longer or more expensive simulations you will need to use the queueing system, which requires a script that runs the code to be placed in the build directory where executables reside or by adding the build directory as working directoty. An example of such a script is shown below:
```
#!/bin/bash
# This is comment
#PBS -q longq
#PBS -l walltime=120:00:00
#PBS -l nodes=1:ppn=192
#PBS -N GEMINI_mooreOK
# Not necessary to set wall time, 
# longq: 120 hours, 42 nodes (each with 2x AMD 96-core processors)
# normalq: 24 hours, 42 nodes

# In you want to send job status to e-mail:
#PBS -M erauID@erau.edu
#PBS -m aeb
cd $PBS_O_WORKDIR

#PBS -e /scratch/zettergm/GEMINI3D/pbs_errors.out
#PBS -o /scratch/zettergm/GEMINI3D/pbs_output.out

# Load modules during job submission!
module load gcc/8.5.0-gcc-8.5.0-cokvw3c
module load openmpi

# Run the program
mpirun -np $PBS_NP ./gemini.bin ../simulations/mooreOK3D_hemis_medres_corrected_control/

