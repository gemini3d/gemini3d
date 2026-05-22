# GEMINI Setup Guide for KTH HPC System

## Installing GEMINI on KTH's Dardel with Bash

Login via ssh and change directory to your workspace. 

You can start by working in the scratch space, this is for temporary files that are deleted if not worked with in 30 days, later move it to the project space. Or directly set it in your project space. 

Dardel uses parallel storaging by the name of klemming and to access them: 

```bash
cd /cfs/klemming/scratch/u/username
cd /cfs/klemming/projects/snic/projectname
```
more information on https://support.pdc.kth.se/doc/data_management/klemming/

### 1. Clone GEMINI

Clone the latest gemini3d repository

```bash
git clone https://github.com/gemini3d/gemini3d.git
cd gemini3d
```
### 2. Check and load softwares present in environment
Staying in the source code directory, start by listing the softwares present in the environment you are working in.

```bash
module list
```

Usually the following modules will already by installed in your environment i.e. the **Cray programming environment** 
<br>
gcc, gfortran, mpich, cmake, libsci<br>. 
Incase any module is missing, check its availability by  

```bash
module avail *insert name of module*
module avail gcc #example
```
Then, load using

```bash
module load cpe/24.11
module load nano # this module will have to be loaded everytime you enter the system 
```
The cpe/24.11 loads all necessary default compilers and environments. Nano module is used to modify scripts.

### 3. Setting variables and building
To ensure the system uses the correct compilers from the loaded environment, set the following variables:

```bash
export CC=cc
export CXX=CC
export FC=ftn
```
GEMINI uses the GNU programming environment instead of the default Cray environment at dardel. So, we will swap the compiler environment before build command:

```bash
module swap PrgEnv-cray PrgEnv-gnu
cmake -B build -DFETCHCONTENT_TRY_FIND_PACKAGE_MODE=NEVER
cmake --build build -j16
```
the **-B build*** generates necessary executing files in the build directory. This takes some time, if any error,check for discrepancy/unavailability in modules. The following cmake flag forces that our system to build all dependencies needed from scratch, and not let our exisiting local system to take any precedence in packages needed for building. 

### 4. Testing built

It is recommended to test the build, however running ctest --test-dir build  on login node causes some failures with standalone MPI test which leads to failure of a few following dependent tests as well. 
However, we can move to running extensive jobs through queuing system if the compilation has been successful. 

## Installing PyGEMINI on KTH's Dardel with Bash

The default python is 2.7.9 thus, to load a higher version of python

```bash
module load python3.11
```

### 1. Create Pygemini 

as per information on the https://github.com/gemini3d/pygemini

Pygemini plotting is independent of the gemini3d build. One can directly clone it and plot outputs from simulation files. For other aspects, like running simulations, we will need to build gemini3d first for which **set the root** in python script to map to the built gemini3d. 

```python
import os
os.environ["GEMINI_ROOT"]="~/Projects/gemini3d/build/_deps/msis-build/"
```
The following steps are available: https://github.com/gemini3d/pygemini


## Submitting jobs

Dardel uses the slurm workload manager. Post build, for generating simulation files, we will submit our job to that manager. Run the following from the gemini3d source directory.

```bash
module load nano
nano gemini_example.slurm # your_file_name.slurm or .sh
```
### slurm script 

an example slurm script is as follows aside from the examples here: https://support.pdc.kth.se/doc/basics/quickstart/

``` bash
#!/bin/bash -l
#SBATCH -A naissXXXX-XX-XXXX #enter project name this is present on your supr account
#SBATCH -J gemini_arcs
#SBATCH -p main
#SBATCH -t 10:00:00 

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24

module load PrgEnv-cray
module load cray-mpich

EXE=/cfs/klemming/scratch/u/username/gemini3d/build/gemini.bin # basically the location where gemini.bin in build of gemini3d is stored
SIM=/cfs/klemming/scratch/u/username/pygemini/simulation/arcs_dist/ # the last hash is important you want it to be able to access all the information it needs regarding inputs and calculation. This is the location of the simulation directory

srun $EXE $SIM
```
### Bash commands to implement after creating above files
```bash
sbatch your_file_name.slurm
```
You may check the status of the file via 
```bash
sstat --jobs=your_job_id
``` 
to check on running file execution 
```bash
tail -f slurm-your_job_id.out
```
for past job 
```bash
sacct --jobs=your_job-id
#or
sacct -j you_job-id --format=JobID,State,ExitCode
```
check for the queued jobs 
```bash
squeue -u $USER
```
if ever your job is killed, or finishes early etc. Remember to cancel that job to allow that allocated node to close for the "extra time" not in use. 

```bash
scancel your_job-id
```

