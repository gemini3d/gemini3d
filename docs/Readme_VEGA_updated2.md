# GEMINI Setup Guide for ERAU HPC system

---
---

## Installing GEMINI on ERAU's VEGA with Bash

### 1. Load modules
Start by loading the GCC, OpenMPI, LAPACK, CMake, and Python modules. In `~/.bash_profile` append the following:
```sh
module load gcc/8.5.0-gcc-8.5.0-cokvw3c
module load openmpi/5.0.2-gcc-8.5.0-diludms
module load netlib-lapack/3.11.0-gcc-8.5.0-hlxv33x
module load cmake/3.27.9-gcc-8.5.0-5bfhfkk
module load python/3.11.7-gcc-8.5.0-wfpoppf
```
These modules will load for every login Bash instance (`bash --login`).
They use the older GCC version (8.5.0) as this is the version that includes the LAPACK module.
If your terminal does not start a login shell (default in VS Code), or if you want to avoid logging out and back in, run
```sh
source ~/.bash_profile
```
You can instead do this in `~/.bashrc` but this is not common practise as this will load these modules in every bash instance, interactive or not. Add
```sh
if [ -f ~/.bashrc ]; then
    source ~/.bashrc
fi
```
to `~/.bash_profile` to run `~/.bashrc` on every login Bash.

### 2. Clone and build GEMINI
Clone the latest gemini3d repository and build it:
```sh
git clone https://github.com/gemini3d/gemini3d.git
cd gemini3d
cmake -B build
cmake --build build -j16
```
The build configuration (`-B build`) will generate the necessary files in the build directory.
This script will take a few minutes to complete; it will note multiple missing packages but will download source code for these to be compiled during the build step.
This step can take up to 10 minutes, depending on download times.
A full compile (`--build build`) will take approximately 5 – 10 minutes, depending on which packages need to be compiled; you will see many warnings, but these can be safely ignored.
Executables are placed in the build directory and can be run from there.
To ensure the binaries complied correctly, ensure you read, for example,
```sh
[100%] Built target gemini.bin
```

### 3. Test GEMINI installation
***Do not run ctest on VEGA login node!*** You will get yelled at. First download the applicable tests,
```sh
ctest --test-dir build --preset download
```
and ensure you read
```sh
100% tests passed, 0 tests failed out of 8
```
Next, run an interactive PBS job.
This interactive session requires all processors available in the node (on VEGA this is always 192).
To start the session, run
```sh
qsub -I -l walltime=1:00:00 -l nodes=1:ppn=192
```
and navigate back to the `gemini3d` directory. Then run the tests,
```sh
ctest --test-dir build
```
and ensure you read
```sh
100% tests passed, 0 tests failed out of 68
```
The amount of tests that run might vary, depending on how many are disabled.
Make sure to exit the interactive session when finished.

---
---

## Installing PyGEMINI on ERAU's VEGA with Bash

### 1. Load modules
If not already done so, in `~/.bash_profile` append the following:
```sh
module load python/3.11.7-gcc-8.5.0-wfpoppf
```
and source it,
```sh
source ~/.bash_profile
```

### 2. Create GEMINI Python environment
PyGEMINI cannot (and should not) be installed on the root VEGA python. Create a python environment in a local `.venvs` location:
```sh
mkdir ~/.venvs
python -m venv ~/.venvs/gemini
```
Next, if you expect to be working in this environment most of the time, append the following in `~/.bash_profile`:
```sh
source ~/.venvs/gemini/bin/activate
```
and then run
```sh
source ~/.bash_profile
```
Otherwise, simply source the `activate` file.

### 3. Install PyGEMINI
When inside the GEMINI environment (terminal starts with e.g. `(gemini) [username@vegaln1 ~]$`), simply run
```sh
pip install gemini3d
```

---
---

## Example PBS script and queue submission
### 1. Set up environment variables
It is convenient to export two environment variables in `~/.bashrc` that locate your fortran `gemini3d` and simulations directory, i.e.
```sh
export GEMINI_ROOT='<PATH_TO_GEMINI3D>'
export GEMINI_SIM_ROOT='<PATH_TO_SIMULTATION_DIR>'
```

### 2. Create PBS batch script
Inside your simulation directory, create a PBS batch script (e.g. `pbs.script`) containing the following:
```sh
# Command options:
#PBS -N example_simulation
#PBS -S /bin/bash
#PBS -q normalq
#PBS -l nodes=1:ppn=128
#PBS -l walltime=4:00:00
#PBS -o $HOME/logs/${PBS_JOBNAME}.${PBS_JOBID}.out
#PBS -e $HOME/logs/${PBS_JOBNAME}.${PBS_JOBID}.err
#PBS -V

# Modules to load:
module purge
module load gcc/8.5.0-gcc-8.5.0-cokvw3c
module load openmpi/5.0.2-gcc-8.5.0-diludms
module load netlib-lapack/3.11.0-gcc-8.5.0-hlxv33x

# Commands to run:
cp -r $GEMINI_SIM_ROOT/$PBS_JOBNAME $PBS_O_HOME/scratch
cp $GEMINI_ROOT/build/gemini.bin $PBS_O_HOME/scratch/$PBS_JOBNAME
cd $PBS_O_HOME/scratch/$PBS_JOBNAME
mpiexec gemini.bin . > $PBS_JOBNAME.out 2> $PBS_JOBNAME.err
cp -nr $PBS_O_HOME/scratch/$PBS_JOBNAME $GEMINI_SIM_ROOT
```
If all is well, you can remove the simulation directory from `~/scratch` when finished.

### 3. Submit your job to the queue
Submit your job to the queue,
```sh
qsub pbs.script
```
and check on its status with
```sh
qstat -u $USER
```
You can check the job status with more detail using `checkjob <JOB_ID>` or by reading
```sh
~/scratch/example_simulation/example_simulation.out
```
