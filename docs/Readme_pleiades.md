# Running GEMINI on the NASA Pleiades system

Logging in to Pleiades-the SSH server will ask for 
[NAS password](https://www.nas.nasa.gov/hecc/support/kb/obtaining-and-changing-your-nas-password_127.html) 
and RSA passcode (fob or app).

```
ssh <username>@sfe6.nas.nasa.gov
ssh pfe
```


1. Set up compilers, wrappers, and environment variables

```
module load gcc/13.2
module use /nasa/modulefiles/testing
module load openmpi/4.1.6-toss4-gnu
export MPI_ROOT=/nasa/openmpi/4.1.6-gnu/bin
export CC=gcc CXX=g++ FC=gfortran
```

2. Pull the repository, configure, and build

```
git clone https://github.com/gemini3d/gemini3d/
cd gemini3d
git submodule update --init --recursive
cmake -B build
cmake --build build -j
```

3. example PBS script

```
#PBS -lselect=40:ncpus=28:mpiprocs=28:model=bro
#PBS -q normal
#PBS -N gemini3d
#PBS -l walltime=4:00:00
module use /nasa/modulefiles/testing
module load gcc/13.2
module load openmpi/4.1.6-toss4-gnu
cd $PBS_O_WORKDIR
mpiexec -np 1120 ./gemini.bin ./GDIround/ > log.out
```

4. submit the job

```
qsub pbs_script
qstat -u <username>
qdel <PID>
```

5. copy data to Pleiades

```
scp -oProxyCommand='ssh <nas_username>@sfe6.nas.nasa.gov ssh-proxy %h' <file1> <nas_username>@pfe.nas.nasa.gov:<file2>
```

6. copy data to long term storage

From the sfe:

```
ssh lou
cp -rv <file> /nobackup/<username>/
```
