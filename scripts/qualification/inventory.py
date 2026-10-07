"""Record the actual local build/Python/data inventory without claiming approval."""
import argparse
import hashlib
import importlib.metadata
import json
import platform
from pathlib import Path
import subprocess
import sys


def digest(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True);p.add_argument('--build',type=Path,required=True)
    p.add_argument('--package-cache',type=Path);p.add_argument('--data-root',type=Path,action='append',default=[])
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    cache={}
    for line in (a.build/'CMakeCache.txt').read_text().splitlines():
        if line and not line.startswith(('#','//')) and '=' in line and ':' in line.split('=',1)[0]:
            k,v=line.split('=',1);cache[k.split(':')[0]]=v
    libkeys=[k for k in cache if k.startswith(('HDF5_','BLAS_','LAPACK_','SCALAPACK_','MPI_')) and ('LIBRARY' in k or 'VERSION' in k)]
    linked=[]
    for key in sorted(libkeys):
        value=cache[key];entry=dict(key=key,value=value)
        path=Path(value)
        if path.is_file():entry.update(sha256=digest(path),resolved_path=str(path.resolve()))
        linked.append(entry)
    packages=[]
    if a.package_cache:
        if not a.package_cache.is_dir():raise ValueError('Requested package cache is missing')
        for path in sorted(a.package_cache.glob('*.deb')):
            proc=subprocess.run(['dpkg-deb','--field',str(path),'Package','Version','Architecture'],text=True,capture_output=True,check=True)
            packages.append(dict(archive=path.name,sha256=digest(path),metadata=proc.stdout.strip()))
    data=[]
    for root in a.data_root:
        if not root.exists():raise ValueError('Requested model/input path is missing: '+str(root))
        paths=[root] if root.is_file() else sorted(root.rglob('*'))
        for path in paths:
            if path.is_file():data.append(dict(root=str(root),path=path.name if root.is_file() else str(path.relative_to(root)),
                                              size=path.stat().st_size,sha256=digest(path)))
    toolchains=[]
    for key in ['CMAKE_C_COMPILER','CMAKE_CXX_COMPILER','CMAKE_Fortran_COMPILER','MPIEXEC_EXECUTABLE','CMAKE_COMMAND']:
        if key not in cache:continue
        proc=subprocess.run([cache[key],'--version'],capture_output=True,text=True,timeout=30)
        toolchains.append(dict(key=key,path=cache[key],returncode=proc.returncode,version=(proc.stdout+proc.stderr).strip()))
    try:
        import h5py
        python_hdf5=h5py.version.hdf5_version
    except ImportError:python_hdf5=None
    versions=sorted([dict(name=d.metadata['Name'],version=d.version) for d in importlib.metadata.distributions()],key=lambda x:x['name'].lower())
    rev=subprocess.run(['git','-C',str(a.source),'rev-parse','HEAD'],check=True,capture_output=True,text=True).stdout.strip()
    dirty=bool(subprocess.run(['git','-C',str(a.source),'status','--porcelain'],check=True,capture_output=True,text=True).stdout.strip())
    payload=dict(schema='gemini.qualification.inventory.1',scope='local test environment; deployment inventory requires target execution',
                 source_commit=rev,source_dirty=dirty,platform=platform.platform(),python=sys.version,python_packages=versions,
                 toolchains=toolchains,python_h5py_hdf5_version=python_hdf5,
                 build_libraries=linked,package_archives=packages,data_files=data,
                 release_license_approval=False,security_advisory_review_complete=False)
    a.output.write_text(json.dumps(payload,indent=2)+'\n')
    print(json.dumps(dict(python_packages=len(versions),library_records=len(linked),package_archives=len(packages),data_files=len(data))))
if __name__=='__main__':main()
