"""Bounded, self-contained HDF5 numeric data for qualification inputs.

This profile does not make a native HDF5 parser a sandbox. Inspect files from
untrusted sources in a separate constrained process/host before admission.
"""
import os
os.environ['HDF5_PLUGIN_PRELOAD']='::'
import h5py

MAX_DATASET_BYTES=256*1024**2
BUILTIN_FILTERS={1,2,3}  # deflate, shuffle, Fletcher32; no dynamic plugins
# h5py exposes the plugin search path, not H5PLset_loading_state. Remove search
# paths as well when h5py was imported earlier by an embedding process. Dataset
# inspection still rejects every non-allowlisted filter before reading payload.
for _ in range(h5py.h5pl.size()):
    h5py.h5pl.remove(0)


def inspect_dataset(obj,max_bytes=MAX_DATASET_BYTES):
    if not isinstance(obj,h5py.Dataset):raise ValueError('Expected HDF5 dataset')
    if obj.is_virtual:raise ValueError('Virtual HDF5 datasets are unsupported')
    props=obj.id.get_create_plist()
    if props.get_external_count():raise ValueError('External HDF5 storage is unsupported')
    if obj.size*obj.dtype.itemsize>max_bytes:raise ValueError('HDF5 dataset exceeds allocation budget')
    for i in range(props.get_nfilters()):
        if props.get_filter(i)[0] not in BUILTIN_FILTERS:raise ValueError('Unsupported HDF5 filter')
    return obj


def dataset(group,name):
    # Check every path component before dereferencing; an intermediate group
    # link is just as capable of leaving the file as a linked final dataset.
    current=group
    for key in name.split('/'):
        if not key or key not in current or not isinstance(current.get(key,getlink=True),h5py.HardLink):
            raise ValueError('Missing or indirect HDF5 path: '+name)
        current=current[key]
    return inspect_dataset(current)


def numeric(group,name,shape=None):
    obj=dataset(group,name)
    if obj.dtype.kind not in 'biuf':raise ValueError('Expected real numeric HDF5 data: '+name)
    if shape is not None and obj.shape!=shape:raise ValueError('HDF5 shape mismatch: '+name)
    return obj[...]
