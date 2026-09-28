"""Immutable run provenance and exact-reload numerical storage."""
import csv
import hashlib
import json
import platform
import subprocess
import uuid
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import scipy
from . import __version__


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def json_default(x):
    if isinstance(x,np.generic): return x.item()
    if isinstance(x,Path): return str(x)
    raise TypeError(type(x).__name__)


def write_json(path,obj):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    payload=json.dumps(obj,indent=2,sort_keys=True,default=json_default,allow_nan=True)
    path.write_text(payload+'\n')
    if json.loads(path.read_text()) != json.loads(payload):
        raise AssertionError('JSON reload mismatch')


def save_npz(path,**arrays):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists(): raise FileExistsError(path)
    np.savez_compressed(path,**arrays)
    with np.load(path,allow_pickle=False) as reloaded:
        if set(reloaded.files)!=set(arrays):raise AssertionError('NPZ keys mismatch')
        for key,val in arrays.items():
            a=np.asarray(val)
            equal=np.array_equal(a,reloaded[key],equal_nan=True) if np.issubdtype(a.dtype,np.inexact) else np.array_equal(a,reloaded[key])
            if not equal:raise AssertionError(f'Exact NPZ reload failed: {path}/{key}')


def write_csv(path,rows):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=keys);wr.writeheader();wr.writerows(rows)
    with path.open() as f:
        read=list(csv.DictReader(f))
    if len(read)!=len(rows):raise AssertionError('CSV reload length')
    for original,loaded in zip(rows,read):
        for k,v in original.items():
            if isinstance(v,(float,np.floating)) and not np.isnan(v) and float(loaded[k])!=v:
                raise AssertionError(f'CSV numeric roundtrip: {k}')


def read_csv(path):
    with Path(path).open() as f: rows=list(csv.DictReader(f))
    for row in rows:
        for k,v in row.items():
            try: row[k]=float(v)
            except ValueError: pass
    return rows

