"""A changed canonical table must stop validation before results are used."""
import importlib.util
import json
from pathlib import Path
import pytest


def test_checksum_mutation_stops_release(tmp_path):
    path=Path(__file__).resolve().parents[1]/'scripts/release.py'
    spec=importlib.util.spec_from_file_location('release_gate',path)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    target=tmp_path/'data/processed/response.csv';target.parent.mkdir(parents=True)
    target.write_text('H\n1.0\n')
    (tmp_path/'data/manifest.json').write_text(json.dumps({'data/processed/response.csv':{'sha256':mod.sha(target)}}))
    target.write_text('H\n2.0\n')
    with pytest.raises(RuntimeError,match='checksum mismatch'):
        mod.verify_data(tmp_path)
