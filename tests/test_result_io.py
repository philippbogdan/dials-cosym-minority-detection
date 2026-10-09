import gzip
import json

from cf.result_io import glob_results, open_result, result_exists


def test_gzip_records_are_read_without_materialising(tmp_path):
    path = tmp_path / "shard_000.jsonl"
    compressed = tmp_path / "shard_000.jsonl.gz"
    compressed.write_bytes(gzip.compress(b'{"id": "trial", "value": 1.23456789012345}\n'))
    assert result_exists(path)
    with open_result(path) as fh:
        assert json.load(fh) == {"id": "trial", "value": 1.23456789012345}
    assert not path.exists()


def test_fresh_records_take_precedence(tmp_path):
    native = tmp_path / "shard_000.jsonl"
    compressed = tmp_path / "shard_000.jsonl.gz"
    compressed.write_bytes(gzip.compress(b'{"value": 1}\n'))
    with open_result(native, "w") as fh:
        fh.write('{"value": 2}\n')
    assert glob_results(str(tmp_path / "shard_*.jsonl")) == [str(native)]
    with open_result(native) as fh:
        assert json.load(fh) == {"value": 2}
