from pathlib import Path

from orbit_benchmark.fair import sha256_file, write_manifest


def test_manifest_contains_checksum(tmp_path: Path):
    target = tmp_path / "table.csv"
    target.write_text("a,b\n1,2\n", encoding="utf-8")
    manifest = write_manifest(tmp_path, tmp_path / "manifest")
    text = manifest.read_text(encoding="utf-8")
    assert sha256_file(target) in text
    assert "table.csv" in text

