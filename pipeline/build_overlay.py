"""Build a patched DIALS Python package under a cross-platform file lock."""

import fcntl
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent.parent


def patches():
    pattern = os.environ.get("PATCHES", "*.patch")
    maximum = int(os.environ.get("PATCH_MAX", "15"))
    return [p for p in sorted((ROOT / "dials-patch").glob(pattern)) if int(p.name[:4]) <= maximum]


def main():
    overlay = Path(os.environ.get("OVERLAY", ROOT / ".venv-dials" / "overlay"))
    python = ROOT / ".venv-dials" / "env" / "bin" / "python"
    if not python.exists():
        raise SystemExit("DIALS environment missing: run bash reproduce.sh env first")
    digest = hashlib.sha256()
    selected = patches()
    if not selected:
        raise SystemExit("No DIALS patches selected")
    for p in selected:
        digest.update(p.name.encode())
        digest.update(p.read_bytes())
    checksum = digest.hexdigest()
    overlay.parent.mkdir(parents=True, exist_ok=True)
    with open(str(overlay) + ".lock", "w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        stamp = overlay / ".patches.sha256"
        if stamp.exists() and stamp.read_text().strip() == checksum:
            return
        source = subprocess.check_output([str(python), "-c", "import dials; print(dials.__path__[0])"], text=True).strip()
        temporary = Path(str(overlay) + ".building")
        if temporary.exists():
            shutil.rmtree(temporary)
        shutil.copytree(source, temporary / "dials", ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
        for p in selected:
            subprocess.run(["patch", "--batch", "-p2", "--forward", "-i", str(p)], cwd=temporary, check=True, stdout=subprocess.DEVNULL)
        (temporary / ".patches.sha256").write_text(checksum + "\n")
        if overlay.exists():
            shutil.rmtree(overlay)
        temporary.rename(overlay)


if __name__ == "__main__":
    main()
