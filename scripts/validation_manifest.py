#!/usr/bin/env python3
"""Create or verify the SHA-256 inventory of a completed validation campaign."""
import argparse
import hashlib
import json
from datetime import datetime
from pathlib import Path
import subprocess
from zoneinfo import ZoneInfo


def git(repo, *args):
    return subprocess.check_output(["git", "-C", str(repo), *args], text=True).strip()


def inventory(root):
    return [{"path": p.relative_to(root).as_posix(), "bytes": p.stat().st_size,
             "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
            for p in sorted(root.rglob("*"))
            if p.is_file() and p != root / "manifest.json"]


def verify(root):
    data = json.loads((root / "manifest.json").read_text())
    actual = inventory(root)
    if data["files"] != actual:
        expected = {r["path"]: r for r in data["files"]}
        found = {r["path"]: r for r in actual}
        differences = [p for p in sorted(expected.keys() | found.keys())
                       if expected.get(p) != found.get(p)]
        raise ValueError("Manifest mismatch: " + ", ".join(differences))
    return len(actual)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--efg-repo", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--fem-repo", type=Path)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    root = args.campaign.resolve()
    if args.verify:
        print(f"PASS: {verify(root)} file hashes verified")
        return
    sources = {}
    for name, repo in (("efg", args.efg_repo), ("fem", args.fem_repo)):
        if repo is not None:
            sources[name] = {"commit": git(repo, "rev-parse", "HEAD"),
                             "remote": git(repo, "remote", "get-url", "origin"),
                             "branch": git(repo, "branch", "--show-current")}
    cases = {}
    for family in ("efg", "fem"):
        for meta in sorted((root / family).glob("*/run_metadata.json")):
            cases[meta.parent.relative_to(root).as_posix()] = json.loads(meta.read_text())
    data = {"schema_version": 1,
            "created_at": datetime.now(ZoneInfo("America/Sao_Paulo")).isoformat(),
            "sources": sources, "cases": cases,
            "hash_policy": "SHA-256 of every campaign file except manifest.json itself",
            "files": inventory(root)}
    (root / "manifest.json").write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n")
    print(f"Manifest written: {len(data['files'])} files")


if __name__ == "__main__":
    main()
