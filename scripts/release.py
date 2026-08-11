#!/usr/bin/env python3
"""Prepare, publish, and push a smafa release."""

import argparse
import re
import subprocess
import sys
from pathlib import Path


def run(*command: str) -> None:
    print(f"+ {' '.join(command)}")
    subprocess.run(command, check=True)


def capture(*command: str) -> str:
    return subprocess.check_output(command, text=True).strip()


def fail(message: str) -> None:
    raise SystemExit(f"error: {message}")


def update_version(version: str) -> None:
    manifest = Path("Cargo.toml")
    updated, count = re.subn(
        r'(?m)^version\s*=\s*"[^"]+"',
        f'version = "{version}"',
        manifest.read_text(),
        count=1,
    )
    if count != 1:
        fail("could not find the package version in Cargo.toml")
    manifest.write_text(updated)


def update_changelog(version: str) -> None:
    path = Path("CHANGELOG.md")
    text = path.read_text()
    if re.search(rf"(?m)^## Version {re.escape(version)}$", text):
        fail(f"CHANGELOG.md already contains Version {version}")
    match = re.search(r"(?ms)^## Unreleased\s*\n(.*?)(?=^## )", text)
    if not match or not match.group(1).strip():
        fail("CHANGELOG.md has no Unreleased entries")
    replacement = f"## Unreleased\n\n## Version {version}\n\n{match.group(1).strip()}\n\n"
    path.write_text(text[: match.start()] + replacement + text[match.end() :])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", required=True, help="release version, e.g. 1.2.3")
    parser.add_argument("--no-commit", action="store_true")
    args = parser.parse_args()

    if not re.fullmatch(r"\d+\.\d+\.\d+(?:[-+][A-Za-z0-9.-]+)?", args.version):
        fail(f"invalid version: {args.version}")
    if capture("git", "status", "--porcelain"):
        fail("git working tree is not clean")

    tag = f"v{args.version}"
    if capture("git", "tag", "--list", tag):
        fail(f"tag already exists: {tag}")

    answer = input("Are the Unreleased changelog entries ready? [y/N] ").lower()
    if answer not in {"y", "yes"}:
        fail("update CHANGELOG.md first")

    update_changelog(args.version)
    update_version(args.version)
    run("cargo", "update", "--workspace")
    run("cargo", "test")
    run("dist", "plan")

    if args.no_commit:
        print("Stopped before commit, tag, publish, and push.")
        return

    run("git", "add", "Cargo.toml", "Cargo.lock", "CHANGELOG.md")
    run("git", "commit", "-m", f"Release {tag}")
    run("git", "tag", tag)
    run("cargo", "publish")
    run("git", "push", "origin", "HEAD", "--tags")
    print("Release published. Next update the bioconda recipe.")


if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as error:
        sys.exit(error.returncode)
