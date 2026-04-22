#!/usr/bin/env python3
"""
Check GH Actions status of last run and get last commit hash for each subproject.
Compare to the last commit hash in libraries.json. If they differ, print a warning.
"""

import subprocess
import json
import functools
from pathlib import Path
import argparse

root = Path(__file__).parents[1]


@functools.cache
def get_repo_name(repo_archive_url: str) -> str:
    """
    Extract username/repo from the repo URL. For example, for
    https://github.com/gemini3d/glow/archive/abc123.tar.gz
    it will return username/repo.
    """

    return get_repo_url(repo_archive_url).partition("github.com/")[2]

@functools.cache
def get_repo_url(repo_archive_url: str) -> str:
    """
    Extract the base repo URL from the archive URL. For example, for
    https://github.com/gemini3d/glow/archive/abc123.tar.gz
    it will return https://github.com/gemini3d/glow.
    """
    return repo_archive_url.partition("/archive/")[0]


def get_last_commit_hash(repo_archive_url: str) -> str | None:
    repo = get_repo_name(repo_archive_url)

    result = subprocess.run(
        ["gh", "api", f"repos/{repo}/commits", "--jq", ".[0].sha"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if result.returncode != 0:
        print(f"Error fetching last commit hash for {repo}")
        return None
    return result.stdout.strip() or None


@functools.cache
def get_default_branch(repo_archive_url: str) -> str | None:
    repo = get_repo_name(repo_archive_url)
    result = subprocess.run(
        ["gh", "api", f"repos/{repo}", "--jq", ".default_branch"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if result.returncode != 0:
        print(f"Error fetching default branch for {repo}")
        return None
    return result.stdout.strip() or None


def get_pinned_commit_hash(repo_archive_url: str) -> str | None:
    archive_suffix = repo_archive_url.partition("/archive/")[2]
    if not archive_suffix:
        return None
    return archive_suffix.partition(".")[0] or None


def check_github_actions_status(repo_archive_url: str) -> bool:
    repo = get_repo_name(repo_archive_url)
    branch = get_default_branch(repo_archive_url)
    if branch is None:
        return False

    result = subprocess.run(
        [
            "gh",
            "run",
            "list",
            "--repo",
            repo,
            "--branch",
            branch,
            "--limit",
            "1",
            "--json",
            "status,conclusion",
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if result.returncode != 0:
        print(f"Error fetching GitHub Actions status for {repo}")
        return False

    runs = json.loads(result.stdout)
    if not runs:
        print(f"No GitHub Actions runs found for {repo}")
        return False

    last_run = runs[0]
    return last_run["status"] == "completed" and last_run["conclusion"] == "success"


def main():
    parser = argparse.ArgumentParser(description="Check subproject CI status and commit hashes.")
    parser.add_argument("config", nargs="?", default=root / "cmake/libraries.json", help="Path to libraries.json")
    args = parser.parse_args()

    config = Path(args.config)
    with config.open("r") as f:
        libraries = json.load(f)

    failed = False

    for lib_name, repo_url in libraries.items():
        lib_failed = False

        stat = check_github_actions_status(repo_url)
        if stat:
            ci_check = "✅"
        else:
            ci_check = "❌"
            lib_failed = True

        pinned_commit_hash = get_pinned_commit_hash(repo_url)
        if pinned_commit_hash is None:
            print(f"❌ {lib_name}: unable to parse pinned hash from libraries.json URL")
            pinned_display = "invalid"
            hash_check = "❌"
            remote_display = "unavailable"
            lib_failed = True
        else:
            pinned_display = pinned_commit_hash
            remote_commit_hash = get_last_commit_hash(repo_url)
            if remote_commit_hash is None:
                hash_check = "❌"
                remote_display = "unavailable"
                lib_failed = True
            else:
                remote_display = remote_commit_hash
                hash_match = remote_commit_hash == pinned_commit_hash
                if hash_match:
                    hash_check = "✅"
                else:
                    hash_check = "❌"
                    lib_failed = True

        print(
            f"{ci_check} CI {hash_check} HASH {lib_name}: "
            f"pinned={pinned_display} remote={remote_display}"
        )

        if lib_failed:
            failed = True

    if failed:
        print(f"\n❌ Some subprojects have failing CI or mismatched hashes")
        raise SystemExit(1)
    else:
        print("\n✅ All subprojects have passing CI and matching hashes")


if __name__ == "__main__":
    main()
