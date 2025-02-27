#!/usr/bin/env python3

import requests
import argparse
import os
import re
import json
from pprint import pprint
from typing import NamedTuple, Optional, TextIO

# GitHub repository details
REPO = "kyclark/sufr"
API_URL = f"https://api.github.com/repos/{REPO}/releases"


class ReleaseInfo(NamedTuple):
    """ Release info """
    os: str
    arch: str


class Args(NamedTuple):
    """ Command-line args """
    version: str
    json: Optional[TextIO]


# --------------------------------------------------
def get_args() -> Args:
    parser = argparse.ArgumentParser(
        description="Update a specific release in GitHub."
    )

    parser.add_argument(
        "version", help="The tag name of the release to update"
    )

    parser.add_argument(
        "-j",
        "--json",
        type=argparse.FileType("rt"),
        help="Local JSON file"
    )

    args = parser.parse_args()

    return Args(version=args.version, json=args.json)


# --------------------------------------------------
def main() -> None:
    args = get_args()
    releases = get_releases_data(args.json)

    if release := find_release_by_version(releases, args.version):
        print("Release '{}' has {} assets".format(args.version, len(release["assets"])))
        markdown_table = generate_markdown_table(release)
        update_release_body(release["id"], markdown_table)
        print(f"Release '{args.version}' updated successfully.")
    else:
        print(f"Release version '{args.version}' not found.")


# --------------------------------------------------
def get_releases_data(file: Optional[TextIO]):
    if file:
        return json.loads(file.read())

    headers = {
        "Authorization": f'token {os.getenv("GITHUB_TOKEN")}',
        "Accept": "application/vnd.github.v3+json",
    }
    response = requests.get(API_URL, headers=headers)
    response.raise_for_status()
    return response.json()


# --------------------------------------------------
def find_release_by_version(releases, version):
    for release in releases:
        if release["tag_name"] == version:
            return release
    return None


# --------------------------------------------------
def extract_os_arch_from_filename(filename) -> Optional[ReleaseInfo]:
    pattern = re.compile(r"^(.+)(?:\.tar\.gz|\.zip)$")

    if match := pattern.search(filename):
        stem = match.group(1)
        _, arch, os = stem.split("-", 2)

        if os == "macos-latest" or os == "apple-darwin":
            os = "MacOS"
        elif os == "ubuntu-latest":
            os = "Ubuntu"
        elif re.search("linux", os):
            os = "Linux"
        elif re.search("windows", os):
            os = "Windows"

        if arch == "x64" or arch == "x86_64":
            arch = "Intel/AMD 64-bit"
        elif arch == "386":
            arch = "Intel/AMD 32-bit"
        elif arch == "arm64":
            if os == "MacOS":
                arch = "M1/M2/M3 (ARM 64-bit)"
            else:
                arch = "ARM 64-bit"
        elif arch == "arm":
            arch = "ARM 32-bit"

        return ReleaseInfo(os, arch)


# --------------------------------------------------
def generate_markdown_table(release) -> str:
    table = "### Release Assets\n"
    table += "| OS | Architecture | Link |\n"
    table += "|---------|----------|-------------|\n"

    for asset in release["assets"]:
        if info := extract_os_arch_from_filename(asset["name"]):
            print(">>> Asset {}".format(asset["name"]))
            download_url = asset["browser_download_url"]
            table += f"| {info.os}  | {info.arch}  | [Download]({download_url}) |\n"

    # Add note about Mac binary signing restriction
    table += (
        "\n***(For Macs) To address the Mac binary signing restriction, use the following command:\n"
        "\n```\n"
        "sudo xattr -dr com.apple.quarantine <path to file>/my-binary-amd64\n"
        "```\n"
    )

    return table


# --------------------------------------------------
def update_release_body(release_id, new_body):
    headers = {
        "Authorization": f'token {os.getenv("GITHUB_TOKEN")}',
        "Accept": "application/vnd.github.v3+json",
    }
    update_url = f"https://api.github.com/repos/{REPO}/releases/{release_id}"
    data = {"body": new_body}

    response = requests.patch(update_url, headers=headers, json=data)
    response.raise_for_status()


# --------------------------------------------------
if __name__ == "__main__":
    main()
