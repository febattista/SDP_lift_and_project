# Copyright (C) 2024 Federico Battista, Fabrizio Rossi, Stefano Smriglio
# SPDX-License-Identifier: GPL-3.0-or-later
# Downloads and builds third-party dependencies: SDPNAL+ and cliquer.

import os
import shutil
import subprocess
import tarfile
import zipfile

dependencies = {
    "SDPNALv1.0": "https://blog.nus.edu.sg/mattohkc/files/2023/11/SDPNALv1.0.zip",
    "cliquer-1.21": "http://users.aalto.fi/~pat/cliquer/cliquer-1.21.tar.gz",
}

target = "ThirdParty"
os.makedirs(target, exist_ok=True)


def download_file(url, dest_folder):
    filename = os.path.join(dest_folder, url.split("/")[-1])

    if shutil.which("curl"):
        cmd = [
            "curl",
            "-L",              # follow redirects
            "--fail",          # fail on HTTP errors
            "--retry", "5",    # retry transient failures
            "--retry-delay", "2",
            "-o", filename,
            url,
        ]
    elif shutil.which("wget"):
        cmd = [
            "wget",
            "-O", filename,
            "--tries=5",
            url,
        ]
    else:
        raise RuntimeError("Neither curl nor wget is available on this system.")

    subprocess.run(cmd, check=True)
    return filename


def extract_file(filepath, dest_folder):
    if filepath.endswith(".zip"):
        with zipfile.ZipFile(filepath, "r") as zip_ref:
            zip_ref.extractall(dest_folder)
    elif filepath.endswith(".tar.gz") or filepath.endswith(".tgz"):
        with tarfile.open(filepath, "r:gz") as tar_ref:
            tar_ref.extractall(dest_folder)
    else:
        raise ValueError(f"Unsupported file type: {filepath}")


for dep, url in dependencies.items():
    print(f"> Downloading {dep} ...", end="", flush=True)
    archive_path = download_file(url, target)
    print("Done!")

    print(f"> Extracting {dep} ...", end="", flush=True)
    extract_file(archive_path, target)
    os.remove(archive_path)
    print("Done!")

target = os.path.join(target, "cliquer-1.21")
print(f"> Moving to {target}")
os.chdir(target)

try:
    print("> Building cliquer-1.21 ...", flush=True)
    subprocess.run(["make", "cl"], check=True)
    print("> Done!")
except Exception:
    print("> Warning: Error building Cliquer: LP formulation will not work!")

print("> Installation completed successfully.")