import os
import requests
import zipfile
import tarfile
import subprocess

# URLs to download
dependencies = {
    "SDPNALv1.0" : "https://blog.nus.edu.sg/mattohkc/files/2023/11/SDPNALv1.0.zip",
    "cliquer-1.21" : "http://users.aalto.fi/~pat/cliquer/cliquer-1.21.tar.gz"
}

# Directory to create
target = "ThirdParty"

# Ensure the ThirdParty directory exists
os.makedirs(target, exist_ok=True)

# Function to download a file
def download_file(url, dest_folder):
    local_filename = os.path.join(dest_folder, url.split('/')[-1])
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(local_filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
    return local_filename

# Function to extract an archive
def extract_file(filepath, dest_folder):
    if filepath.endswith(".zip"):
        with zipfile.ZipFile(filepath, 'r') as zip_ref:
            zip_ref.extractall(dest_folder)
    elif filepath.endswith(".tar.gz") or filepath.endswith(".tgz"):
        with tarfile.open(filepath, 'r:gz') as tar_ref:
            tar_ref.extractall(dest_folder)
    else:
        raise ValueError(f"Unsupported file type: {filepath}")

# Download and extract files
for dep in dependencies:
    print("> Downloading %s ..." % dep, end="")
    archive_path = download_file(dependencies[dep], target)
    print("Done!")
    print("> Exctracting %s ..." % dep, end="")
    extract_file(archive_path, target)
    os.remove(archive_path)
    print("Done!")

# Move to target directory
target = os.path.join(target, 'cliquer-1.21')
print("> Moving to %s" % target)
os.chdir(target)

# Run "make cl" command
try:
    print("> Building cliquer-1.21 ...", flush=True)
    subprocess.run(["make", "cl"], check=True)
    print("> Done!")
except subprocess.CalledProcessError as e:
    print(f"> Error running building Cliquer: {e}")

print("> Installation completed successfully.")
