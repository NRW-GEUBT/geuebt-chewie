#!/usr/bin/env bash
set -Eeu

# Examople of post-deploy script downloading a database for local use


# URL of the .tar archive
url="<URL>"

# Local directory to save the file
local_dir="<PATH>"

# Name of the downloaded file
file_name="<NAME>"

# Remote file (includes .tar.gz extension)
remote_name=$(basename "$url")

# Create local directory if it doesn't exist
mkdir -p "$local_dir"

# Check if the file or directory already exists locally
if [[ -e "$local_dir/$file_name" ]]; then

  echo "The file already exists, skipping download"

else

  # Download the file and save it to the local directory
  echo "Downloading file..."
  curl -L -o "$local_dir/$remote_name" "$url"

  # Calculate the hash of the downloaded file
  download_hash=$(openssl dgst -r -sha256 "$local_dir/$remote_name")
  
  echo "$download_hash" > "$local_dir/$file_name.sha256"
  date --iso-8601='minutes' > "$local_dir/$file_name.timestamp"
  echo "$url" > "$local_dir/$file_name.source"
  # No directory in archive for kraken!
  mkdir -p "$local_dir/$file_name"
  tar -xzvf "$local_dir/$remote_name" -C "$local_dir/$file_name" && rm "$local_dir/$remote_name"

fi
