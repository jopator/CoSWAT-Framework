import os, requests

url = "https://springernature.figshare.com/ndownloader/files/28919850"

# directory where this script is running
file_dir = os.path.dirname(os.path.abspath(__file__))
out = os.path.join(file_dir, "Bathymetry_Rasters.zip")

r = requests.get(url, stream=True)
r.raise_for_status()

with open(out, "wb") as f:
    for chunk in r.iter_content(8192):
        if chunk:
            f.write(chunk)