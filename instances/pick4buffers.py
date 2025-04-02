import os
import shutil

in_folder = "generated-small"
out_folder = "4-buffers"

files = os.listdir(in_folder)  # Lists all files and folders
files = [f for f in files if os.path.isfile(os.path.join(in_folder, f))]  # Filters only files

for file in files:
    with open(in_folder + "/" + file, 'r') as f:
        if f.readline().startswith("4 "):
            shutil.copy(in_folder + "/" + file, out_folder)
