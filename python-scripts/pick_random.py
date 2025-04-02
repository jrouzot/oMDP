import os
import shutil
import random

def copy_random_files(source_folder, destination_folder, n):
    # Create the destination folder if it doesn't exist
    if not os.path.exists(destination_folder):
        os.makedirs(destination_folder)

    # Get a list of all files in the source folder
    all_files = os.listdir(source_folder)

    # Choose n random files
    random_files = random.sample(all_files, n)

    # Copy each random file to the destination folder
    for file_name in random_files:
        source_path = os.path.join(source_folder, file_name)
        destination_path = os.path.join(destination_folder, file_name)
        shutil.copy2(source_path, destination_path)
        print(f"Copied {file_name} to {destination_folder}")

# Example usage
source_folder = '../instances/expes-cp-2025'
destination_folder = '../instances/expes-test'
number_of_files_to_copy = 100

copy_random_files(source_folder, destination_folder, number_of_files_to_copy)

