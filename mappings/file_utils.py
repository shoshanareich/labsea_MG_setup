import glob
import os

def get_files(run_dir, file_name):

    stems = set()  # use a set to avoid duplicates from .meta/.data pairs
    
    for f in glob.glob(os.path.join(run_dir, file_name + "_*.*.*")):
        name = os.path.basename(f)
    
        # skip "_lo." files
        if "_lo." in name:
            continue
    
        parts = name.split(".")
        if (
            len(parts) == 3 and
            parts[1].isdigit() and len(parts[1]) == 10 and
            parts[2] in ("meta", "data")
        ):
            stem = ".".join(parts[:2])  # e.g. xx_theta.0000000001
            stems.add(stem)
    
    return stems
