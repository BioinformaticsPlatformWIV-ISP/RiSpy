import os
import shutil

def check_exists(path: str, label: str):
    if not os.path.exists(path):
        raise FileNotFoundError(f"{label} does not exist: {path}")

def check_is_dir(path: str, label: str):
    if not os.path.isdir(path):
        raise NotADirectoryError(f"{label} is not a directory: {path}")

def check_tool(tool: str):
    if shutil.which(tool) is None:
        raise RuntimeError(f"Required tool not found in PATH: {tool}")
