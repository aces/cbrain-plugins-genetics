import subprocess
import os
import argparse
import sys
import datetime

def print_log(message):
    """Prints a log message with a timestamp."""
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    subprocess.run(["echo", f"{timestamp}: {message}"], shell=False)

def print_header(tool_name, tool_version, wrapper_version, wrapper_author):
    """Prints a header with version and author information."""
    print_log(f"Starting {tool_name} {tool_version} wrapper script")
    print_log(f"By: {wrapper_author}")
    print_log(f"Version: {wrapper_version}")

def is_dir(path):
    """Checks if a path is a directory for argparse."""
    if os.path.isdir(path):
        return path
    else:
        raise argparse.ArgumentTypeError(f"{path} is not a valid path")

def is_file(path):
    """Checks if a path is a file for argparse."""
    if os.path.isfile(path):
        return path
    else:
        raise argparse.ArgumentTypeError(f"{path} is not a valid file")

def run_command(command, shell=True):
    """Runs a shell command, logs it, and handles errors."""
    print_log(f"Executing command: {' '.join(command) if isinstance(command, list) else command}")
    try:
        # using check=True to raise CalledProcessError on non-zero exit status
        result = subprocess.run(command, shell=shell, check=True, text=True)
        return result
    except subprocess.CalledProcessError as e:
        print_log(f"Error executing command: {command}")
        print_log(f"Return code: {e.returncode}")
        print_log(f"Error message: {e.stderr}")
        print_log(f"Output: {e.stdout}")
        sys.exit(e.returncode)

