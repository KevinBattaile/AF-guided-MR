#!/usr/bin/env python3
import os
import sys

# Determine the absolute path of this script
current_dir = os.path.dirname(os.path.abspath(__file__))
# Add the src/ directory to sys.path
src_dir = os.path.join(current_dir, "src")
if src_dir not in sys.path:
    sys.path.insert(0, src_dir)

# Import the CLI function and execute it
from af_guided_mr.cli import main

if __name__ == "__main__":
    main()
