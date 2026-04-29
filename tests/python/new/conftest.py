"""Pytest configuration for new interface tests.

Adds the parent tests/python/ directory to sys.path so that utils.py
(shared with the old interface tests) can be imported.
"""

import os
import sys

# Add tests/python/ to sys.path for shared utilities (utils.py)
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
