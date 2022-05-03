"""
Tests that propkatraj can be imported.
"""

import sys
import mdachecker


def test_import_basic():
    assert "propkatraj" in sys.modules
