"""
Tests that mdachecker can be imported.
"""

import sys
import mdachecker


def test_import_basic():
    assert "mdachecker" in sys.modules
