"""
Basic tests for StructureChecker
"""
from importlib import resources

import pytest
import MDAnalysis as mda
from MDAnalysis.exceptions import NoDataError

from mdachecker.checkers import CheckStructure


@pytest.fixture(scope='module')
def nolig_u():
    with resources.path('mdachecker.tests.datafiles',
                        '2RBN.nolig.pdb') as fn:
        u = mda.Universe(str(fn))
    return u


@pytest.fixture(scope='module')
def lig_u():
    with resources.path('mdachecker.tests.datafiles',
                        '2RBN.pdb') as fn:
        u = mda.Universe(str(fn))
    return u


@pytest.fixture(scope='module')
def longbond_u():
    with resources.path('mdachecker.tests.datafiles',
                        '2RBN.nolig.pdb') as fn:
        u = mda.Universe(str(fn))

    u.add_TopologyAttr('bonds', [(0, 1648),])
    return u


def test_nobonds_error(nolig_u):
    with pytest.raises(NoDataError, match="Universe has no bonds"):
        CheckStructure(nolig_u.atoms).run()


def test_nobonds_warning(lig_u):
    ag = lig_u.select_atoms('protein')

    with pytest.warns(UserWarning, match="atomgroup has no bonds"):
        CheckStructure(ag).run()


def test_altloc_overlaps(lig_u):

    checker = CheckStructure(lig_u.atoms, bondcut=None).run()

    assert len(checker.results.badframes) == 1
    assert checker.results.badframes[0][0] == 0
    assert checker.results.badframes[0][2] is None
    assert len(checker.results.badframes[0][1]) == 27
    for pairs in checker.results.badframes[0][1]:
        at1_loc = lig_u.atoms[pairs[0]].altLoc
        at2_loc = lig_u.atoms[pairs[1]].altLoc
        assert (at1_loc == 'B') or (at2_loc == 'B')


def test_longbond(longbond_u):

    checker = CheckStructure(longbond_u.atoms, atomcut=None).run()

    assert len(checker.results.badframes) == 1
    assert checker.results.badframes[0][0] == 0
    assert checker.results.badframes[0][1] is None
    assert len(checker.results.badframes[0][2]) == 1
    dist = checker.results.badframes[0][2][0][2]
    assert pytest.approx(20.49389, 0.01) == dist
