import numpy as np
import warnings

from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.lib.distances import self_capped_distance
from MDAnalysis.exceptions import NoDataError


class CheckStructure(AnalysisBase):
    r"""Simple structure checking analysis

    The trajectory is read, frame by frame, and the structure of the system is
    analysed for atomic overlaps and excessively long bonds.

    Parameters
    ----------
    atomgroup : AtomGroup or UpdatingAtomGroup
        Group of atoms to analyse the structure for. Please note that structure
        checking can be quite expensive for large groups of atoms as all-to-all
        distances need to be calculated. Adequately choosing a sub-section of
        the system can be quite useful.
    atomcut : float (optional)
        The minimum allowed distance between atoms in a system. Any atom pairs
        less than this distance away from each other will be considered as
        overlapping.
    bondcut : float (optional)
        The maximum allowed bond length between atoms in a system. Any bond
        greater than this value will be considered as a bond violation.
        Note: only bonds between atoms in the atomgroup will be analysed.


    Attributes
    ----------
    results.badframes : list
        A list containing an entry for each frame which violates the structure
        checking criteria. Each entry contains the following:
        ``[frame_id, overlap_info, bond_info]`` where ``overlap_info` and
        ``bond_info`` are 2-D ``dtype=object`` :class:`numpy.ndarray` with
        each row containing ``[atom1 index, atom2 index, distance]``. Please
        note that ``index`` here is defined at the 0-based internal MDAnalysis
        index. Should either ``atomcut`` or ``bondcut`` be set to ``None``,
        the respective results entry (i.e. ``overlap_info`` or ``bond_info``)
        will be set to ``None``.


    Raises
    ------
    NoDataError
        If no bonds are present in the input Universe.


    Notes
    -----
        * Bond analysis requires the presence of bonds. If no bonds are present
          in the Universe, a NoDataError will be thrown at class creation.
          Otherwise, should there be no bonds present in a given frame, then
          a warning will be thrown at that given frame. Please set bondcut
          to ``None`` if you do not wish to analyse bonds.
    """
    def __init__(self, atomgroup, atomcut=0.8, bondcut=2.5, **kwargs):
        super(CheckStructure, self).__init__(atomgroup.universe.trajectory,
                                             **kwargs)
        self.ag = atomgroup
        self.atomcut = atomcut

        if bondcut is not None:
            if not hasattr(self.ag.universe.atoms, "bonds"):
                errmsg = ("Universe has no bonds, cannot carry out "
                          "bond analysis, please set bondcut to None")
                raise NoDataError(errmsg)

        self.bondcut = bondcut

    def _prepare(self):
        self.results.badframes = []

    def _single_frame(self):

        overlaps = None
        longbonds = None

        if self.atomcut is not None:
            pairs, distances = self_capped_distance(
                    self.ag.atoms.positions, max_cutoff=self.atomcut,
                    box=self.ag.universe.dimensions, return_distances=True)
            if pairs.size > 0:
                # convert to ix indices
                ixs = []
                for pair in pairs:
                    ixs.append([self.ag.atoms[pair[0]].ix,
                                self.ag.atoms[pair[1]].ix])

                # concatenate the overlaps
                overlaps = np.concatenate(
                        (ixs, np.reshape(distances, (-1, 1))),
                        axis=1, dtype=object)

        if self.bondcut is not None:
            try:
                bond_ids = np.where(self.ag.intra_bonds.values() > self.bondcut)[0]
            except AttributeError:
                wmsg = f"atomgroup has no bonds for frame {self._frame_index}"
                warnings.warn(wmsg)

            if bond_ids.size > 0:
                # get the bonds
                bonds = self.ag.intra_bonds[bond_ids]
                ixs = bonds.indices
                vals = bonds.values()

                longbonds = np.concatenate(
                        (bonds.indices, np.reshape(bonds.values(), (-1, 1))),
                        axis=1, dtype=object)

        if (overlaps is not None) or (longbonds is not None):
            self.results.badframes.append((self._frame_index, overlaps, longbonds))
