from MDAnalysis.analysis.base import AnalysisBase


class CheckStructure(AnalysisBase):
    def __init__(self, atomgroup, atomcut=0.8, bondcut=2.5, **kwargs):
        super(CheckStructure, self).__init__(atomgroup.universe.trajectory,
                                             **kwargs)

        self.ag = atomgroup
        self.atomcut = atomcut
        self.bondcut = bondcut

    def _prepare(self):
        self.badframes = []

    def _single_frame(self):
        self.ag.
