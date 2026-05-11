#!/usr/bin/env python3
# chemeleon_fingerprint.py
#
# this file contains the class CheMeleonFingerprint which can be instantiated
# and called to generate the CheMeleon learned embeddings for a list of SMILES
# strings and/or RDKit Mols. you may wish to simply copy or download this file directly for use,
# or adapt the code for your own purposes. No other files are required for it
# to work, though you must `pip install 'chemprop>=2.2.0'` for this to run.
#
# run `python chemeleon_fingerprint.py` for a quick usage demo, otherwise you
# should `import` the CheMeleonFingerprint class into your other code and use
# it there (following the example at the bottom of this file) to generate
# your learned fingerprints
#
# MIT License
#
# Copyright (c) 2025 Jackson Burns
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

from pathlib import Path
from urllib.request import urlretrieve

import sys
import numpy as np
import torch # type: ignore
from chemprop import featurizers, nn # type: ignore
from chemprop.data import BatchMolGraph # type: ignore
from chemprop.models import MPNN # type: ignore
from chemprop.nn import RegressionFFN # type: ignore
from rdkit.Chem import Mol, MolFromSmiles

class CheMeleonFingerprint:
    def __init__(self, device: str | torch.device | None = None):
        self.featurizer = featurizers.SimpleMoleculeMolGraphFeaturizer()
        agg = nn.MeanAggregation()
        ckpt_dir = Path().home() / ".chemprop"
        ckpt_dir.mkdir(exist_ok=True)
        mp_path = ckpt_dir / "chemeleon_mp.pt"
        if not mp_path.exists():
            print('INFO: Downloading chemeleon_mp.pt to %s...' % mp_path, file=sys.stderr)
            urlretrieve(
                r"https://zenodo.org/records/15460715/files/chemeleon_mp.pt",
                mp_path,
            )
        chemeleon_mp = torch.load(mp_path, weights_only=True)
        mp = nn.BondMessagePassing(**chemeleon_mp["hyper_parameters"])
        mp.load_state_dict(chemeleon_mp["state_dict"])
        self.model = MPNN(
            message_passing=mp,
            agg=agg,
            predictor=RegressionFFN(input_dim=mp.output_dim),  # not actually used
        )
        self.model.eval()
        if device is not None:
            self.model.to(device=device)

    # one 2048 floats ndarray per molecule
    def __call__(self, mol: Mol) -> np.ndarray:
        bmg = BatchMolGraph([self.featurizer(mol)])
        bmg.to(device=self.model.device)
        with torch.no_grad():
            fps = self.model.fingerprint(bmg).numpy(force=True)
            return fps[0]

# simple test one one molecule
if __name__ == "__main__":
    chemeleon_fingerprint = CheMeleonFingerprint()
    fp = chemeleon_fingerprint(MolFromSmiles("CCC"))
    assert(len(fp) == 2048)
    for i, x in enumerate(fp):
        if i > 0:
            print('\t%f' % x, end='')
        else:
            print('%f' % x, end='')
    print() # EOL
