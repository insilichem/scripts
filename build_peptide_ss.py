import random
from BuildStructure import placePeptide as _placePeptide

"""
Builds peptides with secondary structure assignments.

Sequence must be provided in 1-letter FASTA code.

Predicted secondary structures must be provided with 1-letter sequences:
    H -> alpha-helix
    I -> 3/10-helix
    P -> pi-helix
    A -> antiparallel beta-sheet (strand)
    S -> parallel beta-sheet (strand)

For example:

    sequence   = ERPKLFEELGKQIQQYMKIISSFKNKEDQYDHLDAADMTKVEKST
    prediction = CCHHHHHHHHHHHHHHHHHHHHHCCCCCCCCCCCHHHHHHHHHHH

If you use any other character (C, L, *, -, etc), it will be understood as a
free loop or coil. Since for coils there is no pre-assigned phi-psi angles,
we take a random guess out of a small library of coil angles obtained from
experimental structures. This can have a huge impact on the resulting
structure, so make sure to use your own values if needed and also, create
several candidates for further refinement.

This can be used from Chimera's commandline with:

    runscript build_peptide_ss.py <sequence> <prediction> [number of copies]
"""

_PHIPSI_ANGLES = {
    'H': (-57, -47),
    'I': (-49, -26),
    'P': (-57, -70),
    'A': (-139, 135),
    'S': (-119, 113),
}


# Randomly taken from a small subset. Feel free to add your own phi, psi angles
_COIL_ANGLES = [
    (-71.6011579949, 138.787346719),
    (-148.13178271, -8.52669209117),
    (-33.2816140035, -54.4979616768),
    (-104.48083984, 141.364437684),
    (-78.8640184608, 170.676947226),
    (-79.2485435476, 23.1742643632),
    (35.5995841589, 39.9507484959),
    (114.185411549, -15.9762813484),
    (-55.2126828065, 128.709122999),
    (-58.7044412711, 125.316529112),
    (-109.945735902, 142.471218644),
]


def place_peptide_ss(sequence, prediction, name='predicted', **kwargs):
    """
    Wrapper around Chimera's BuildStructure.placePeptide that computes
    the phi-psi angles on the fly following the predicted secondary
    structure assignaments.
    """
    angles = [_PHIPSI_ANGLES.get(p, random.choice(_COIL_ANGLES))
              for p in prediction]
    return _placePeptide(sequence, angles, name, **kwargs)


if __name__ == '__main__':
    import sys
    if len(sys.argv) not in (3,4):
        sys.exit('Usage: python build_peptide_ss.py <sequence> <prediction> [<n>]')

    sequence = sys.argv[1]
    prediction = sys.argv[2]

    n = 1 if len(sys.argv) == 3 else int(sys.argv[3])
    for _ in range(n):
        place_peptide_ss(sequence, prediction)
