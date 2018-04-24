#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import pychimera
pychimera.patch_environ()
pychimera.enable_chimera()


import os
import zipfile
import tempfile
import chimera
from glob import glob
from Measure.center import atoms_center_of_mass as center_of_mass

tempdir = tempfile.mkdtemp('gaudiview')


def parse_zip(path):
    """
    GAUDI zips its results files. We extract them on the fly to temp
    directories and feed those to the corresponding parsers (Chimera
    and YAML, as of now).
    """
    index = -1
    with zipfile.ZipFile(path) as z:
        index += 1
        tmp = os.path.join(tempdir, os.path.splitext(os.path.basename(path))[0])
        try:
            os.mkdir(tmp)
        except OSError:  # Assume it exists
            pass
        z.extractall(tmp)
        molecules = sorted(entry for entry in os.listdir(tmp)
                           if entry.endswith('.mol2'))
        for molecule in molecules:
            absname = os.path.join(tmp, molecule)
            subid = 0
            for m in chimera.openModels.open(absname, baseId=index, temporary=True,
                                             subid=subid, shareXform=True):
                yield m
                subid += 1


def render(counter=0, output_dir='.'):
    # ligand = sorted(chimera.openModels.list(), key=lambda m: m.numAtoms)[0]
    # chimera.runCommand('show #{}.{} zr < 5'.format(ligand.id, ligand.subid))
    chimera.runCommand('show @V zr < 5')
    chimera.runCommand('copy file {}/img{:03d}.png png width 1024'.format(output_dir, counter))


def main(mask='*.zip', output_directory='.'):
    try:
        os.mkdir(output_directory)
    except:
        pass
    tan_color = chimera.colorTable.getColorByName('tan')
    files = sorted(glob(mask), reverse=True)
    # distances = []
    # print('Presorting...')
    # for i, zip_file in enumerate(files):
    #     centers = []
    #     for molecule in parse_zip(zip_file):
    #         centers.append(center_of_mass(molecule.atoms))
    #     distances.append(chimera.Point(*centers[0]).distance(chimera.Point(*centers[1])))
    #     chimera.closeSession()
    
    # distance_sorted = [f for (d, f) in sorted(zip(distances, files), reverse=True)]
    for i, zip_file in enumerate(files):
        print('{}/{} Rendering {}...'.format(i+1, len(files), zip_file))
        molecules = list(parse_zip(zip_file))
        protein = sorted(molecules, key=lambda m: m.numAtoms)[-1]
        protein.color = tan_color
        if not i:
            chimera.runCommand('focus')
        render(counter=i, output_dir=output_directory)
        chimera.closeSession()


if __name__ == '__main__':
    import sys
    main(*sys.argv[1:4])
    