#!/usr/bin/env python
# -*- coding: utf-8 -*-

##############
# GaudiMM: Genetic Algorithms with Unrestricted
# Descriptors for Intuitive Molecular Modeling
#
# https://github.com/insilichem/gaudi
#
# Copyright 2017 Jaime Rodriguez-Guerra, Jean-Didier Marechal
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
##############

"""
Rotamers for metal coordination searches

- Preparse the protein space and look for potential metal sites
- For each possible site:
    - Calculate potential Rotamers
    - Modify the input template to include found rotamers and new search place
    - Run a separate GAUDI instance

- Collect all results in a single output
"""

from __future__ import print_function, absolute_import, division
import sys
from multiprocessing import Pool, cpu_count
from collections import OrderedDict
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform, pdist
import chimera
from CGLutil.AdaptiveTree import AdaptiveTree
from gaudi.box import open_models_and_close
from gaudi.cli.gaudi_run import main as gaudi_run
from gaudi.parse import Settings


COORDINATING_ATOM_TYPES = set('Npl O3 O2- S3 O2 N2'.split())
COORDINATING_ATOM_NAMES = set('NE2 ND1 OD1 OD2 OE1 OE2 O SG OG1 NE1 OH OXT'.split())
# COORDINATING_RESIDUES = set('ASP GLU ASN GLN HIS ARG CYS'.split())
COORDINATING_RESIDUES = set('ASP GLU'.split())


def find_metal_binding_sites(protein, tree=None, min_coordinators=2, radius=2.5, verbose=True):
    """
    Retrieve potential binding sites in a protein.

    Parameters
    ----------
    protein : chimera.Molecule
        The protein to scan for potential metal binding sites.

    Returns
    -------
    np.array
        A (n,3) array with the coordinates of the n sites found.

    Notes
    -----
    The algorithm could be implemented as:
        1. Fill the protein bounding box with probes
        2. For each probe, scan for potentially coordinating residues
        3. If a cluster of >3 probes is found and the ligand fits,
           the centroid of those can be considered a metal binding site.
    """
    good_probes = OrderedDict()
    good_residues = set()
    if tree is None:
        tree = AdaptiveTree(protein.atomCoordinatesArray().tolist(), protein.atoms, 1.0)
    grid = _grid(protein, p=0.5, r=1.0)
    for i, probe in enumerate(grid):
        residues = find_coordinating_residues(tree, probe, radius,)
        coordinating_res = [r for r in residues for a in r.atoms
                            if a.name in COORDINATING_ATOM_NAMES]
        coordinating_num = len(coordinating_res)
        if coordinating_num >= min_coordinators:
            good_probes[i] = coordinating_num
            good_residues.update(coordinating_res)
    if verbose:
        chimera.selection.setCurrent(good_residues)
        for res in good_residues:
            print(res)
    good_grid = grid[good_probes.keys()]
    distances = pdist(good_grid)
    linkaged = linkage(distances, method='average')
    flat_cluster = fcluster(linkaged, 10, criterion='distance')
    return grid, flat_cluster, good_probes, good_residues


def plot_binding_sites(grid, clusters, grid_coordinators, residues=None, radius=3.5, clear=True,
                       confidence=0.5):
    base_radius = radius
    good_grid = grid[grid_coordinators.keys()]
    coordinators = np.array(grid_coordinators.values())
    if clear:
        for model in chimera.openModels.list():
            if model.name == 'sphere':
                chimera.openModels.close([model])
    uniq, count = np.unique(clusters, return_counts=True)
    scale = count.max()
    for i, size in zip(uniq, count):
        if size < confidence*scale:
            continue
        cluster = good_grid[clusters==i]
        x, y, z = np.average(cluster, axis=0, weights=coordinators[clusters==i])
        radius = (pdist(cluster).max() / 2.0) + base_radius
        alpha = float(size)/(2*scale)
        chimera.runCommand('shape sphere mesh true '
                                        'center {},{},{} '
                                        'radius {} '
                                        'color 1,1,1,{}'.format(x, y, z, radius, alpha))
    if residues is not None:
        for r in residues:
            for a in r.atoms:
                a.display = True


def _grid(protein, p=0.1, r=0.3):
    bbox = protein.bbox()[1]
    (xmi, ymi, zmi), (xma, yma, zma) = bbox.llf.data(), bbox.urb.data()
    numx = int((xma - xmi + 2 * p) / r) + 1
    numy = int((yma - ymi + 2 * p) / r) + 1
    numz = int((zma - zmi + 2 * p) / r) + 1
    x = np.linspace(xmi - p, xma + p, numx)
    y = np.linspace(ymi - p, yma + p, numy)
    z = np.linspace(zmi - p, zma + p, numz)

    grid = []
    for i in xrange(len(x)):
        for j in xrange(len(y)):
            for k in xrange(len(z)):
                grid.append([x[i], y[j], z[k]])
    return np.reshape(grid, (-1, 3))


def find_coordinating_residues(tree, site, radius=2.5):
    """
    Given a protein, search for potential coordinating residues
    within `radius` A from `site` coordinates.
    """
    residues = set()
    sqradius = radius*radius
    nearby_atoms = tree.searchTree(site, radius)
    for atom in nearby_atoms:
        # if atom.idatmType in COORDINATING_ATOM_TYPES:
        if (#atom.name in COORDINATING_ATOM_NAMES and
             atom.residue.type in COORDINATING_RESIDUES):
            try:
                cb_atom = atom.residue.atomsMap['CB'][0]
            except KeyError:
                continue
            if 2.0*radius > cb_atom.coord().distance(chimera.Point(*site)) > radius:
                residues.add(atom.residue)
    return residues


def prepare_input(template, search, rotamers):
    """
    Given an options dict (template), modify the center parameter of
    a Search gene and the list of residues in Rotamers
    """
    cfg = template.copy()
    search_gene, rotamers_gene = None, None
    for gene in cfg.genes:
        if not search_gene and gene.module == 'gaudi.genes.search':
            search_gene = gene
        elif not rotamers_gene and gene.module == 'gaudi.genes.rotamers':
            rotamers_gene = gene

    search_gene.center = search
    rotamers_gene.residues = rotamers
    cfg.output.path += '_'.join([str(int(s)) for s in search])

    return cfg


def run(template, n_processes=None):
    cfg = Settings(template, validation=False)
    protein_path = cfg.genes['Protein']['path']

    with open_models_and_close(protein_path) as models:
        protein = models[0]
        tree = AdaptiveTree(protein.atomCoordinatesArray(), protein.atoms, 1.0)
        sites, _ = find_metal_binding_sites(protein)
        rotamers = [find_coordinating_residues(tree, site) for site in sites]

    templates = [prepare_input(cfg, site, rots) for (site, rots) in zip(sites, rotamers)]

    _parallel_run(gaudi_run, templates, n_processes=int(n_processes))


def _parallel_run(func, args, n_processes=None, verbose=True):
    if n_processes > cpu_count() or n_processes < 1:
        n_processes = cpu_count()
    pool = Pool(processes=n_processes, maxtasksperchild=1)
    if verbose:
        print("Running with {} processes...".format(n_processes))
    try:
        pool.map(func, args, chunksize=1)
    except KeyboardInterrupt:
        pool.terminate()
    except Exception as e:
        print(e)
        pool.terminate()
    finally:
        pool.close()
        pool.join()


def main():
    run(*sys.argv[1:])


if __name__ == '__main__':
    main()
