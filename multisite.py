#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
multisite.py
------------

This script can find potentially metal coordinating centers in
a protein structure.

It uses the following approach:
1. Flood the protein bounding box with a grid of points
2. For each point in the grid, evaluate if:
    a) There is a potentially coordinating residue within r A
    b) If so, if its beta-carbon its within r and 2*r
    c) If both criteria are satisfied, this point constitutes
       a "good" probe
3. Gather all good probes and cluster them by distance
4. For each cluster, compute the centroid and count the number
   of probes contained.
5. If the number of probes is greater than cutoff*max(num_probes),
   this centroid is considered a potentially coordinating center
   suitable for docking essays.
   Disclaimer: if the cluster is big enough, potentially coordinating
   residues might be further apart than the specified radius threshold,
   but since the reported point is meant to be origin of a search process,
   this should not matter.

Depending on the CLI arguments (run ./multisite.py -h for details), the script
will report the found sites on a plain-text table or run parallel GaudiMM
jobs for each site.

Atom types, names and residues can be overriden with the corresponding
environment variables, which should contain a list of blank-separated words:

    export COORDINATING_ATOM_TYPES="S3 O2 N2"
    export COORDINATING_ATOM_NAMES="NE2 ND1 NE1"
    export COORDINATING_RESIDUES="ASP GLU ASN GLN"

"""

from __future__ import print_function, absolute_import, division
import pychimera
pychimera.patch_environ()
pychimera.enable_chimera()
import sys
import os
import argparse
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

COORDINATING_ATOM_TYPES = set((os.environ.get('COORDINATING_ATOM_TYPES') or
                              'Npl O3 O2- S3 O2 N2').split())
COORDINATING_ATOM_NAMES = set((os.environ.get('COORDINATING_ATOM_NAMES') or
                               'NE2 ND1 OD1 OD2 OE1 OE2 O SG OG1 NE1 OH OXT').split())
COORDINATING_RESIDUES = set((os.environ.get('COORDINATING_RESIDUES') or
                             'ASP GLU ASN GLN HIS ARG CYS').split())
# COORDINATING_RESIDUES = set('ASP GLU'.split())


def find_metal_binding_sites(protein, tree=None, min_coordinators=2, radius=2.5, verbose=True,
                             backbone=True):
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
    grid = _grid(protein)
    for i, probe in enumerate(grid):
        residues = find_coordinating_residues(tree, probe, within=(radius, 2*radius),
                                              backbone=backbone)
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


def process_binding_sites(grid, clusters, grid_coordinators, residues=None,
                          cutoff=0.5, radius=3.5, clear=True, plot=True):
    base_radius = radius
    good_grid = grid[grid_coordinators.keys()]
    coordinators = np.array(grid_coordinators.values())
    if plot and clear:
        for model in chimera.openModels.list():
            if model.name == 'sphere':
                chimera.openModels.close([model])
    uniq, count = np.unique(clusters, return_counts=True)

    scale = count.max()
    centers, scores = [], []
    for i, size in zip(uniq, count):
        if size < cutoff*scale:
            continue
        cluster = good_grid[clusters==i]
        xyz = np.average(cluster, axis=0, weights=coordinators[clusters==i])
        radius = (pdist(cluster).max() / 2.0) + base_radius
        alpha = float(size)/(2*scale)
        if plot:
            chimera.runCommand('shape sphere mesh true '
                                            'center {xyz[0]},{xyz[1]},{xyz[2]} '
                                            'radius {} '
                                            'color 1,1,1,{}'.format(xyz, radius, alpha))
        centers.append(xyz)
        scores.append(size)

    if plot and residues is not None:
        for r in residues:
            for a in r.atoms:
                a.display = True
    return centers, scores


def _grid(protein, p=0.5, r=1.0):
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


def find_coordinating_residues(tree, site, within=(3.5, 7), criteria='restype',
                               strict_atom='CB', backbone=True):
    """
    Given a protein, search for potential coordinating residues
    within `radius` A from `site` coordinates. This is, residues
    featuring a beta-carbon within r and 2*r from site.
    """
    residues = set()
    within_top, within_bottom = max(within), min(within)
    nearby_atoms = tree.searchTree(site, within_bottom)
    site_point = chimera.Point(*site)
    for atom in nearby_atoms:
        if backbone and atom.name == 'O' and any([a.name == 'C' for a in atom.neighbors]):
            residues.add(atom.residue)
        elif atom.residue.type in COORDINATING_RESIDUES:
            if strict_atom:
                try:
                    cb_atom = atom.residue.atomsMap[strict_atom][0]
                except KeyError:
                    continue
                if within_top >= cb_atom.coord().distance(site_point) >= within_bottom:
                    residues.add(atom.residue)
            else:
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


def run(inputfile, n_processes=None, dry_run=False, cutoff=0.5,
        min_coordinators=2, radius=2.5, backbone=True, **kwargs):
    try:
        chimera.runCommand('open ' + inputfile)
        protein = chimera.openModels.list()[0]
        GAUDIMM_TPL = False
    except:
        cfg = Settings(inputfile, validation=False)
        protein = chimera.openModels.open(cfg.genes['Protein']['path'])[0]
        GAUDIMM_TPL = True

    print('Generating tree...')
    tree = AdaptiveTree(protein.atomCoordinatesArray().tolist(), protein.atoms, 1.0)
    print('Probing protein space...')
    sites, clusters, coordinators, residues = find_metal_binding_sites(protein, tree=tree,
        min_coordinators=min_coordinators, radius=radius, backbone=backbone, verbose=False)
    print('Post-processing', sites.shape[0], 'sites with cutoff', cutoff)
    centers, scores = process_binding_sites(sites, clusters, coordinators, residues,
                                            cutoff=cutoff, plot=False)
    rotamers = [find_coordinating_residues(tree, site, within=(radius, radius*2), strict_atom=None,
                backbone=True) for site in centers]

    if GAUDIMM_TPL:
        chimera.openModels.close([protein])
        templates = [prepare_input(cfg, site, rots) for (site, rots) in zip(sites, rotamers)]
        for i, template in enumerate(templates, 1):
                with open('template_{}.yaml'.format(i), 'w') as f:
                    f.write(template.toYAML())
        if not dry_run:
            _parallel_run(gaudi_run, templates, n_processes=n_processes)

        return

    lines = []
    center_width, score_width, residue_width = len('XYZ'), len('Probes'), len('Residues around centroid')
    pos_width = 1
    sorted_data = sorted(zip(centers, scores, rotamers), key=lambda e: e[1], reverse=True)
    for pos, (center, score, residues) in enumerate(sorted_data, 1):
        resnames = ','.join([str(r) for r in rotamers])
        pos_str, center_str, score_str  = str(pos), str(center), str(score)
        residue_str = ', '.join(['{}-{}'.format(r.type, r.id.position) for r in residues])
        if len(pos_str) > pos_width:
            pos_width = len(pos_str)
        if len(center_str) > center_width:
            center_width = len(center_str)
        if len(score_str) > score_width:
            score_width = len(score_str)
        if len(residue_str) > residue_width:
            residue_width = len(residue_str)
        lines.append((pos, center_str, score_str, residue_str))
    print(' {:>{pos_width}} | {:^{center_width}} | {:{score_width}} | {:{residue_width}}'.format(
            '#', 'XYZ', 'Probes', 'Residues around centroid', pos_width=pos_width, center_width=center_width,
            score_width=score_width, residue_width=residue_width))
    print('-{}-+-{}-+-{}-+-{}-'.format('-'*pos_width, '-'*center_width, '-'*score_width, '-'*residue_width))
    for line in lines:
        print(' {:>{pos_width}} | {:^{center_width}} | {:>{score_width}} | {:<{residue_width}}'.format(
                line[0], line[1], line[2], line[3], pos_width=pos_width,
                center_width=center_width, score_width=score_width, residue_width=residue_width))
    chimera.openModels.close([protein])


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


def parse_cli():
    p = argparse.ArgumentParser()
    p.add_argument('inputfile', type=str,
        help='Chimera-compatible Molecule or GaudiMM input file that should be used as a '
             'template to create a separate input file for each potential binding site found. '
             'If a  GaudiMM job is supplied, it expects a Molecule gene called `Protein`.')
    p.add_argument('-n', '--n_processes', default=None, type=int,
        help='How many cores should be used in the parallel calculation')
    p.add_argument('--dry_run', action='store_true',
        help="Don't run the GaudiMM jobs. Just write the patched input files.")
    p.add_argument('--backbone', action='store_true',
        help="Always consider backbone oxygen atoms as coordinators, no matter the residue.")
    p.add_argument('--cutoff', type=float, default=0.5,
        help='Cluster density cutoff. A binding site will only be reported if its '
             'cluster had more probes than cutoff*max(cluster_population). It should '
             'be a value within (0.0, 1.0]')
    p.add_argument('--min_coordinators', type=int, default=2,
        help='Minimum number of coordinating residues whose beta-carbon should be '
             'within radius a of a given grid probe for that probe to be considered '
             'potentially coordinating.')
    p.add_argument('--radius', type=float, default=2.5,
        help='Any potentially coordinating residue should have a beta-carbon within '
             'radius and 2*radius of a given grid probe.')
    return p.parse_args()


def main():
    args = parse_cli()
    run(**vars(args))


if __name__ == '__main__':
    main()
