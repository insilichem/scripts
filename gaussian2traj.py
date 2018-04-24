#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import cclib as cc
import mdtraj as md


def gaussian_traj(path):
    g = cc.parser.Gaussian(path)
    parsed = g.parse()
    top = md.Topology()
    chain = top.add_chain()
    residue = top.add_residue('UNK', chain)
    for i, z in enumerate(parsed.atomnos):
        elem = md.element.Element.getByAtomicNumber(z)
        a = top.add_atom('A{}'.format(i+1), elem, residue, i)
    traj = md.Trajectory(parsed.atomcoords/10., top)
    traj = traj.center_coordinates()
    traj = traj.superpose(traj, frame=0)
    output_path = path + '.pdb'
    traj.save_pdb(output_path)
    print('Written {} to {}...'.format(traj, output_path))
    return traj

if __name__ == '__main__':
    import sys
    try:
        gaussian_traj(sys.argv[1])
    except IndexError:
        sys.exit('Usage: python gaussian2traj.py <gaussian_output_file>')