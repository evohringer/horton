#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
# Copyright (C) 2011-2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>
#
# This file is part of Horton.
#
# Horton is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# Horton is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
#--


import sys, argparse, os

import h5py as h5, numpy as np
from horton import System, angstrom, setup_weights, ESPCost, log
from horton.scripts.common import reduce_data
from horton.scripts.espfit import parse_wdens, parse_wnear, parse_wfar, load_rho, save_weights


def parse_args():
    parser = argparse.ArgumentParser(prog='horton-esp-cost.py',
        description='Construct a cost function and fit charges to ESP.')

    parser.add_argument('cube',
        help='The cube file.')
    parser.add_argument('--reduce', '-r', default=1, type=int,
        help='Reduce the grid by subsamping with the given stride in all three '
             'directions. Zero and negative values are ignored.')
    parser.add_argument('--overwrite', default=False, action='store_true',
        help='Overwrite existing output in the HDF5 file')
    parser.add_argument('--suffix', default=None, type=str,
        help='Add an additional suffix to the HDF5 output group.')

    parser.add_argument('--rcut', default=10.0, type=float,
        help='The real-space cutoff for the electrostatic interactions in '
             'angstrom.')
    parser.add_argument('--alpha-scale', default=3.0, type=float,
        help='The alpha scale (alpha = alpha_scale/rcut) for the separation '
             'between short-range and long-range electrostatic interactions.')
    parser.add_argument('--gcut-scale', default=1.1, type=float,
        help='The gcut scale (gcut = gcut_scale*alpha) for the reciprocal '
             'space constribution to the electrostatic interactions.')

    parser.add_argument('--wdens', default=None, type=str,
        help='Define weights based on an electron density. The argument has '
             'the following format: "dens.cube:rho0:alpha". The last or the '
             'last two fields are optional and have as default values 2e-4 and '
             '1.0 respectively. The density is loaded from the file in the '
             'first field. The second field is the density at which the weight '
             'function switches and the third field is the half width of the '
             'switching function in log10(bohr**-3) units.')
    parser.add_argument('--wnear', default=None, type=str, nargs='+',
        help='Define weights that go to zero near the nuclei. Multiple '
             'arguments are allowed and have the following format: '
             '"number:r0:gamma". The first field is the elemenet number and '
             'the last two fields define the switching function and are in '
             'Angstrom. The last field is optional and is 0.5*A by default. '
             'The second field is the middle of the switching function. '
             'The third field is the half width of the switching function.')
    parser.add_argument('--wfar', default=None, type=str,
        help='Define weights that go to zero far from the nuclei. The '
             'argument has the following format: "r0:gamma". '
             'The first field is the distance (to the closest nucleus) at '
             'which the weight switches to zero. The second field is the half '
             'width of the switching function. To avoid artifacts, a smooth '
             'version of "distance to the closest nucleus" is used. The second '
             'field is optional and has 1.0 as default value')
    parser.add_argument('--wsave', default=None, type=str,
        help='Save the weights array to the given cube file.')

    # TODO: add argument to chop of last slice(s)

    return parser.parse_args()


def main():
    args = parse_args()

    # check if the folder is already present in the output file
    fn_h5 = args.cube + '.h5'
    grp_name = 'espfit_r%i' % args.reduce
    if args.suffix is not None:
        grp_name += '_'+args.suffix

    if os.path.isfile(fn_h5):
        with h5.File(fn_h5, 'r') as f:
            if 'espfit/%s' % grp_name in f and not args.overwrite:
                if log.do_warning:
                    log.warn('Skipping because "%s" is already present in the output.' % grp_name)
                return

    # Load the system
    if log.do_medium:
        log('Loading potential array')
    sys = System.from_file(args.cube)
    ui_grid = sys.props['ui_grid']
    esp = sys.props['cube_data']

    # Reduce the grid if required
    if args.reduce > 1:
        esp, ui_grid = reduce_data(esp, ui_grid, args.reduce)

    # Construct the weights for the ESP Cost function.
    wdens = parse_wdens(args.wdens)
    if wdens is not None:
        if log.do_medium:
            log('Loading density array')
        rho = load_rho(wdens[0], args.reduce, ui_grid)
        wdens = (rho,) + wdens[1:]
    weights = setup_weights(sys, ui_grid,
        dens=wdens,
        near=parse_wnear(args.wnear),
        far=parse_wnear(args.wfar),
    )
    # write the weights to a cube file if requested
    if args.wsave is not None:
        if log.do_medium:
            log('Saving weights array')
        save_weights(args.wsave, sys, ui_grid, weights)
    # rescale weights such that the cost function is the mean-square-error
    if weights.max() == 0.0:
        raise ValueError('No points with a non-zero weight were found')
    wmax = weights.min()
    wmin = weights.max()
    weights /= weights.sum()

    # Construct the cost function
    if log.do_medium:
        log('Setting up cost function (may take a while)')
    cost = ESPCost.from_grid_data(
        sys, ui_grid, esp, weights, args.rcut*angstrom, args.alpha_scale,
        args.gcut_scale)

    # Store command line arguments
    results = {}
    results['reduce'] = args.reduce
    results['rcut'] = args.rcut
    results['alpha_scale'] = args.alpha_scale
    results['gcut_scale'] = args.gcut_scale

    # Store cost function info
    results['A'] = cost._A
    results['B'] = cost._B
    results['C'] = cost._C

    # Store cost function properties
    results['evals'] = np.linalg.eigvalsh(cost._A)
    abs_evals = abs(results['evals'])
    if abs_evals.min() == 0.0:
        results['cn'] = 0.0
    else:
        results['cn'] = abs_evals.max()/abs_evals.min()

    # Report some on-screen info
    if log.do_medium:
        log('Important parameters:')
        log.hline()
        log('Number of grid points:        %12i' % esp.size)
        log('Used number of grid points:   %12i' % (weights>0).sum())
        log('Lowest weight:                %12.5e' % wmin)
        log('Highest weight:               %12.5e' % wmax)
        log('Lowest abs eigen value:       %12.5e' % abs_evals.min())
        log('Highest abs eigen value:      %12.5e' % abs_evals.max())
        log('Condition number:             %12.5e' % results['cn'])
        log.hline()

    # Store the results in an HDF5 file
    with h5.File(fn_h5) as f:
        # Store essential system info
        # TODO: this should be implemented with an improved implementation of System.to_file
        sys_grp = f.require_group('system')
        coordinates = sys_grp.require_dataset('coordinates', sys.coordinates.shape, float, exact=True)
        coordinates[:] = sys.coordinates
        numbers = sys_grp.require_dataset('numbers', sys.numbers.shape, long, exact=True)
        numbers[:] = sys.numbers
        rvecs = sys_grp.require_dataset('rvecs', (sys.cell.nvec, 3), float, exact=True)
        rvecs[:] = sys.cell.rvecs

        # Store grid rvecs
        grp_espfit = f.require_group('espfit')
        if grp_name in grp_espfit:
            del grp_espfit[grp_name]
        grp = grp_espfit.create_group(grp_name)
        grp['grid_rvecs'] = ui_grid.grid_cell.rvecs

        # Store results
        for key, value in results.iteritems():
            grp[key] = value

        if log.do_medium:
            log('Results written to %s:espfit/%s' % (fn_h5, grp_name))


if __name__ == '__main__':
    main()