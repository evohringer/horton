# -*- coding: utf-8 -*-
# Horton is a development platform for electronic structure methods.
# Copyright (C) 2011-2013 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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
#pylint: skip-file


import numpy as np

from horton import *
from horton.test.common import tmpdir


def test_load_ext_charges():
    fn = context.get_fn('test/water_ext_charges.pc')
    mol = IOData.from_file(fn)
    check_water_ext_charges(mol)


def check_water_ext_charges(mol):
    assert mol.charges[0] == -0.82
    assert mol.charges[1] == 0.41
    assert mol.charges[1] == 0.41
    assert abs(np.linalg.norm(mol.coordinates[0] - mol.coordinates[1])/angstrom - 0.944390) < 1e-5
    assert abs(np.linalg.norm(mol.coordinates[0] - mol.coordinates[2])/angstrom - 0.947998) < 1e-5
