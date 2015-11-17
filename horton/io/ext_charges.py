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
'''External charges File Format for QM/MM'''


import numpy as np

from horton.units import angstrom
from horton.periodic import periodic


__all__ = ['load_ext_charges']


def load_ext_charges(filename):
    ''' Load external charges from a .pc file.

       **Argument:**

       filename
            The file to load the charges from

       **Returns:** dictionary with ``coordinates`` of the charges and
                    the ``charges``.
    '''
    f = file(filename)
    size = int(f.next())
    coordinates = np.empty((size, 3), float)
    charges = np.empty(size, float)
    for i in xrange(size):
        words = f.next().split()
        charges[i] = float(words[0])
        coordinates[i,0] = float(words[1])*angstrom
        coordinates[i,1] = float(words[2])*angstrom
        coordinates[i,2] = float(words[3])*angstrom
    f.close()
    return {
        'coordinates': coordinates,
        'charges': charges
    }
