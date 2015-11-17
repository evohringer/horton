#!/usr/bin/env python

from horton import *
import numpy as np

# Load the coordinates of one water molecule from file.
# Use the XYZ file from HORTON's test data directory.
mol = Molecule.from_file('water.xyz')

# read in the charges from the second water molecule from file pointcharges.pc
ext_charges = Molecule.from_file('pointcharges.pc')

# Create a Gaussian basis set
obasis = get_gobasis(mol.coordinates, mol.numbers, 'cc-pVDZ')

# Create a linalg factory
lf = CholeskyLinalgFactory(obasis.nbasis)

# Compute Gaussian integrals
olp = obasis.compute_overlap(lf)
kin = obasis.compute_kinetic(lf)
na = obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, lf)
er = obasis.compute_electron_repulsion(lf)
# Calculate the interaction of the electrons with the external point charges
ext = obasis.compute_nuclear_attraction(ext_charges.coordinates, ext_charges.charges, lf)

# Create alpha orbitals
exp_alpha = lf.create_expansion()

# Initial guess
guess_core_hamiltonian(olp, kin, na, exp_alpha)


# Construct the restricted HF effective Hamiltonian
# The function compute_nucnuc_qmmm calculates the interaction between the nuclei
# and the interaction between the nuclei and the external charges
external = {'nn': compute_nucnuc_qmmm(mol.coordinates, mol.pseudo_numbers, ext_charges.coordinates, ext_charges.charges)}
terms = [
    RTwoIndexTerm(kin, 'kin'),
    RDirectTerm(er, 'hartree'),
    RExchangeTerm(er, 'x_hf'),
    RTwoIndexTerm(na, 'ne'),
    RTwoIndexTerm(ext, 'elec'),
]
ham = REffHam(terms, external)

# Decide how to occupy the orbitals (5 alpha electrons)
occ_model = AufbauOccModel(5)

# Converge WFN with plain SCF
occ_model.assign(exp_alpha)
dm_alpha = exp_alpha.to_dm()
scf_solver = ODASCFSolver(1e-6)
scf_solver(ham, lf, olp, occ_model, dm_alpha)

print 'energy',ham.cache['energy']
# Assign results to the molecule object and write it to a file, e.g. for
# later analysis
mol.title = 'RHF computation on water'
mol.energy = ham.cache['energy']

# The energy obtained from an ORCA calculation (same basis set) was: -76.034053013178
error = mol.energy - -76.034053013178


print 'The error is %+1.4e' % (error)
