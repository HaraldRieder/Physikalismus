# Planck's constant over 2 pi in eVs units.
function ret = h_bar_eV()
  ret = 6.582119569e-16;
endfunction

# Calculates 2x2 density matrix of state |psi> in a 2 qubit space
# depending on the view defined by unitary matrix f and the time.
# f is meant relative to the energy representation.
#
# E_4 vector with 4 energy eigenvalues
# d_4 amplitudes of energy eigenvectors
# (|psi> = sum of d_i * |E_i>)
# f_4_4 4x4 transformation matrix to 2x2 split view
# t time scalar value
function ret = rho_4_4(E_4, d_4, f_4_4, t)  
  omega_ij = 
endfunction
