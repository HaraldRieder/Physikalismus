energies=[1,2;5,7]; # eV
amplitudes=[1,2;5,7];
f=0.5*[
  1,0,1,sqrt(2);
  sqrt(2),1,0,-1;
  -1,sqrt(2),1,0;
  0,-1,sqrt(2),-1
];
# check unitarity of f
f*f'

# Planck's constant over 2 pi in eVs units.
function retval = hbar()
  retval = 6.582119569e-16;
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
  ret = [0+0i,0+0i;0+0i,0+0i];
  f_adjugate = f_4_4';
  for k = 1:2
    for m = 1:2
      for i = 1:4
        for j = 1:4
          omega_ij = (E_4(i)-E_4(j))/hbar;
          exp_ij = exp(-i*omega_ij*t)
          for p = 1:2
            index_f = (k-1)*2+(p-1)+1;
            index_f_adj = (m-1)*2+(p-1)+1;
            ret(k,m) += d_4(i)*d_4(j)'*f_4_4(i,index_f)*f_adjugate(index_f_adj,j)*exp_ij;
          endfor
        endfor
      endfor
    endfor
  endfor
endfunction

for t = 0:1
  rho = rho_4_4(energies, amplitudes, f, t)  
endfor
