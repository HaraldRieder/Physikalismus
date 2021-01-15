energies=[1,2;5,7]; # eV
amplitudes=[1,2;5,7]/sqrt(79);
f=0.5*[
  1,0,1,sqrt(2);
  sqrt(2),1,0,-1;
  -1,sqrt(2),1,0;
  0,-1,sqrt(2),-1
];
# check unitarity of f
unity=f*f'

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
  ret = [0,0;0,0];
  f_adjugate = f_4_4';
  for k_ = 1:2
    for m_ = 1:2
      for i_ = 1:4
        for j_ = 1:4
          omega_ij = (E_4(i_)-E_4(j_))/hbar;
          exp_ij = exp(-i*omega_ij*t);
          for p_ = 1:2
            index_f = (k_-1)*2+(p_-1)+1;
            index_f_adj = (m_-1)*2+(p_-1)+1;
            ret(k_,m_) += d_4(i_)*d_4(j_)'*f_4_4(i_,index_f)
                          *f_adjugate(index_f_adj,j_)*exp_ij;
          endfor
        endfor
      endfor
    endfor
  endfor
endfunction

N_points=3;
t=linspace(0,0.8e-15,N_points); # femto seconds 1e-15
S1=linspace(0,0,N_points);
S2=linspace(0,0,N_points);
for t_ = 1:N_points
  rho = rho_4_4(energies, amplitudes, f, t(t_));
  S = -trace(rho*log2(rho))
  S1(t_) = real(S); # suppress imaginary part (numeric errors < 1e-15)
  rho = rho_4_4(energies, amplitudes, unity, t(t_))
  log2(rho)
  S = -trace(rho*log2(rho));
  S2(t_) = real(S);
endfor
plot(t,S1,"*;f;",t,S2,"*;unity matrix;");
xlabel("t");
ylabel("S");
