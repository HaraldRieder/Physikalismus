# Time evolution of the von Neumann entropy of a 1 qubit density matrix
# according to Schrödinger dynamics.

# example one: energy ratios are rational numbers -> periodic entropy
energies=[1,2,5,7]; # eV

# example two: irrational energy ratios -> quasi periodic entropy
energies2=[1,2,5,2.6*e]; # eV, 1 non-rational number

# example amplitudes of our state vector in energy representation
amplitudes=[1,2,5,7];
# normalize |psi>
amplitudes*=1/sqrt(amplitudes*amplitudes');

# example unitary matrix f, defines the split of the 2 qubit space
f=0.5*[
  1,0,1,sqrt(2);
  sqrt(2),1,0,-1;
  -1,sqrt(2),1,0;
  0,-1,sqrt(2),-1
];
# check unitarity of f
identity=f*f'

# Planck's constant over 2 pi in eV*s units.
function retval = hbar()
  retval = 6.582119569e-16; 
endfunction

# period length of periodic example one
printf("period is %10.3e seconds\n",2*pi*hbar); # 1/1, 1/2, 1/3, 1/4, 1/5, 1/6 -> 1/1 common period

# von Neumann entropy of a density matrix with base 2 logarithm
function retval = entropy(matrix)
  retval = -trace(matrix*logm(matrix))*log2(e);
endfunction

#pure_state_density_matrix=amplitudes'*amplitudes
#S=entropy(pure_state_density_matrix)

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
            ret(k_,m_) += d_4(i_)*d_4(j_)'*f_4_4(i_,index_f)*f_adjugate(index_f_adj,j_)*exp_ij;
          endfor
        endfor
      endfor
    endfor
  endfor
endfunction

N_points=500;
t=linspace(0,10.0e-15,N_points); # femto seconds 1e-15
S1=S2=S3=S4=linspace(0,0,N_points);
for t_ = 1:N_points
  rho = rho_4_4(energies, amplitudes, identity, t(t_));
  S = entropy(rho);
  S1(t_) = real(S);
  rho = rho_4_4(energies, amplitudes, f, t(t_));
  #trc=trace(rho) # must be 1
  #trc_square=trace(rho*rho); # must be less or equal 1
  S = entropy(rho);
  S2(t_) = real(S); # suppress imaginary part (numeric errors < 1e-15)
  rho = rho_4_4(energies2, amplitudes, identity, t(t_));
  S = entropy(rho);
  S3(t_) = real(S);
  rho = rho_4_4(energies2, amplitudes, f, t(t_));
  S = entropy(rho);
  S4(t_) = real(S);
endfor
hf=figure(1, 'position', [110 400 1200 400]);
h=plot(t,S1,"-;identity matrix;",t,S2,"-;f;");
xlabel("t/s");
ylabel("S");
print (hf, "periodic_entropy.png", "-dpng","-S1200,400");
hf=figure(2, 'position', [110 40 1200 400]);
h=plot(t,S3,"-;identity matrix;",t,S4,"-;f;");
xlabel("t/s");
ylabel("S");
print (hf, "quasi_periodic_entropy.png", "-dpng","-S1200,400");
