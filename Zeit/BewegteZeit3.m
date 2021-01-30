n=5
# identity matrix t,x
psi=eye(n);
#rearranged as vector of n^2 rows
psi=reshape(psi.', [], 1);
# normalize vector
psi=psi/(psi'*psi);
# density matrix
rho=psi*psi';

function ret = rho_reduced(psi)
  r=rows(psi);
  max=sqrt(r);
  ret = zeros(max);
  for k_ = 1:max
    for m_ = 1:max
      for i_ = 1:max
        ret(k_,m_)+=psi(i_);
      endfor
    endfor
  endfor
endfunction


rho_red = rho_reduced(psi);

# von Neumann entropy of a density matrix with base 2 logarithm
function retval = entropy(matrix)
  retval = -trace(matrix*logm(matrix)*log2(e));
endfunction

S=entropy(rho_red)

graphics_toolkit "gnuplot"

