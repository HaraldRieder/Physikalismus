m_proton = 1.672621898e-27 # kg
m_sun = 1.9884e30 # kg
r_sun = 6.96342e8 # m
V_sun = 4/3*pi()*r_sun^3; # m^3
T_sun_max = 1.6e7 # K
T_sun_min = 5.778e3 # surface K

function k = k()
  k = 1.3806504e-23; # Boltzmann J/K
endfunction 

function h = h()
  h = 6.6260755e-34; # Planck Js
endfunction

# MKSA statistical mechanics ideal gas with Boltzmann constant
# M total mass
# m particle mass
# V volume
function S = S_ideal_gas_classic(T,V,M,m)
  S = log(T);
  N = M/m;
  S = k*3/2*N*log(2*pi()*m*k*T*V.^(2/3)/(h.^2 * N.^(2/3))) + k*5/2*N;
endfunction

# MKSA Bekenstein-Hawking entropy with Boltzmann constant
function S = S_black_hole_Schwarzschild(M)
  C = 2.65e+40;
  S = k * C *(M/1000000000000)^2;
endfunction

# MKSA Hawking temperature of a Schwarzschild black hole
function T = T_black_hole_Schwarzschild(M)
  T = 1.23e+23 / M;
endfunction

# MKSA mass of Schwarzschild black hole with temperature T
function M = M_black_hole_Schwarzschild(T)
  M = 1.23e+23 / T; 
endfunction

# MKSA radius of Schwarzschild black hole with mass M
function R = R_black_hole_Schwarzschild(M)
  R = 1.49e-27 * M;
endfunction

graphics_toolkit('gnuplot');
#T = T_sun_min:100000:T_sun_max
#printf("N=%d", m_sun/m_proton)
#plot(T, S_ideal_gas_classic(T,V_sun,m_sun,m_proton)/k)
#xlabel("T / K");
#ylabel("S");
#title("S ideales Gas - klassisch");

T = 1 # Kelvin
M_bh = M_black_hole_Schwarzschild(T)
S_bh = S_black_hole_Schwarzschild(M_bh)/k
R_bh = R_black_hole_Schwarzschild(M_bh)
r = R_bh:R_bh:R_bh*100; 
V = 4/3*pi*r.^3;
#plot(R_bh,S_bh, r,S_ideal_gas_classic(T,V,m_sun,m_proton)/k)
#plot(R_bh,S_bh, r,S_ideal_gas_classic(T,V,m_sun,m_proton)/k);
S_classic = S_ideal_gas_classic(T,V,m_sun,m_proton)/k
semilogy(R_bh,S_bh, r,S_classic);

xlabel("Radius/m");
ylabel("Entropie");
title("Entropie schwarzes Loch vs. klassisch");
