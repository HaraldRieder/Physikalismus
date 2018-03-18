m_proton = 1.672621898e-27 # kg
m_neutrino_max = 4e-36 # kg
m_sun = 1.9884e30 # kg
r_sun = 6.96342e8 # m
A_sun = 4*pi*r_sun^2; # m^2
V_sun = 4/3*pi*r_sun^3; # m^3
T_sun_max = 1.6e7 # K
T_sun_min = 5.778e3 # surface K

# Boltzmann J/K
function k = k()
  k = 1.3806504e-23; 
endfunction 

# Planck Js
function h = h()
  h = 6.6260755e-34; 
endfunction

# gravitational constant
function G = G()
  G = 6.67408e-11;
endfunction

# speed of light
function c = c()
  c = 299792458;
endfunction  

# volume of a sphere
function V = Volume(r)
  V = 4/3*pi*r^3;
endfunction

# MKSA statistical mechanics ideal gas with Boltzmann constant
# T temperature
# V volume
# M total mass
# m particle mass
function S = S_ideal_gas_classic(T,V,M,m)
  N = M./m;
  S = k*3/2*N*log(2*pi*m*k*T*V.^(2/3)/(h.^2 * N.^(2/3))) + k*5/2*N;
endfunction

# MKSA statistical mechanics ideal gas with Boltzmann constant
# T temperature
# A surface area
# M total mass
# m particle mass
function S = S_ideal_gas_classic_2(T,A,M,m)
  N = M/m;
  S = k*3/2*N*log(2*pi*m*k*T* A/(4*pi) * (4*pi/3)^(2/3) /(h.^2 * N.^(2/3))) + k*5/2*N;
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

# MKSA area of Schwarzschild black hole with mass M
function A = A_black_hole_Schwarzschild(M)
  A = 16*pi*M^2*G^2/(c^4);
endfunction

graphics_toolkit('gnuplot');

T_bh = T_black_hole_Schwarzschild(m_sun)
S_bh = S_black_hole_Schwarzschild(m_sun) / k
A_bh = A_black_hole_Schwarzschild(m_sun)

figure(1)
T = logspace(1,log10(T_sun_max),500);
S_classic = S_ideal_gas_classic(T,V_sun,m_sun,m_proton) / k;
#S_classic_neutrino = S_ideal_gas_classic(T,V_sun,m_sun,m_neutrino_max) / k;
#loglog(T,S_classic,T,S_classic_neutrino);
loglog(T,S_classic,T_bh,S_bh);
xlabel("Temperatur / K");
ylabel("Entropie / k \n");
legend("ideales Gas klassisch","schwarzes Loch");
#title("Entropie einer Sonnenmasse");

figure(2)
T = T_sun_max;
A = logspace(13,log10(4*pi*r_sun^2),500);
S_classic = S_ideal_gas_classic_2(T,A,m_sun,m_proton) / k;
loglog(A,S_classic,A_bh,S_bh);
xlabel("Oberfläche / m^2");
ylabel("Entropie / k \n");
legend("ideales Gas klassisch","schwarzes Loch");
#title("Entropie einer Sonnenmasse");


figure(3)
T = T_bh:T_sun_max/1000:T_sun_max;
S_classic = S_ideal_gas_classic(T,V_sun,m_sun,m_proton) / k;
semilogy(T,S_classic,T_bh,S_bh);
xlabel("Temperatur / K");
ylabel("Entropie / k \n");
legend("ideales Gas klassisch","schwarzes Loch");
#title("Entropie einer Sonnenmasse");

figure(4)
T = T_sun_max;
A = A_bh:A_sun/1000:A_sun;
S_classic = S_ideal_gas_classic_2(T,A,m_sun,m_proton) / k;
semilogy(A,S_classic,A_bh,S_bh);
xlabel("Oberfläche / m^2");
ylabel("Entropie / k \n");
legend("ideales Gas klassisch","schwarzes Loch");
#title("Entropie einer Sonnenmasse");







