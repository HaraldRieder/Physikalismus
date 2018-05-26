t_shortest = 0.7e-18; # seconds, shortest time ever measured 2017
eV = 1.602e-19; # Joule
GeV = 1.602e-10; # Joule
TeV = 1.602e-7; # Joule

function t_planck = t_planck() 
  t_planck = 5.391e-44; # seconds
endfunction

function h_planck = h_planck() 
  h_planck = 1.054571e-34; # h/2pi Joule seconds
endfunction

# Annäherung e-Funktion über Polynom
function e_approximation = e_approximation(x,n) 
  e_approximation = (1 + x./n).^n;
endfunction

# Relativer Fehler der Polynom-Näherung 
function relative_error = relative_error(x,n)
  a = e_approximation(x,n);
  relative_error = abs(1 - a./(e.^x));  
endfunction

# energy in Joule to exponent x
function energy_to_x = energy_to_x(energy) 
  energy_to_x = energy .* t_planck ./ h_planck;
endfunction

graphics_toolkit('gnuplot');

n = t_shortest / t_planck
n=10000
subtitle=sprintf("n=%s",num2str(n));
x = -500:0.02:500;

figure(1)
err = relative_error(x,n);
semilogy(x,err);
title({"relative error";subtitle});

figure(2)
semilogy(x,e.^x,x,e_approximation(x,n));
title({"e^x versus (1+x/n)^n";subtitle});
legend("e^x","(1+x/n)^n");

_1eV_as_x = energy_to_x(eV)
printf("1 eV relative error\n");
err = relative_error(_1eV_as_x,n)

_1TeV_as_x = energy_to_x(TeV)
printf("1 TeV relative error\n");
err = relative_error(_1TeV_as_x,n)

figure(3)
max=_1TeV_as_x*1000
x = 0:max/500:max;
err = relative_error(x,n);
plot(x,err);
title({"relative error up to 1000 TeV";subtitle});
xlabel("x");

