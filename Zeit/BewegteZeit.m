pkg load symbolic

# Lorentz boost in x direction with beta=v/c
function retval = boost(t,x,beta)
  gamma=1/sqrt((1-beta^2));
  retval = gamma * (t + beta * x);
endfunction

# clock state
# approximation of the delta distribution in a limited area
function retval = psi_clock (t,x, beta)
  length = 500; # extension on the light cone
  width = length/10; # extension orthogonal to the light cone
  t_ = boost(t,x,beta);
  x_ = boost(x,t,beta);
  plus = (1/sqrt(2))*(t_+x_);
  minus = (1/sqrt(2))*(t_-x_);
  retval = heaviside(length-plus,0) .* heaviside(plus+length,0);
  retval = retval .* heaviside(width-minus,0) .* heaviside(minus+width,0);
endfunction


#t = x = linspace(-1000, 1000, 200)';
#[tt, xx] = meshgrid(t,x);
#psi = psi_clock(tt,xx,0.99); # 0, 0.7, 0.9, 0.99
#hidden off
##mesh(t,x,psi);
#surf(t,x,psi);
#xlabel ("t");
#ylabel ("x");
#title ("\\psi_{clock}");
#view(-10, 50);

# TODO Normierung von psi_clock
# TODO plus oder minus, das ist hier die Frage! Im boost und im psi_clock!

# 4 lines in the x,t plane
# enclose the area where the clock state is not zero

graphics_toolkit "gnuplot"

beta=0.999;
length=2000;
width=length/1000;
slope=1;

function retval = boost_slope(slope, beta)
  retval = (slope+beta)/(1+slope*beta);
endfunction
function retval = boost_intercept(intercept, beta)
  retval = intercept*sqrt(1-beta^2);
endfunction

slope_1 = 1#boost_slope(1,beta);
slope_2 = -1#boost_slope(-1,beta);
intercept_1 = boost_intercept(width/2*sqrt(2),beta);
intercept_2 = boost_intercept(length/2*sqrt(2),beta);

set(gca, 'xaxislocation', 'origin');
set(gca, 'yaxislocation', 'origin');
x=-1000:1:1000;
grid off
plot (x, slope_1 * x + intercept_1, "b");
hold on
plot (x, slope_1 * x - intercept_1, "b");
plot (x, slope_2 * x + intercept_2, "b");
plot (x, slope_2 * x - intercept_2, "b");
hold off
xlabel ("x");
ylabel ("t");
axis ([-1000,1000,-1000,1000]);
#axis("equal", "nolabel");
axis("equal");
title(sprintf("|\\beta| = %f", beta));

