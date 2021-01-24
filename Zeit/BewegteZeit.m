pkg load symbolic

# Lorentz boost in x direction with beta=v/c
# Returns the transformed (t,x) as row vector
function retval = boost(t,x,beta)
  gamma=1/sqrt((1-beta^2));
  retval = [gamma * (t + beta * x),gamma * (x + beta * t)];
endfunction

# clock state
# approximation of the delta distribution by rectangular funcctions
function retval = psi_rect (t,x, beta)
  # TODO length und width müssen auch transformiert werden
  length = 500; # extension on the light cone
  width = length/10; # extension orthogonal to the light cone
  t_ = boost(t,x,beta);
  x_ = boost(x,t,beta);
  plus = (1/sqrt(2))*(t_+x_);
  minus = (1/sqrt(2))*(t_-x_);
  retval = heaviside(length-plus,0) .* heaviside(plus+length,0);
  retval = retval .* heaviside(width-minus,0) .* heaviside(minus+width,0);
endfunction


function retval = boost_a(a, beta)
  gamma=1/sqrt((1-beta^2));
  retval=a/(gamma^2*(1+beta)^2);  
endfunction

function retval = boost_b(b, beta)
  gamma=1/sqrt((1-beta^2));
  retval=b/(gamma^2*(1-beta)^2);  
endfunction

function retval = psi_gauss(t,x,a,b)
  # normalization 
  N = 2*sqrt(sqrt(a*b)/pi);
  retval = N * exp(-a*(x-t).^2 - b*(x+t).^2);
endfunction

a=0.1
b=0.001
beta1=0;
beta2=0.9;
beta3=0.999;

graphics_toolkit "gnuplot"

t = x = linspace(-20, 20, 50)';

figure(1);
boosted = boost(t,x,beta1);
boosted_t = boosted(:,1);
boosted_x = boosted(:,2);
[tt,xx] = meshgrid(boosted_t,boosted_x);
psi = psi_gauss(tt,xx,a,b);
mesh(t,x,psi);
view (-35, 25);
xlabel ("t");
ylabel ("x");
zlabel ("");
title ({"\\psi";sprintf("a/b=%d, \\beta=%d",a/b,beta1)});
colormap("winter");
figure(2);
boosted = boost(t,x,beta2);
boosted_t = boosted(:,1);
boosted_x = boosted(:,2);
[tt,xx] = meshgrid(boosted_t,boosted_x);
psi = psi_gauss(tt,xx,a,b);
#psi = psi_gauss(tt,xx,boost_a(a,beta2),boost_b(b,beta2));
mesh(t,x,psi);
view (-35, 25);
xlabel ("t");
ylabel ("x");
zlabel ("");
title ({"\\psi";sprintf("a/b=%d, \\beta=%d",a/b,beta2)});
colormap("winter");
figure(3);
boosted = boost(t,x,beta3);
boosted_t = boosted(:,1);
boosted_x = boosted(:,2);
[tt,xx] = meshgrid(boosted_t,boosted_x);
psi = psi_gauss(tt,xx,a,b);
#psi = psi_gauss(tt,xx,boost_a(a,beta3),boost_b(b,beta3));
mesh(t,x,psi);
view (-35, 25);
xlabel ("t");
ylabel ("x");
zlabel ("");
title ({"\\psi";sprintf("a/b=%d, \\beta=%d",a/b,beta3)});
colormap("winter");
#axis("off","tight","square");
grid off


# 4 lines in the x,t plane
# enclose the area where the clock state is not zero

length=2000;
width=length/8;
slope=1;

function retval = boost_slope(slope, beta)
  retval = (slope+beta)/(1+slope*beta);
endfunction
function retval = boost_intercept(intercept, beta)
  retval = intercept*sqrt(1-beta^2);
endfunction


figure(4);
set(gca, 'xaxislocation', 'origin');
set(gca, 'yaxislocation', 'origin');
x=-1000:1:1000;
slope_1 = boost_slope(1,beta1);
slope_2 = boost_slope(-1,beta1);
intercept_1 = boost_intercept(width/2*sqrt(2),beta1);
intercept_2 = boost_intercept(length/2*sqrt(2),beta1);
grid off
h = plot (
  x, slope_1 * x + intercept_1, sprintf("b;\\beta=%d ;",beta1), 
  x, slope_1 * x - intercept_1, "b",
  x, slope_2 * x + intercept_2, "b",
  x, slope_2 * x - intercept_2, "b");
set (h, "linewidth", 1);  
hold on
slope_1 = boost_slope(1,beta2);
slope_2 = boost_slope(-1,beta2);
intercept_1 = boost_intercept(width/2*sqrt(2),beta2);
intercept_2 = boost_intercept(length/2*sqrt(2),beta2);
h = plot (
  x, slope_1 * x + intercept_1, sprintf("g;\\beta=%d ;",beta2), 
  x, slope_1 * x - intercept_1, "g",
  x, slope_2 * x + intercept_2, "g",
  x, slope_2 * x - intercept_2, "g");
set (h, "linewidth", 1);  
slope_1 = boost_slope(1,beta3);
slope_2 = boost_slope(-1,beta3);
intercept_1 = boost_intercept(width/2*sqrt(2),beta3);
intercept_2 = boost_intercept(length/2*sqrt(2),beta3);
h = plot (
  x, slope_1 * x + intercept_1, sprintf("r;\\beta=%d ;",beta3),
  x, slope_1 * x - intercept_1, "r",
  x, slope_2 * x + intercept_2, "r",
  x, slope_2 * x - intercept_2, "r");
set (h, "linewidth", 1);  
hold off
xlabel ("x");
ylabel ("t");
axis ([-1000,1000,-1000,1000]);
#axis("equal", "nolabel");
axis("equal");
title("boosts of a rectangle");

