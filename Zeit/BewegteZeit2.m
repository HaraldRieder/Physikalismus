pkg load symbolic

# Lorentz boost in x direction with beta=v/c
# Returns the transformed (t,x) as row vector
function retval = boost(t,x,beta)
  gamma=1./sqrt((1-beta.^2));
  retval = [gamma .* (t + beta * x),gamma .* (x + beta * t)];
endfunction

graphics_toolkit "gnuplot"

# corners of rectangular function psi
c0=[2,0.5];
c1=[2,-0.5];
c2=-c0;
c3=-c1;
# now rotate the points pi/4
rot=1/sqrt(2)*[1,-1;1,1];
c0=rot*c0'
c1=rot*c1'
c2=rot*c2'
c3=rot*c3'

# beta=v/c 
beta_max=0.99
beta = linspace(0, beta_max, 10)'
tx0 = boost(c0(2),c0(1),beta);
tx1 = boost(c1(2),c1(1),beta);
tx2 = boost(c2(2),c2(1),beta);
tx3 = boost(c3(2),c3(1),beta);

figure(1);
tx00(:,1)'
tx00(:,2)'
plot(tx0(:,1),tx0(:,2));
hold on
plot(tx1(:,1),tx1(:,2));
plot(tx2(:,1),tx2(:,2));
plot(tx3(:,1),tx3(:,2));
hold off
size=12 
axis ([-size, size,-size, size], "square");
xlabel ("x");
ylabel ("t");
title ({"4 corners",sprintf("boosted with \\beta=0 .. %f", beta_max)});
grid off
