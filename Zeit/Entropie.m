graphics_toolkit "gnuplot"

a=10
b=0.001


function retval = psi(t,x,a,b)
  # normalization 
  N = 2*sqrt(sqrt(a*b)/pi);
  retval = N * exp(-a*(x-t).^2 - b*(x+t).^2);
endfunction

figure (1);
t = x = linspace (-20, 20, 80)';
[tt,xx] = meshgrid(t,x);
psi = psi(tt,xx,a,b);
mesh(t,x,psi);
view (-35, 30);
xlabel ("t");
ylabel ("x");
zlabel ("");
title ({"\\psi";sprintf("a/b=%d",a/b)});
axis("off","tight","square");
grid off

function retval = entropy(a,b)
  retval = (1 .- log(8*a.*b./(pi*(a+b))))/2;
endfunction

# TODO leider haben wir hier negative Entropiewerte dabei
figure (2);
a = b = linspace (0.001, 100, 50)';
[aa,bb] = meshgrid(a,b);
S=entropy(aa,bb);
mesh(a,b,S);
xlabel ("a");
ylabel ("b");
title ("S(a,b)");
view (50, 30);
