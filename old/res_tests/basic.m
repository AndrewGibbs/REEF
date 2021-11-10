poles = [0.01 0.001i -0.005];
F = @(x) x;
f = @(x) F(x)./polyval(poly(poles),x);
f_polar = @(theta) f(exp(1i*theta));
I = integral(f_polar,0,2*pi,'AbsTol',1e-14)/(2*pi*1i);

R = 0;
P = 1;

not_poles = [2 3; 3 1; 1 2;];
for n=1:3
   R = R + F(poles(n))*(poles(not_poles(n,1))-poles(not_poles(n,2)));
   P = P*(poles(not_poles(n,1))-poles(not_poles(n,2)));
end

Ires = R/P;

err = abs(I-Ires)