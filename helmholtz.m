% for eqn a_bar * u - b_bar * Delta(u) = f
% given f, how to find u
function u = helmholtz(a_bar, b_bar, f, m, n, Lx, Ly)

k1       = [0:m/2, -(m/2-1):-1];
k2       = [0:n/2, -(n/2-1):-1];
[k1,k2] = meshgrid(k1, k2);
k1_sqr  = k1.^2;
k2_sqr  = k2.^2;


f1  = fft2(f);
u1  = f1./( a_bar + b_bar*(pi/Lx)^2*k1_sqr + b_bar*(pi/Ly)^2*k2_sqr );
u = ifft2(u1);