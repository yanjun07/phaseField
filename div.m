% Divergence of (u1, u2) on domain [-Lx, Lx]*[-Ly, Ly] using
% Spectral method.
function uu = div(u1, u2, m, n, Lx, Ly)

k1 = [0:m/2-1, 0,  -(m/2-1):-1];
k2 = [0:n/2-1, 0,  -(n/2-1):-1];
[k1, k2] = meshgrid(k1, k2);

u11  = fft2(u1);
u22  = fft2(u2);
u1x  = ifft2(u11.*k1*1i)*pi/Lx;
u2y  = ifft2(u22.*k2*1i)*pi/Ly;

uu = u1x + u2y;
