% computing div (phi * grad V)

function f = div_fd( phi, V, m, n, Lx, Ly )

dx = 2*Lx/m;
dy = 2*Ly/n;

phir = circshift(phi,[ 0,-1]);
phil = circshift(phi,[ 0, 1]);
phid = circshift(phi,[-1, 0]);
phiu = circshift(phi,[ 1, 0]);

Vr = circshift(V,[ 0,-1]);
Vl = circshift(V,[ 0, 1]);
Vd = circshift(V,[-1, 0]);
Vu = circshift(V,[ 1, 0]);

Vxr = (Vr-V)/dx;
Vxl = (V-Vl)/dx;
Vyu = (V-Vu)/dy;
Vyd = (Vd-V)/dy;

phiVx_x = (phi+phir)/2.*Vxr - (phi+phil)/2.*Vxl;
phiVx_x = phiVx_x/dx;

phiVy_y = (phi+phid)/2.*Vyd - (phi+phiu)/2.*Vyu;
phiVy_y = phiVy_y/dy;

f = phiVx_x + phiVy_y;

