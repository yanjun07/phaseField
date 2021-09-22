% computing d_x(g*d_x)[f], d_x(g*d_y)[f], d_y(g*d_x)[f], d_y(g*d_y)[f]

function [phiVx_x,phiVy_x,phiVx_y,phiVy_y] = div5_fd( phi, V, m, n, Lx, Ly )

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

Vrd = circshift(V,[ -1,-1]);
Vru = circshift(V,[  1,-1]);
Vld = circshift(V,[ -1, 1]);
Vlu = circshift(V,[  1, 1]);

Vxr = (Vr-V)/dx;
Vxl = (V-Vl)/dx;
Vxd = ( (V+Vr+Vd+Vrd)/4 - (V+Vl+Vd+Vld)/4 )/dx;
Vxu = ( (V+Vr+Vu+Vru)/4 - (V+Vl+Vu+Vlu)/4 )/dx;

Vyr = ( (V+Vr+Vd+Vrd)/4 - (V+Vr+Vu+Vru)/4 )/dy;
Vyl = ( (V+Vl+Vd+Vld)/4 - (V+Vl+Vu+Vlu)/4 )/dy;
Vyu = (V-Vu)/dy;
Vyd = (Vd-V)/dy;

phiVx_x = ( (phi+phir)/2.*Vxr - (phi+phil)/2.*Vxl )/dx;

phiVx_y = ( (phi+phid)/2.*Vxd - (phi+phiu)/2.*Vxu )/dy;

phiVy_x = ( (phi+phir)/2.*Vyr - (phi+phil)/2.*Vyl )/dx;

phiVy_y = ( (phi+phid)/2.*Vyd - (phi+phiu)/2.*Vyu )/dy;