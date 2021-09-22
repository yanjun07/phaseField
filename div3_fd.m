% computing div (f*u, f*v)

function ff = div3_fd( f, u, v, m, n, Lx, Ly )

dx = 2*Lx/m;
dy = 2*Ly/n;

fr = circshift(f,[ 0,-1]);
fl = circshift(f,[ 0, 1]);
fu = circshift(f,[ 1, 0]);
fd = circshift(f,[-1, 0]);
ur = circshift(u,[ 0,-1]);
ul = circshift(u,[ 0, 1]);
vu = circshift(v,[ 1, 0]);
vd = circshift(v,[-1, 0]);

fu_x = (fr+f)/2.*(ur+u)/2 - (fl+f)/2.*(ul+u)/2;
fu_x = fu_x/dx;

fv_y = (fd+f)/2.*(vd+v)/2 - (fu+f)/2.*(vu+v)/2;
fv_y = fv_y/dy;

ff = fu_x + fv_y;
