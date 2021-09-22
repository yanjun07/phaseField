% computing div (f*g*u, f*g*v)

function ff = div4_fd( f, g, u, v, m, n, Lx, Ly )

dx = 2*Lx/m;
dy = 2*Ly/n;

fr = circshift(f,[ 0,-1]);
fl = circshift(f,[ 0, 1]);
fu = circshift(f,[ 1, 0]);
fd = circshift(f,[-1, 0]);

gr = circshift(g,[ 0,-1]);
gl = circshift(g,[ 0, 1]);
gu = circshift(g,[ 1, 0]);
gd = circshift(g,[-1, 0]);

ur = circshift(u,[ 0,-1]);
ul = circshift(u,[ 0, 1]);
vu = circshift(v,[ 1, 0]);
vd = circshift(v,[-1, 0]);

fgu_x = (fr+f)/2.*(gr+g)/2.*(ur+u)/2 - (fl+f)/2.*(gl+g)/2.*(ul+u)/2;
fgu_x = fgu_x/dx;

fgv_y = (fd+f)/2.*(gd+g)/2.*(vd+v)/2 - (fu+f)/2.*(gu+g)/2.*(vu+v)/2;
fgv_y = fgv_y/dy;

ff = fgu_x + fgv_y;
