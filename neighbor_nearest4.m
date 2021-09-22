% given meshgrid [-Lx:dx:Lx-dx]*[-Ly:dy:Ly-dy], and u on meshgrid. For
% vectors x and y, which might not on meshgrid, evaluate u(x,y).
% Notice that (x,y) might not be on the meshgrid, one need find the nearest
% neighbor of (x,y), and its u-value.
% size(u) = n*m, x,y = 1*N

function u_av = neighbor_nearest4( u, x, y, n, dx, dy, Lx, Ly )

u = u(:);

x_grid = getpos_largestlessthan(x, dx, -Lx);
y_grid = getpos_largestlessthan(y, dy, -Ly);

u_av = 0.25* ( u( (x_grid-1)*n + y_grid ) + u( (x_grid-1)*n + y_grid + 1 ) + u( x_grid*n + y_grid ) + u( x_grid*n + y_grid + 1 ) );

u_av = u_av';
