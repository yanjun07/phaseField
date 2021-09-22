% given vector x, find the mesh points nearest to x from
% the meshgrid -Lx:dx:Lx-dx

function x_grid = getpos( x, dx, x_left )

x_grid = round((x-x_left)/dx) + 1;
