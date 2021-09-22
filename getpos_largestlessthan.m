% given vector x, find the Largest mesh point less than x from 
% the meshgrid -Lx:dx:Lx-dx

function x_grid = getpos_largestlessthan( x, dx, x_left )

x_grid = floor((x-x_left)/dx) + 1;

