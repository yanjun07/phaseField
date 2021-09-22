function [fx, fy] = grad_fd(f, m, n, Lx, Ly)

dx = 2*Lx/m;
dy = 2*Ly/n;

fx = ( circshift(f,[0,-1]) - circshift(f,[0,1]) )/(2*dx);
fy = ( circshift(f,[-1,0]) - circshift(f,[1,0]) )/(2*dy);


