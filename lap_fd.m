% find Delta(f) = f_xx + f_yy

function f_lap = lap_fd(f,dx,dy)

fr = circshift(f,[0,-1]);
fl = circshift(f,[0, 1]);  
fd = circshift(f,[-1,0]);
fu = circshift(f,[ 1,0]); 
  
fxx = ( fr - 2*f + fl )/(dx^2);
fyy = ( fd - 2*f + fu )/(dy^2);
f_lap = fxx + fyy;


