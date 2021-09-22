function [phi,rhoa,rhom,u,v,randx0,randy0,randxn,randyn] = periodicshift(phi,rhoa,rhom,u,v,randx0,randy0,randxn,randyn,dx,dy,m,n,Lx,Ly,lambda)

     if (max(phi(n-n/16,:))>=0.5)  
%          periodic_shifter_y = periodic_shifter_y + (n*dy/4);          
         box   = (phi>=lambda);         
         u2    = zeros(size(phi));    u2(box) =    u(box);    u2 = circshift(   u2,-n/4);
         v2    = zeros(size(phi));    v2(box) =    v(box);    v2 = circshift(   v2,-n/4);
         rhoa2 = zeros(size(phi)); rhoa2(box) = rhoa(box); rhoa2 = circshift(rhoa2,-n/4);
         rhom2 = zeros(size(phi)); rhom2(box) = rhom(box); rhom2 = circshift(rhom2,-n/4);
                                                             phi = circshift(  phi,-n/4);
         u = u2; v = v2; rhoa = rhoa2; rhom = rhom2;         
         randyn = randyn - 0.5*Ly;
         randy0 = randy0 - 0.5*Ly;
         
     elseif (max(phi(n/16,:))>=0.5)  
%          periodic_shifter_y = periodic_shifter_y - (n*dy/4);          
         box   = (phi>=lambda);         
         u2    = zeros(size(phi));    u2(box) =    u(box);    u2 = circshift(   u2,n/4);
         v2    = zeros(size(phi));    v2(box) =    v(box);    v2 = circshift(   v2,n/4);
         rhoa2 = zeros(size(phi)); rhoa2(box) = rhoa(box); rhoa2 = circshift(rhoa2,n/4);
         rhom2 = zeros(size(phi)); rhom2(box) = rhom(box); rhom2 = circshift(rhom2,n/4);
                                                             phi = circshift(  phi,n/4);
         u = u2; v = v2; rhoa = rhoa2; rhom = rhom2;         
         randyn = randyn + 0.5*Ly;
         randy0 = randy0 + 0.5*Ly; 
         
     elseif (max(phi(:,m/16))>=0.5)   
%          periodic_shifter_x = periodic_shifter_x - (n*dx/4);
         box   = (phi>=lambda);         
         u2    = zeros(size(phi));    u2(box) =    u(box);    u2 = circshift(   u2,[0,m/4]);
         v2    = zeros(size(phi));    v2(box) =    v(box);    v2 = circshift(   v2,[0,m/4]);
         rhoa2 = zeros(size(phi)); rhoa2(box) = rhoa(box); rhoa2 = circshift(rhoa2,[0,m/4]);
         rhom2 = zeros(size(phi)); rhom2(box) = rhom(box); rhom2 = circshift(rhom2,[0,m/4]);
                                                             phi = circshift(  phi,[0,m/4]);
         u = u2; v = v2; rhoa = rhoa2; rhom = rhom2;         
         randxn = randxn + 0.5*Lx;
         randx0 = randx0 + 0.5*Lx; 

     elseif (max(phi(:,m-m/16))>=0.5)   
%          periodic_shifter_x = periodic_shifter_x + (n*dx/4);
         box   = (phi>=lambda);         
         u2    = zeros(size(phi));    u2(box) =    u(box);    u2 = circshift(   u2,[0,-m/4]);
         v2    = zeros(size(phi));    v2(box) =    v(box);    v2 = circshift(   v2,[0,-m/4]);
         rhoa2 = zeros(size(phi)); rhoa2(box) = rhoa(box); rhoa2 = circshift(rhoa2,[0,-m/4]);
         rhom2 = zeros(size(phi)); rhom2(box) = rhom(box); rhom2 = circshift(rhom2,[0,-m/4]);
                                                             phi = circshift(  phi,[0,-m/4]);
         u = u2; v = v2; rhoa = rhoa2; rhom = rhom2;         
         randxn = randxn - 0.5*Lx;
         randx0 = randx0 - 0.5*Lx;         
     end   
     
end