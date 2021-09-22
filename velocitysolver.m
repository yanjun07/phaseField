function [u,v,Fgrip,randxn,randyn,Fgrip1_mx,Fgrip2_mx,Fslip1_mx,Fslip2_mx,tadh] = velocitysolver(F1,F2,tadh,dt,phi,m,n,dx,dy,u0,v0,randxo,randyo,Lx,Ly,kgrip0,kslip0,randx0,randy0,label_slip,label_grip,Nadh,nu0,trick,xi,phi_box,step,mod_N)
  
  tadh = tadh + dt;  

  max_err     = 10;
  it_steps    = 0;
  error_limit = 0.01;
  max_steps   = 20;
 
  while ( (max_err > error_limit) && ( it_steps < max_steps ) ) 
            
            % find the speed (u,v) at the Nadh adhesion sites
       randu0 = neighbor_nearest4(u0,randxo,randyo,n, dx, dy, Lx, Ly);
       randv0 = neighbor_nearest4(v0,randxo,randyo,n, dx, dy, Lx, Ly);
            % Adhension force calculation   
       randxn = randxo + dt*randu0;
       randyn = randyo + dt*randv0;
            
       kgrip  = kgrip0.*tadh; 
       kslip  = kslip0.*tadh;
       
       Fgrip1 = -kgrip.*(randxn-randx0); Fgrip1(label_slip) = 0;
       Fgrip2 = -kgrip.*(randyn-randy0); Fgrip2(label_slip) = 0;
       Fslip1 = -kslip.*randu0;          Fslip1(label_grip) = 0;
       Fslip2 = -kslip.*randv0;          Fslip2(label_grip) = 0;
        
       x_grid = getpos(randxn,dx,-Lx); 
       y_grid = getpos(randyn,dy,-Ly); 
        
       Fgrip1_mx = full( sparse(y_grid,x_grid,Fgrip1,n,m,Nadh) );     
       Fgrip2_mx = full( sparse(y_grid,x_grid,Fgrip2,n,m,Nadh) );     
       Fslip1_mx = full( sparse(y_grid,x_grid,Fslip1,n,m,Nadh) );    
       Fslip2_mx = full( sparse(y_grid,x_grid,Fslip2,n,m,Nadh) ); 
   
       [vx,vy]  = grad(v0,m,n,Lx,Ly);
       [ux,uy]  = grad(u0,m,n,Lx,Ly);
   
       phivx = phi.*vx;
       phivy = phi.*vy;
       phiux = phi.*ux;
       phiuy = phi.*uy;

       [~,phivx_y] = grad(phivx,m,n,Lx,Ly);
       [~,phivy_y] = grad(phivy,m,n,Lx,Ly);
       [phiux_x,~] = grad(phiux,m,n,Lx,Ly);
       [phiuy_x,~] = grad(phiuy,m,n,Lx,Ly);
 
       Fadhn1 = Fgrip1_mx + Fslip1_mx;
       Fadhn2 = Fgrip2_mx + Fslip2_mx;
   
       u_rhs = nu0*( phivx_y + phiux_x ) + F1 + Fadhn1;
       v_rhs = nu0*( phiuy_x + phivy_y ) + F2 + Fadhn2;  
    
       trickux = (phi-trick).*ux;
       trickuy = (phi-trick).*uy;
       trickvx = (phi-trick).*vx;
       trickvy = (phi-trick).*vy;
   
       u_rhs = u_rhs + nu0*div(trickux,trickuy,m,n,Lx,Ly);
       v_rhs = v_rhs + nu0*div(trickvx,trickvy,m,n,Lx,Ly);

       u = helmholtz(xi,nu0*trick,u_rhs,m,n,Lx,Ly);
       v = helmholtz(xi,nu0*trick,v_rhs,m,n,Lx,Ly);   
   
       u2 = zeros(size(u)); u2(phi_box) = u(phi_box); u = u2;
       v2 = zeros(size(v)); v2(phi_box) = v(phi_box); v = v2;
  
       u_max = max( abs(u(:)) );
       v_max = max( abs(v(:)) );
       u_rel = max( ( abs(u(:)-u0(:)) ) /u_max );
       v_rel = max( ( abs(v(:)-v0(:)) ) /v_max );
       max_err = max(u_rel, v_rel);
       it_steps = it_steps + 1;
       u0 = u; v0 = v;
  
  end
  
  Fgrip = sqrt(Fgrip1.^2 + Fgrip2.^2); 
  
  if (mod(step, mod_N)==0) 
  fprintf('Iteration %1.f\v Error %1.5f\n',it_steps, max_err);
  end