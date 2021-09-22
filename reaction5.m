% solving 
%          phi_t = -u*grad(phi)+Gamma*(...)
%          (phi*rhoa)_t = RHS
% use forward forward Euler scheme to solve phi equation,
%
%          phi^new - phi^old = dt*[ -u*grad(phi)+Gamma*(...) ]^old
%
% for rhoa equation, we use the following discretization
%
%  rhoa^new = rhoa^old - dt/phi^old*[ -u*grad(phi)+Gamma*(...) ]^old + dt/phi^old*RHS^old
% 
% which is the same as
% 
%   rhoa^new = (2*phi^old-phi^new)/phi^old*rhoa^old + dt/phi^old*RHS^old

% But the D_m term goes inside of the divergence

function [rhoa, rhom] = reaction5(rhoa0,rhom0,phi0,phi,u0,v0,dt,dx,dy,m,n,Lx,Ly,lambda,rhoa_tot,ka,kb,kc,Ka,KD,Da,Dm0)



     phi_box = (phi0>=lambda);
     
     phirhoa = phi0.*rhoa0;

     rhoa_area = sum(phirhoa(:))*dx*dy;
     phi0_area = sum(phi0(:))*dx*dy;
     rhoa_cyt  = (rhoa_tot - rhoa_area)/phi0_area;
     
     Dm       = Dm0./(1+rhoa0*KD);
     Da_term  = div_fd(Da *phi0,rhoa0,m,n,Lx,Ly);
     Dm_term  = div_fd(Dm.*phi0,rhom0,m,n,Lx,Ly); 
     f_rhoa   = kb*(rhoa0.^2./(Ka^2+rhoa0.^2)+ka)*rhoa_cyt - kc*rhoa0;
  
     uterm = -div4_fd(phi0,rhoa0,u0,v0,m,n,Lx,Ly);
     vterm = -div4_fd(phi0,rhom0,u0,v0,m,n,Lx,Ly);
 
     rhoa_rhs = Da_term + phi0.*f_rhoa + uterm;       
     prod1    = (2*phi0-phi).*rhoa0 + dt*rhoa_rhs;
     rhoa     = phi0.*rhoa0; rhoa(phi_box) = prod1(phi_box)./phi0(phi_box);
        
     rhom_rhs = Dm_term + vterm;
     prod2    = (2*phi0-phi).*rhom0 + dt*rhom_rhs;
     rhom     = phi0.*rhom0; rhom(phi_box) = prod2(phi_box)./phi0(phi_box);
     
     
     rhomtot0 = phi0.*rhom0; rhomtot0 = sum( rhomtot0(:) );
     rhomtot  = phi .*rhom ; rhomtot  = sum( rhomtot (:) );
     rhom = rhom*rhomtot0/rhomtot;