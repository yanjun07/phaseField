% given phi0, find curvature
function [curv_term,Fx,Fy] = force(phi,dphix,dphiy,grad_phi_sqr,rhoa,rhom,m,n,dx,dy,Lx,Ly,etaa0,etam0,epsilon,gamma0,kappa0,Peri_critical,phi_tilde,g,sigma)

dG   = 36*phi.*(1-phi).*(1-2*phi);
ddG  = 36*(1-6*phi + 6*phi.^2);
grad_phi_abs = sqrt(grad_phi_sqr);
box = (grad_phi_abs > 0.01);
lap_phi = lap(phi,m,n,Lx,Ly);

Peri = sum(grad_phi_abs(:))*dx*dy; % cell perimeter to determine the cell tension
  
[phirhomx,phirhomy] = grad(phi.*rhom,m,n,Lx,Ly);
Fmyox = etam0*phirhomx;
Fmyoy = etam0*phirhomy;

  
% if nargin == 21
%    phirhoa = phi.*rhoa.*band;
% else
%    phirhoa = phi.*rhoa;
% end

phirhoa = phi.*rhoa;

phirhoadphixx = phirhoa.*dphix.*dphix;
phirhoadphixy = phirhoa.*dphix.*dphiy;
phirhoadphiyy = phirhoa.*dphiy.*dphiy;
  
Fpolymx = -etaa0*epsilon*div(phirhoadphixx,phirhoadphixy,m,n,Lx,Ly);
Fpolymy = -etaa0*epsilon*div(phirhoadphixy,phirhoadphiyy,m,n,Lx,Ly);

% grad_phi_norm1 = zeros(n,m);
% grad_phi_norm2 = zeros(n,m);
% grad_phi_norm1(box) = dphix(box)./grad_phi_abs(box);
% grad_phi_norm2(box) = dphiy(box)./grad_phi_abs(box);
% grad_phi_norm1r = circshift(grad_phi_norm1,[0,-1]);
% grad_phi_norm1l = circshift(grad_phi_norm1,[0, 1]);
% grad_phi_norm2d = circshift(grad_phi_norm2,[-1,0]);
% grad_phi_norm2u = circshift(grad_phi_norm2,[ 1,0]);
% curv = - ( grad_phi_norm1r - grad_phi_norm1l )/(2*dx) -...
%            ( grad_phi_norm2d - grad_phi_norm2u )/(2*dy); 


phi_i_j = phi;
phi_ip1_j   = circshift(phi,[ 0,-1]);
phi_im1_j   = circshift(phi,[ 0, 1]);
phi_i_jp1   = circshift(phi,[-1, 0]);
phi_i_jm1   = circshift(phi,[ 1, 0]);
phi_ip1_jp1 = circshift(phi,[-1,-1]);
phi_ip1_jm1 = circshift(phi,[ 1,-1]);
phi_im1_jp1 = circshift(phi,[-1, 1]);
phi_im1_jm1 = circshift(phi,[ 1, 1]);

phix_iphalf_j = (phi_ip1_j - phi_i_j  )/dx;
phix_imhalf_j = (phi_i_j   - phi_im1_j)/dx;
phiy_i_jphalf = (phi_i_jp1 - phi_i_j  )/dy;
phiy_i_jmhalf = (phi_i_j   - phi_i_jm1)/dy;

phiy_iphalf_j = (phi_ip1_jp1 + phi_i_jp1   - phi_ip1_jm1 - phi_i_jm1  )/(4*dy);
phiy_imhalf_j = (phi_i_jp1   + phi_im1_jp1 - phi_i_jm1   - phi_im1_jm1)/(4*dy);
phix_i_jphalf = (phi_ip1_jp1 + phi_ip1_j   - phi_im1_jp1 - phi_im1_j  )/(4*dx);
phix_i_jmhalf = (phi_ip1_j   + phi_ip1_jm1 - phi_im1_j   - phi_im1_jm1)/(4*dx);


grad_phi_abs_iphalf_j = sqrt( phix_iphalf_j.^2 + phiy_iphalf_j.^2 );
grad_phi_abs_imhalf_j = sqrt( phix_imhalf_j.^2 + phiy_imhalf_j.^2 );
grad_phi_abs_i_jphalf = sqrt( phix_i_jphalf.^2 + phiy_i_jphalf.^2 );
grad_phi_abs_i_jmhalf = sqrt( phix_i_jmhalf.^2 + phiy_i_jmhalf.^2 );

curv = zeros(n,m);
curv(box) = - ( phix_iphalf_j(box)./grad_phi_abs_iphalf_j(box) - phix_imhalf_j(box)./grad_phi_abs_imhalf_j(box) )/dx -...
              ( phiy_i_jphalf(box)./grad_phi_abs_i_jphalf(box) - phiy_i_jmhalf(box)./grad_phi_abs_i_jmhalf(box) )/dy;
curv_term = curv.*grad_phi_abs + lap_phi;     


ten  = epsilon*lap_phi - dG/epsilon;
bend = lap(ten,m,n,Lx,Ly) - 1/epsilon^2*ddG.*ten;

if nargin <= 19 % single-cell case

   gamma = gamma0*(1 + 0.0*( max(Peri_critical,Peri) - Peri_critical ) ) ;
   Fmem = -gamma*ten + kappa0*bend;     
   Fmemx = Fmem.*dphix;
   Fmemy = Fmem.*dphiy;  
   
else % two-cell case
    
   lap_phi_tilde = lap(phi_tilde,m,n,Lx,Ly);
   Fcellcell = g/2*phi_tilde - sigma*epsilon*lap_phi_tilde;
   
   gamma = gamma0*(1 + 0.0*( max(Peri_critical,Peri) - Peri_critical ) ) ;
   Fmem = -gamma*ten + kappa0*bend + Fcellcell; 
   Fmemx = Fmem.*dphix;
   Fmemy = Fmem.*dphiy;
   
end

  
Fx = Fmyox + Fpolymx + Fmemx;
Fy = Fmyoy + Fpolymy + Fmemy;   