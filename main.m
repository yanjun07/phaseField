clear all;
clc
close all
%% parameter values
%  global lambda rhoa_tot ka kb kc Ka KD Da Dm0 % for reaction5



nu0     = VIS_PARA; %1e3;   % viscosity
etam0   = CONTR_PARA; %100;   % contraction
etaa0   = 560;
Fgrip0  = 20/2.0; %10;  % force, smaller -> kgrip0 easier to break
kgrip0  = 1.0 * (log(nu0) - log(500) + 0.1) ; %25/10.0;   % hook coef. friction ~ x 
kslip0  = 25/100.0; %0.25; % friction ~ v
xi      = 0.5;
ron     = 5e-3;
roff0   = 2e-3;
rdie    = 0.2; % larger -> kslip0 easier to break
ka      = 0.01;
kb      = 10;
Ka      = 1;
kc      = 10;
Da      = 0.2;
Dm0     = 2;
KD      = 2;
gamma0   = 25;
kappa0   = 25;
rhoa_tot = 800;
rhom_0   = 0.3;
Nadh    = 4000;
N       = Nadh*100;
epsilon = 1;
Gamma   = 0.4;

T         = 200;
dt        = 4e-3;
lambda    = 1e-4;
mod_N     = 10;
trick     = 2;

where2save = sprintf('contra_%d_fric_%ddot%d.dat',etam0, floor(kgrip0), rem(100*kgrip0,100));

%% Initial set-up
m      = 2^9; 
n      = 2^9;
Lx     = 50;
Ly     = 50;
x      = linspace(-Lx, Lx, m+1); x = x(1:end-1); dx = 2*Lx/m;
y      = linspace(-Ly, Ly, n+1); y = y(1:end-1); dy = 2*Ly/n;
[x,y]  = meshgrid(x,y);
%% Initial phi, rhoa, rhom, u, v, 
r0 = 10; x0 =  -8; y0 =  -10; dis = sqrt((x-x0).^2+(y-y0).^2); 
phi0 = 0.5+0.5*tanh(3.0*(r0-dis)/epsilon);
rhoa0 = 1.8*ones(size(x)); rhoa0(y<y0) = 0; rhom0 = rhom_0*ones(size(x));
u0 = zeros(size(x)); v0 = zeros(size(x));
%% Initial cell perimeter to stablize the cell shape
dG            = 36*phi0.*(1-phi0).*(1-2*phi0);
[dphix,dphiy] = grad(phi0,m,n,Lx,Ly);
grad_phi_sqr  = dphix.^2+dphiy.^2;
grad_phi_abs  = sqrt(grad_phi_sqr);
Peri_critical = sum(grad_phi_abs(:))*dx*dy;
%% Initial  adhesion site
phirhoa = phi0.*rhoa0;
phirhoa_max = max(phirhoa(:));
f = phirhoa/phirhoa_max;
% generating N random points (rand_x, rand_y) in [0,Lx-dx]*[-Ly,Ly-dy]
randx = -Lx + (2*Lx-dx)*rand(1,N); 
randy = -Ly + (2*Ly-dy)*rand(1,N);
% pick randomly-generated points which inside the cell
randr    = rand(1,N);
randf    = neighbor_nearest4( f,    randx, randy, n, dx, dy, Lx, Ly );
randphi0 = neighbor_nearest4( phi0, randx, randy, n, dx, dy, Lx, Ly );
box0_adh = (randr < randf & randphi0>0.5);

randx = randx(box0_adh); randx = randx(1:Nadh); randx0=randx; randxo=randx;
randy = randy(box0_adh); randy = randy(1:Nadh); randy0=randy; randyo=randy;
% initial age of adhesion sites
tadh = zeros(1,Nadh);
% initial index for gripping and slipping sites
label_grip =  true(1,Nadh); 
label_slip = false(1,Nadh);
% initial Fgrip, Fslip
Fslip1 = zeros(1,Nadh);
Fslip2 = zeros(1,Nadh);

% data saving
saveCount = 1;
saveFreq  = 1000;
image_directory = 'Figures';
save_file_dirct = 'Data'; %_name

if ~exist(image_directory)
	mkdir(image_directory);
end
if ~exist(save_file_dirct)
	mkdir(save_file_dirct);
end

Xcenter0 = sum(sum(x.*phi0))/sum(sum(phi0));
Ycenter0 = sum(sum(y.*phi0))/sum(sum(phi0));

Vcenter_rec = zeros(floor(T/(dt*mod_N)),3);

%% for-loop iteration
for step = 1 : T/dt
 
% STEP 1: pre-calculation 
  dG            = 36*phi0.*(1-phi0).*(1-2*phi0);  
  [dphix,dphiy] = grad(phi0,m,n,Lx,Ly);
  grad_phi_sqr  = dphix.^2+dphiy.^2;
  grad_phi_abs  = sqrt(grad_phi_sqr);  
  
% STEP 2: phase field function: phi 
   [curv_term,Fx,Fy]  =  force(phi0,dphix,dphiy,grad_phi_sqr,rhoa0,rhom0,m,n,dx,dy,Lx,Ly,etaa0,etam0,epsilon,gamma0,kappa0,Peri_critical); 
   Multiplier = Gamma*(-dG/epsilon + epsilon*curv_term);
   phi        = phi0 + dt*(-u0.*dphix - v0.*dphiy) + dt*Multiplier;  
   phi_box    = (phi>=lambda);
   
% STEP 3: rhoa and rhom 
  [rhoa, rhom] = reaction5(rhoa0,rhom0,phi0,phi,u0,v0,dt,dx,dy,m,n,Lx,Ly,lambda,rhoa_tot,ka,kb,kc,Ka,KD,Da,Dm0);
  
% STEP 4: velocity field: u and v  
  [u,v,Fgrip,randxn,randyn,Fgripx_mx,Fgripy_mx,Fslipx_mx,Fslipy_mx,tadh] = velocitysolver(Fx,Fy,tadh,dt,phi,m,n,dx,dy,u0,v0,randxo,randyo,Lx,Ly,kgrip0,kslip0,randx0,randy0,label_slip,label_grip,Nadh,nu0,trick,xi,phi_box,step,mod_N);
   
%STEP 6: periodic shift
[phi,rhoa,rhom,u,v,randx0,randy0,randxn,randyn] = periodicshift(phi,rhoa,rhom,u,v,randx0,randy0,randxn,randyn,dx,dy,m,n,Lx,Ly,lambda); 
% STEP 7: update data
   phi0  = phi; rhoa0 = rhoa; rhom0 = rhom; u0 = u; v0 = v;   
% STEP 8: update adhesion-related terms tadh, label_grip, label_slip
  [label_grip,label_slip,tadh,randxo,randyo,randx0,randy0,box_adh,Ndie] = adhesion_update(Fgrip,Fgrip0,Nadh,dt,label_grip,phi,rhoa,randx0,randy0,randxn,randyn,n,dx,dy,Lx,Ly,ron,roff0,rdie,N,tadh);
 
%  plotting data
     if (mod(step, mod_N)==0)  
         
         N_grip = sum(label_grip);
         N_slip = sum(label_slip);
         fprintf('Time %1.4f\v Grip %1.f\v Slip %1.f\v box_adn %1.f\v  Ndie %1.f\n',step*dt,N_grip,N_slip,sum(box_adh),Ndie);
         
         Xcenter = sum(sum(x.*phi))/sum(sum(phi));
         Ycenter = sum(sum(y.*phi))/sum(sum(phi));
         
         Vx_center = ( Xcenter - Xcenter0 ) / ( dt * mod_N );
         Vy_center = ( Ycenter - Ycenter0 ) / ( dt * mod_N );
         
         Xcenter0 = Xcenter;
         Ycenter0 = Ycenter;
         Vcenter_rec(step/mod_N,:) = [ Vx_center, Vy_center, sqrt(Vx_center * Vx_center +  Vy_center * Vy_center) ];

	 clf; color_plot_single(x,y,phi,rhoa,rhom,u,v,Fgripx_mx,Fgripy_mx,Fslipx_mx,Fslipy_mx,xi,Fgrip0); drawnow; 
          
          if (mod(step,saveFreq)==0)
            if(~isempty(image_directory))
               savename = sprintf('%s/%d.png',image_directory,saveCount);
               print('-dpng',savename);     
            end
            if(rem(saveCount,5)==0) 
               save_file_name = sprintf('Data/%d.mat',saveCount); 
               save(save_file_name);
            end
            saveCount = saveCount + 1;
          end 
        
                 
     end 
         
end

csvwrite(where2save, Vcenter_rec)
