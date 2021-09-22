function [label_grip,label_slip,tadh,randxo,randyo,randx0,randy0,box_adh,Ndie] = adhesion_update(Fgrip,Fgrip0,Nadh,dt,label_grip,phi,rhoa,randx0,randy0,randxn,randyn,n,dx,dy,Lx,Ly,ron,roff0,rdie,N,tadh)

%      Fgrip0_rand = neighbor_nearest4(Fgrip0,randxn,randyn,n, dx, dy, Lx, Ly);
%      roff  = roff0*exp(Fgrip./Fgrip0_rand);     
     roff  = roff0*exp(Fgrip/Fgrip0); 
     
     roff1 = rand(1,Nadh);
     label_grip(roff1<roff*dt) = 0; 
     label_slip = 1 - label_grip;
     
     % For checking r_on and r_die, it starts with the same random sequence
     % ron1.
     ron1 = rand(1,Nadh);
     label_slip(ron1<ron*dt) = 0;
     label_grip = 1 - label_slip;
     
     label_slip = logical(label_slip);
     label_grip = logical(label_grip);
        
     randphi  = neighbor_nearest4(phi,randxn,randyn,n,dx,dy,Lx,Ly);
     Box_die1 = (label_slip & ron1<(ron+rdie)*dt);
     Box_die2 = (randphi<0.5);
     Box_die3 = (randxn<-Lx | randxn>Lx);
     Box_die  = (Box_die1 | Box_die2 | Box_die3);
     Ndie = sum(Box_die); Ndie1 = N;
         
     randx_reborn = -Lx + (2*Lx-dx)*rand(1,Ndie1);  
     randy_reborn = -Ly + (2*Ly-dy)*rand(1,Ndie1);
     
     randr_reborn = rand(1,Ndie1);
     phirhoa = phi.*rhoa;
     phirhoa_max = max(phirhoa(:));
     f_reborn = phirhoa/phirhoa_max;
     randf_reborn   = neighbor_nearest4(f_reborn,randx_reborn,randy_reborn,n,dx,dy,Lx,Ly);
     randphi_reborn = neighbor_nearest4(phi,     randx_reborn,randy_reborn,n,dx,dy,Lx,Ly);

     box_adh = (randr_reborn < randf_reborn & randphi_reborn>0.5);
     randx1_reborn = randx_reborn(box_adh); 
     randy1_reborn = randy_reborn(box_adh);
     randx1_reborn = randx1_reborn(1: Ndie); 
     randy1_reborn = randy1_reborn(1: Ndie);    
     
     tadh(Box_die) = 0;
     label_grip(Box_die) = 1; 
     label_slip = 1 - label_grip;
     label_grip = logical(label_grip);
     label_slip = logical(label_slip);
     randxn(Box_die) = randx1_reborn;
     randyn(Box_die) = randy1_reborn;
     randx0(Box_die) = randx1_reborn;
     randy0(Box_die) = randy1_reborn;
     randxo = randxn;
     randyo = randyn; 