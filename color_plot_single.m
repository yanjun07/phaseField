function h = color_plot_single(x,y,phi,rhoa,rhom,u,v,Fgripx_mx,Fgripy_mx,Fslipx_mx,Fslipy_mx,xi,Fgrip0)

colormap(summer);

Fgrip = sqrt(Fgripx_mx.^2+Fgripy_mx.^2);
Fslip = sqrt(Fslipx_mx.^2+Fslipy_mx.^2);
Ffric = xi*sqrt(u.^2+v.^2);
Ftot  = Fgrip + Fslip + Ffric;

imagesc([min(x(:)) max(x(:))],[min(y(:)) max(y(:))],Ftot/100,[0,1.5]);
shading('interp');
axis image;
axis off;
colorbar('SouthOutside');
hold on;

contour(x,y,phi,0.5*[1,1],'linewidth',1.5,'color',[67,190,216]/256);hold on;

phirhoa = phi.*rhoa; phirhoamax = max(phirhoa(:));
contour(x,y,phirhoa,phirhoamax/2*[1,1],'--','linewidth',3,'color',[247,168,174]/256);

phi_cut = zeros(size(phi));
phi_cut(phi>=1/2) = 1;
phi_cut(phi<1/2) = 0;
quiver(x(1:10:end,1:10:end),y(1:10:end,1:10:end),...
       phi_cut(1:10:end,1:10:end).*u(1:10:end,1:10:end),...
       phi_cut(1:10:end,1:10:end).*v(1:10:end,1:10:end),'-w');
   

h = 0;


