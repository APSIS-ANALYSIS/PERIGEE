clear all;
clc;

syms x y z tt lambda mu rho aa w;

ux = aa*x*sin(w*x)*sin(w*y)*sin(w*z)*tt;
uy = aa*y*sin(w*x)*sin(w*y)*sin(w*z)*tt;
uz = aa*z*sin(w*x)*sin(w*y)*sin(w*z)*tt;

ux_x = diff(ux,x);
ux_y = diff(ux,y);
ux_z = diff(ux,z);

uy_x = diff(uy,x);
uy_y = diff(uy,y);
uy_z = diff(uy,z);

uz_x = diff(uz,x);
uz_y = diff(uz,y);
uz_z = diff(uz,z);

ux_tt = diff(ux,tt,2);
uy_tt = diff(uy,tt,2);
uz_tt = diff(uz,tt,2);

ep_xx = ux_x;
ep_yy = uy_y;
ep_zz = uz_z;
ep_yz = uy_z+uz_y;
ep_xz = ux_z+uz_x;
ep_xy = ux_y+uy_x;

sig_xx = (lambda+2*mu)*ep_xx + lambda*(ep_yy+ep_zz);
sig_yy = (lambda+2*mu)*ep_yy + lambda*(ep_xx+ep_zz);
sig_zz = (lambda+2*mu)*ep_zz + lambda*(ep_xx+ep_yy);
sig_yz = mu*ep_yz;
sig_xz = mu*ep_xz;
sig_xy = mu*ep_xy;

fx = rho*ux_tt-diff(sig_xx,x)-diff(sig_xy,y)-diff(sig_xz,z);
fy = rho*uy_tt-diff(sig_xy,x)-diff(sig_yy,y)-diff(sig_yz,z);
fz = rho*uz_tt-diff(sig_xz,x)-diff(sig_yz,y)-diff(sig_zz,z);

fx = simplify(fx);
fy = simplify(fy);
fz = simplify(fz);
