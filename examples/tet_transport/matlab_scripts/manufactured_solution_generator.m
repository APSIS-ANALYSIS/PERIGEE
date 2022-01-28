clear all; clc;

syms x y z t rho cap kappa;

u = t*t*t*t * x*(x-1)*y*(y-1)*z*(z-1);

u_x = diff(u,x);
u_y = diff(u,y);
u_z = diff(u,z);

f = rho * cap * diff(u,t) - kappa * diff(u_x, x) - kappa * diff(u_y, y) - kappa * diff(u_z, z);

f = simplify(f);