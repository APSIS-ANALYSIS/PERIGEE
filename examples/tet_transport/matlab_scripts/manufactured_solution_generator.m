clear all; clc;

syms x y z t rho cap kappa;

u = t*t * sin(pi*x) * sin(pi*y) * sin(pi*z);

u_x = diff(u,x);
u_y = diff(u,y);
u_z = diff(u,z);

f = rho * cap * diff(u,t) - kappa * diff(u_x, x) - kappa * diff(u_y, y) - kappa * diff(u_z, z);

f = simplify(f);