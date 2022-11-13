% Obtain manufactured solution force for the incompressible Navier-Stokes
% Equations.
clear all; clc;

syms x y z t pi nu;

phi_1 = t*t*x*(x-1)*y^2*(y-1)^2*z^2*(z-1)^2;
phi_3 = t*t*x^2*(x-1)^2*y^2*(y-1)^2*z*(z-1);

u = diff(phi_3, y);

v = diff(phi_1, z) - diff(phi_3,x);

w = -diff(phi_1, y);

u = simplify(u);

v = simplify(v);

w = simplify(w);

div_vel = diff(u,x) + diff(v,y) + diff(w,z);
div_vel = simplify(div_vel);

p = t*t*sin(pi*x) * sin(pi*y)*sin(pi*z);


u_t = diff(u,t);
u_x = diff(u,x);
u_y = diff(u,y);
u_z = diff(u,z);
u_xx = diff(u_x,x);
u_yy = diff(u_y,y);
u_zz = diff(u_z,z);

v_t = diff(v,t);
v_x = diff(v,x);
v_y = diff(v,y);
v_z = diff(v,z);
v_xx = diff(v_x,x);
v_yy = diff(v_y,y);
v_zz = diff(v_z,z);

w_t = diff(w,t);
w_x = diff(w,x);
w_y = diff(w,y);
w_z = diff(w,z);
w_xx = diff(w_x,x);
w_yy = diff(w_y,y);
w_zz = diff(w_z,z);

p_x = diff(p,x);
p_y = diff(p,y);
p_z = diff(p,z);


f1 = u_t + u_x * u + u_y * v + u_z * w + p_x - nu * (u_xx+u_yy+u_zz);

f2 = v_t + v_x * u + v_y * v + v_z * w + p_y - nu * (v_xx+v_yy+v_zz);

f3 = w_t + w_x * u + w_y * v + w_z * w + p_z - nu * (w_xx+w_yy+w_zz);

f1 = simplify(f1);

f2 = simplify(f2);

f3 = simplify(f3);

% EOF