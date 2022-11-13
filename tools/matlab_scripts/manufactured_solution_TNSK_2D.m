% This is a code I compute manufactured solution forcing term for TNSK 2D
clear all; clc; 

syms x y t;
syms invRe invWe pi fac8_27 inv_gam fac4_3 fac2_3 kappa;

% ------- Exact Solution -------
rho = 0.6;

u = sin(pi*t)*sin(2*pi*x)*sin(2*pi*y);
v = 0.0;

theta = 0.85;
% ------------------------------

rho_x = diff(rho,x);
rho_y = diff(rho,y);

lap_rho = diff(rho_x, x) + diff(rho_y, y);
lap_rho = simple(lap_rho);


p = fac8_27 * theta * rho / (1-rho) - rho * rho;

iota = -rho + fac8_27 * theta * inv_gam;


rhoE = rho * iota + 0.5 * rho * (u*u + v*v) + 0.5 * invWe * (rho_x * rho_x + rho_y * rho_y);

rhoE = simple(rhoE);

% ------ f1 = mass eqn. -------
f1 = diff(rho,t) + diff(rho*u,x) + diff(rho*v, y);
f1 = simple(f1);
% -----------------------------

u_x = diff(u,x); u_x = simple(u_x);
v_x = diff(v,x); v_x = simple(v_x);
u_y = diff(u,y); u_y = simple(u_y);
v_y = diff(v,y); v_y = simple(v_y);

tau11 = invRe * (fac4_3 * u_x - fac2_3 * v_y); tau11 = simple(tau11);
tau12 = invRe * (u_y + v_x); tau12 = simple(tau12);
tau21 = tau12;
tau22 = invRe * (fac4_3 * v_y - fac2_3 * u_x); tau22 = simple(tau22);


sig11 = invWe * ( rho * lap_rho - 0.5 * rho_x * rho_x + 0.5 * rho_y * rho_y );
sig12 = -invWe * rho_x * rho_y;
sig21 = -invWe * rho_x * rho_y;
sig22 = invWe * ( rho * lap_rho + 0.5 * rho_x * rho_x - 0.5 * rho_y * rho_y );

sig11 = simple(sig11); sig12 = simple(sig12); sig21 = simple(sig21); sig22 = simple(sig22);

% ------- f2 = mom in x -------
f2 = diff(rho*u, t) + diff(rho*u*u, x) + diff(rho*u*v, y) + diff(p,x) ...
  - diff(tau11, x) - diff(tau12, y) - diff(sig11,x) - diff(sig12,y);
f2 = simple(f2);
% -----------------------------

% ------- f3 = mom in y -------
f3 = diff(rho*v, t) + diff(rho*u*v, x) + diff(rho*v*v, y) + diff(p,y) ...
  - diff(tau21, x) - diff(tau22, y) - diff(sig21, x) - diff(sig22, y);
f3 = simple(f3);
% -----------------------------

q_x = - kappa * diff(theta,x);
q_y = - kappa * diff(theta,y);

PI_x = invWe * rho * (u_x + v_y) * rho_x;
PI_y = invWe * rho * (u_x + v_y) * rho_y;


% ------- f4 = energy eqn. -------
f4 = diff(rhoE, t) + diff(rhoE*u, x) + diff(rhoE*v, y) + diff(p*u,x) + diff(p*v, y) ...
  - diff(tau11*u,x) - diff(tau12*v, x) - diff(tau21*u, y) - diff(tau22*v, y) ...
  - diff(sig11*u,x) - diff(sig12*v, x) - diff(sig21*u, y) - diff(sig22*v, y) ...
  + diff(q_x, x) + diff(q_y, y) + diff(PI_x , x) + diff(PI_y , y);

f4 = simple(f4);
% --------------------------------

syms a b c d e f g;
f1 = subs(f1, sin(2*pi*x), a);
f1 = subs(f1, sin(2*pi*y), b);
f1 = subs(f1, sin(pi*t), c);
f1 = subs(f1, cos(2*pi*x), e);
f1 = subs(f1, cos(2*pi*y), f);
f1 = subs(f1, cos(pi*t), g);

f2 = subs(f2, sin(2*pi*x), a);
f2 = subs(f2, sin(2*pi*y), b);
f2 = subs(f2, sin(pi*t), c);
f2 = subs(f2, cos(2*pi*x), e);
f2 = subs(f2, cos(2*pi*y), f);
f2 = subs(f2, cos(pi*t), g);

f3 = subs(f3, sin(2*pi*x), a);
f3 = subs(f3, sin(2*pi*y), b);
f3 = subs(f3, sin(pi*t), c);
f3 = subs(f3, cos(2*pi*x), e);
f3 = subs(f3, cos(2*pi*y), f);
f3 = subs(f3, cos(pi*t), g);

f4 = subs(f4, sin(2*pi*x), a);
f4 = subs(f4, sin(2*pi*y), b);
f4 = subs(f4, sin(pi*t), c);
f4 = subs(f4, cos(2*pi*x), e);
f4 = subs(f4, cos(2*pi*y), f);
f4 = subs(f4, cos(pi*t), g);

syms pi2;
f1 = subs(f1, pi^2, pi2);
f2 = subs(f2, pi^2, pi2);
f3 = subs(f3, pi^2, pi2);
f4 = subs(f4, pi^2, pi2);

f1 = simple(f1);
f2 = simple(f2);
f3 = simple(f3);
f4 = simple(f4);