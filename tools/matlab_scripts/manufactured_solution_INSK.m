% This is a symbolic computing code to get manufactured solution for INSK 3D

clear all; clc;

syms x y z t;

syms invRe Ca pi theta;

rho = 0.6 + 0.1*t*cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z);

u = t^3 * sin(2*pi*x) * sin(pi*y) * sin(pi*z);

v = t^3 * sin(pi*x) * sin(2*pi*y) * sin(pi*z);

w = sin(pi*t) * sin(pi*x) * sin(pi*y) * sin(2*pi*z);

p = 8*theta*rho/(27*(1-rho)) - rho * rho;

div_velo = diff(u,x) + diff(v,y) + diff(w,z);

tau11 = (2*diff(u,x) - 2*div_velo / 3) * invRe;
tau12 = (diff(u,y) + diff(v,x)) * invRe;
tau13 = ( diff(u,z) + diff(w,x) ) * invRe;

tau21 = tau12;
tau22 = ( 2*diff(v,y) - 2*div_velo / 3 ) * invRe;
tau23 = ( diff(v,z) + diff(w,y) ) * invRe;

tau31 = tau13;
tau32 = tau23;
tau33 = ( 2*diff(w,z) - 2*div_velo / 3 ) * invRe;

lap_rho = diff(diff(rho, x), x) + diff(diff(rho,y),y) + diff(diff(rho, z), z);

% derivatives
f1 = diff(rho, t) + diff(rho*u, x) + diff(rho*v, y) + diff(rho*w, z);

f2 = diff(rho*u , t) + diff(rho*u*u, x) + diff(rho*u*v, y) + diff(rho*u*w , z) ...
  + diff(p, x) - diff(tau11, x) - diff(tau12, y) - diff(tau13, z) - Ca*Ca* rho * diff(lap_rho, x);

f3 = diff(rho*v , t) + diff(rho*u*v, x) + diff(rho*v*v, y) + diff(rho*v*w , z) ...
  + diff(p, y) - diff(tau21, x) - diff(tau22, y) - diff(tau23, z) - Ca*Ca* rho * diff(lap_rho, y);

f4 = diff(rho*w , t) + diff(rho*u*w, x) + diff(rho*w*v, y) + diff(rho*w*w , z) ...
  + diff(p, z) - diff(tau31, x) - diff(tau32, y) - diff(tau33, z) - Ca*Ca* rho * diff(lap_rho, z);

syms a b c d e f;
f1 = subs(f1, sin(2*pi*x), a);
f1 = subs(f1, sin(2*pi*y), b);
f1 = subs(f1, sin(2*pi*z), c);

f1 = subs(f1, cos(2*pi*x), d);
f1 = subs(f1, cos(2*pi*y), e);
f1 = subs(f1, cos(2*pi*z), f);

f2 = subs(f2, sin(2*pi*x), a);
f2 = subs(f2, sin(2*pi*y), b);
f2 = subs(f2, sin(2*pi*z), c);

f2 = subs(f2, cos(2*pi*x), d);
f2 = subs(f2, cos(2*pi*y), e);
f2 = subs(f2, cos(2*pi*z), f);

f3 = subs(f3, sin(2*pi*x), a);
f3 = subs(f3, sin(2*pi*y), b);
f3 = subs(f3, sin(2*pi*z), c);

f3 = subs(f3, cos(2*pi*x), d);
f3 = subs(f3, cos(2*pi*y), e);
f3 = subs(f3, cos(2*pi*z), f);

f4 = subs(f4, sin(2*pi*x), a);
f4 = subs(f4, sin(2*pi*y), b);
f4 = subs(f4, sin(2*pi*z), c);

f4 = subs(f4, cos(2*pi*x), d);
f4 = subs(f4, cos(2*pi*y), e);
f4 = subs(f4, cos(2*pi*z), f);

syms g h k l m n;
f1 = subs(f1, sin(pi*x), g);
f1 = subs(f1, sin(pi*y), h);
f1 = subs(f1, sin(pi*z), k);

f1 = subs(f1, cos(pi*x), l);
f1 = subs(f1, cos(pi*y), m);
f1 = subs(f1, cos(pi*z), n);

f2 = subs(f2, sin(pi*x), g);
f2 = subs(f2, sin(pi*y), h);
f2 = subs(f2, sin(pi*z), k);

f2 = subs(f2, cos(pi*x), l);
f2 = subs(f2, cos(pi*y), m);
f2 = subs(f2, cos(pi*z), n);

f3 = subs(f3, sin(pi*x), g);
f3 = subs(f3, sin(pi*y), h);
f3 = subs(f3, sin(pi*z), k);

f3 = subs(f3, cos(pi*x), l);
f3 = subs(f3, cos(pi*y), m);
f3 = subs(f3, cos(pi*z), n);

f4 = subs(f4, sin(pi*x), g);
f4 = subs(f4, sin(pi*y), h);
f4 = subs(f4, sin(pi*z), k);

f4 = subs(f4, cos(pi*x), l);
f4 = subs(f4, cos(pi*y), m);
f4 = subs(f4, cos(pi*z), n);


f1 = simple(f1);
f2 = simple(f2);
f3 = simple(f3);
f4 = simple(f4);