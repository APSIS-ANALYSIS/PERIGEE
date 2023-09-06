% This is the code that compute the manufactured solution forcing for 
% Hyperelastic 3D model
clear all; clc;

syms x y z t;
syms kappa mu rho0;
syms pi;
syms J0d67;

u = sin(pi*t)*sin(3*pi*x) * sin(3*pi*y) * sin(3*pi*z);
v = 0; %sin(pi*t) * sin(pi*x) * sin(pi*y) * sin(pi*z);
w = 0; %sin(pi*t) * sin(pi*x) * sin(pi*y) * sin(pi*z);

ut = diff(u, t); utt = diff(ut, t); utt = simplify(utt);
vt = diff(v, t); vtt = diff(vt, t); vtt = simplify(vtt);
wt = diff(w, t); wtt = diff(wt, t); wtt = simplify(wtt);

F = [diff(u,x), diff(u,y), diff(u,z);
     diff(v,x), diff(v,y), diff(v,z);
     diff(w,x), diff(w,y), diff(w,z)];
 
I = [1, 0, 0; 
    0, 1, 0; 
    0, 0, 1];

F = F + I;

F = simplify(F);

J = det(F); J = simplify(J);

Jm0d67 = J^(-2/3);

Ft = transpose(F);

C = Ft * F; C = simplify(C);

trC = trace(C);

Cinv = inv(C); Cinv = simplify(Cinv);

S = 0.5 * kappa * (J*J - 1.0) * Cinv + mu * J0d67 * I ...
    - mu * J0d67 * trC * (1/3) * Cinv;

P = F * S;

fx = rho0 * utt  - diff(P(1,1), x) - diff(P(1,2), y) - diff(P(1,3), z);
fx = simplify(fx);
fy = rho0 * vtt  - diff(P(2,1), x) - diff(P(2,2), y) - diff(P(2,3), z);
fy = simplify(fy);
fz = rho0 * wtt  - diff(P(3,1), x) - diff(P(3,2), y) - diff(P(3,3), z);
fz = simplify(fz);


% EOF