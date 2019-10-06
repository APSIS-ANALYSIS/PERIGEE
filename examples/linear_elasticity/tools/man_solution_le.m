% This is the code that compute the manufactured solution forcing for 
% Hyperelastic 3D model
clear all; clc;

syms x y z;
syms lambda mu;

u = sin(pi*x) * sin(2*pi*y) * sin(pi*z);
v = sin(2*pi*x) * sin(pi*y) * sin(pi*z);
w = sin(pi*x) * sin(pi*y) * sin(2*pi*z);

ux = diff(u,x); ux = simplify(ux);
uy = diff(u,y); uy = simplify(uy);
uz = diff(u,z); uz = simplify(uz);

vx = diff(v,x); vx = simplify(vx);
vy = diff(v,y); vy = simplify(vy);
vz = diff(v,z); vz = simplify(vz);

wx = diff(w,x); wx = simplify(wx);
wy = diff(w,y); wy = simplify(wy);
wz = diff(w,z); wz = simplify(wz);

divu = ux + vy + wz;

%sigma = zeros(3,3);

sigma(1,1) = 2 * mu * ux + lambda * divu;
sigma(1,2) = mu * (uy + vx);
sigma(1,3) = mu * (uz + wx);

sigma(2,1) = sigma(1,2);
sigma(2,2) = 2 * mu * vy + lambda * divu;
sigma(2,3) = mu * (vz + wy);

sigma(3,1) = sigma(1,3);
sigma(3,2) = sigma(2,3);
sigma(3,3) = 2 * mu * wz + lambda * divu;

fx = - diff(sigma(1,1), x) - diff(sigma(1,2), y) - diff(sigma(1,3), z);
fx = simplify(fx);
fy = - diff(sigma(2,1), x) - diff(sigma(2,2), y) - diff(sigma(2,3), z);
fy = simplify(fy);
fz = - diff(sigma(3,1), x) - diff(sigma(3,2), y) - diff(sigma(3,3), z);
fz = simplify(fz);


% Get the sigma l2 norm
sigma2 = 0.0;
for ii = 1 : 3
  for jj = 1 : 3
    sigma2 = sigma2 + sigma(ii,jj) * sigma(ii,jj);
  end
end

sigma2x = int(sigma2, x, 0, 1);
sigma2xy = int(sigma2x, y, 0, 1);
sigma2xyz = int(sigma2xy, z, 0, 1); 

% Get the tractions on the boundary
%H_top = sigma * [0;0;1];  H_top = subs(H_top, z, 1.0); H_top = simplify(H_top);
H_fro = sigma * [1;0;0];  H_fro = subs(H_fro, x, 1.0); H_fro = simplify(H_fro);
H_bac = sigma * [-1;0;0]; H_bac = subs(H_bac, x, 0.0); H_bac = simplify(H_bac);
H_rig = sigma * [0;1;0];  H_rig = subs(H_rig, y, 1.0); H_rig = simplify(H_rig);
H_lef = sigma * [0;-1;0]; H_lef = subs(H_lef, y, 0.0); H_lef = simplify(H_lef);

% EOF