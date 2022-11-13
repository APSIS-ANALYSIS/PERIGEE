% -------------------------------------------------------------
% NURBS-3D shape routine
% This is an evaluation routine that gives the volumetric
% NURBS basis function and its derivative values.
% June 15 2015 Ju Liu
% -------------------------------------------------------------
clear all; clc;

% Identify the element number we want to evaluate in each direction
ex = 3; ey = 4; ez = 3;

nqp = 4; % number of evaluation points in each direction

qua_x = 1; qua_y = 4; qua_z = 2; % evaluation point id in each direction

quaindex = qua_x - 1 + (qua_y-1)*nqp + (qua_z-1)*nqp*nqp;

sDegree = 1; tDegree = 2; uDegree = 2;

nLocBas = (sDegree+1) * (tDegree+1) * (uDegree + 1);

% Specify the knot vectors
U = [0.0, 0.0, 1/3, 2/3, 1.0, 1.0];
V = [0.0, 0.0, 0.0, 1/6, 1/3, 1/2, 2/3, 5/6, 1.0, 1.0, 1.0];
W = [0.0, 0.0, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.0, 1.0];

% Specify the weights at the element (ex, ey, ez)
Weight = dlmread('cpw.txt');

% Specify the control points for basis functions in the element
cp_x = dlmread('cpx.txt');
cp_y = dlmread('cpy.txt');
cp_z = dlmread('cpz.txt');


% Generate the number of elements in each direction
nElemX = length(U) - 2*sDegree - 1;
nElemY = length(V) - 2*tDegree - 1;
nElemZ = length(W) - 2*uDegree - 1;

e_x_start = U(sDegree + ex); e_x_end = U(sDegree + ex + 1);
e_y_start = V(tDegree + ey); e_y_end = V(tDegree + ey + 1);
e_z_start = W(uDegree + ez); e_z_end = W(uDegree + ez + 1);

eh_x = e_x_end - e_x_start;
eh_y = e_y_end - e_y_start;
eh_z = e_z_end - e_z_start;

eindex = ex-1 + (ey-1) * nElemX + (ez-1) * nElemX * nElemY;

% Generate the evaluation points
[qp, wq] = Gauss(nqp, 0, 1);

qp_x = e_x_start + 0.999999 * eh_x;
qp_y = e_y_start + qp * eh_y;
qp_z = e_z_start + qp* eh_z;

qua = [qp_x(qua_x), qp_y(qua_y), qp_z(qua_z)];

% Now we calculate the B-spline basis functions within the
% element (ex, ey, ez) at (qua(1), qua(2), qua(3)).

N_x = zeros(sDegree+1, 2); % the local basis function evaluation
N_y = N_x; N_z = N_x;

for ii = 0 : sDegree
  index_x = ex + ii - 1;
  
  N_x(ii+1,:) = DerOneBasisFun(sDegree, 1, index_x, U, qua(1));
end

for jj = 0 : tDegree
  index_y = ey + jj - 1;
  
  N_y(jj+1,:) = DerOneBasisFun(tDegree, 1, index_y, V, qua(2));
end

for kk = 0 : uDegree
  index_z = ez + kk - 1;
  
  N_z(kk+1,:) = DerOneBasisFun(uDegree, 1, index_z, W, qua(3));
end
 
N = zeros(nLocBas,1);
dN_ds = N; dN_dt = N; dN_du = N;
w = 0.0; dw_ds = 0.0; dw_dt = 0.0; dw_du = 0.0;

for kk = 0 : uDegree
  for jj = 0 : tDegree
    for ii = 0 : sDegree
      index = ii + jj * (sDegree+1) ...
        + kk * (sDegree+1) * (tDegree + 1) + 1;
      N(index) = N_x(ii+1,1) * N_y(jj+1,1) * N_z(kk+1,1) * Weight(index);
      w = w + N(index);
      
      dN_ds(index) = N_x(ii+1,2) * N_y(jj+1,1) * N_z(kk+1,1) * Weight(index);
      dN_dt(index) = N_x(ii+1,1) * N_y(jj+1,2) * N_z(kk+1,1) * Weight(index);
      dN_du(index) = N_x(ii+1,1) * N_y(jj+1,1) * N_z(kk+1,2) * Weight(index);
      
      dw_ds = dw_ds + dN_ds(index);
      dw_dt = dw_dt + dN_dt(index);
      dw_du = dw_du + dN_du(index);
    end
  end
end

R = zeros(nLocBas, 1);
dR_ds = R; dR_dt = R; dR_du = R;

for index = 1 : nLocBas
  R(index) = N(index) / w;
  dR_ds(index) = (dN_ds(index) - R(index) * dw_ds) / w;
  dR_dt(index) = (dN_dt(index) - R(index) * dw_dt) / w;
  dR_du(index) = (dN_du(index) - R(index) * dw_du) / w;
end

dxds = zeros(3,3);

dxds(1,1) = cp_x' * dR_ds; 
dxds(1,2) = cp_x' * dR_dt; 
dxds(1,3) = cp_x' * dR_du;
dxds(2,1) = cp_y' * dR_ds; 
dxds(2,2) = cp_y' * dR_dt; 
dxds(2,3) = cp_y' * dR_du;
dxds(3,1) = cp_z' * dR_ds; 
dxds(3,2) = cp_z' * dR_dt; 
dxds(3,3) = cp_z' * dR_du;

dsdx = inv(dxds);

dR_dx = dR_ds * dsdx(1,1) + dR_dt * dsdx(2,1) + dR_du * dsdx(3,1);
dR_dy = dR_ds * dsdx(1,2) + dR_dt * dsdx(2,2) + dR_du * dsdx(3,2);
dR_dz = dR_ds * dsdx(1,3) + dR_dt * dsdx(2,3) + dR_du * dsdx(3,3);

Jacobi = det(dxds);

% EOF