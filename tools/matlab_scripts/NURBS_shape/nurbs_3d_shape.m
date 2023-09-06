% This is a function that we use to calculate the 3D NURBS basis
clear all; clc;

ex = 2; ey = 3; ez = 2;
nqp = 3;
qua_x = 1; qua_y = 2; qua_z = 3;

sDegree = 2; tDegree = 2; uDegree = 2;

U = [0.0, 0.0, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.0, 1.0];
V = [0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0];
W = [0.0, 0.0, 0.0, 0.333333, 0.666667, 1.0, 1.0, 1.0];

nElemX = length(U) - 2*sDegree - 1;
nElemY = length(V) - 2*tDegree - 1;
nElemZ = length(W) - 2*uDegree - 1;


e_x_start = U(sDegree + ex); e_x_end = U(sDegree + ex + 1);
e_y_start = V(tDegree + ey); e_y_end = V(tDegree + ey + 1);
e_z_start = W(uDegree + ez); e_z_end = W(uDegree + ez + 1);

eh_x = e_x_end - e_x_start;
eh_y = e_y_end - e_y_start;
eh_z = e_z_end - e_z_start;

[qp, wq] = Gauss(nqp, 0, 1);
qp = fliplr(qp');
wq = fliplr(wq');

qp_x = e_x_start + qp * eh_x;
qp_y = e_y_start + qp * eh_y;
qp_z = e_z_start + qp * eh_z;


qua = [qp_x(qua_x), qp_y(qua_y), qp_z(qua_z)];

% we calculate the (sDegree+1) * (tDegree+1) * (uDegree +1) functions
% within element (ex, ey, ez) at (qua_x, qua_y, qua_z).
nLocBas = (sDegree+1) * (tDegree+1) * (uDegree + 1);
R = zeros(nLocBas, 1);
dR_dx = zeros(nLocBas, 1);
dR_dy = zeros(nLocBas, 1);
dR_dz = zeros(nLocBas, 1);

for ii = 0 : sDegree
  index_x = ex + ii - 1;
  
  f_x = DerOneBasisFun(sDegree, 1, index_x, U, qua(1));
  
  for jj = 0 : tDegree
    index_y = ey + jj - 1;
    
    f_y = DerOneBasisFun(tDegree, 1, index_y, V, qua(2));
    
    for kk = 0 : uDegree
      index_z = ez + kk - 1;
      
      f_z = DerOneBasisFun(uDegree, 1, index_z, W, qua(3));
      
      index = ii + jj * (sDegree+1) + kk * (sDegree+1) * (tDegree + 1) + 1;
      
      R(index) = f_x(1) * f_y(1) * f_z(1);
      dR_dx(index) = f_x(2) * f_y(1) * f_z(1);
      dR_dy(index) = f_x(1) * f_y(2) * f_z(1);
      dR_dz(index) = f_x(1) * f_y(1) * f_z(2);      
    end
  end
end
