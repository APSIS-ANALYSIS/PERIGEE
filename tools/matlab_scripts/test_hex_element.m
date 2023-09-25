clear;
clc;
format long;

u = [1.0, 1.0, 1.0];
r = [0.0, 0.5, 1.0];
s = [0.0, 0.5, 1.0];
t = [0.0, 0.5, 1.0];

ctrl_x = genMap(r,u,u);
ctrl_y = genMap(u,s,u);
ctrl_z = genMap(u,u,s);

fid = fopen('ctrlpts.txt','wt');
rand_ratio = 0.1;
for ii = 1 : 27
    randx = rand()*rand_ratio-0.5*rand_ratio;
    randy = rand()*rand_ratio-0.5*rand_ratio;
    randz = rand()*rand_ratio-0.5*rand_ratio;
    ctrl_x(ii) = ctrl_x(ii) + randx;
    ctrl_y(ii) = ctrl_y(ii) + randy;
    ctrl_z(ii) = ctrl_z(ii) + randz;
end

for ii = 1 : 27
    fprintf(fid,'%.16f ', ctrl_x(ii)');
end
fprintf(fid,'\n');
for ii = 1 : 27
    fprintf(fid,'%.16f ', ctrl_y(ii)');
end
fprintf(fid,'\n');
for ii = 1 : 27
    fprintf(fid,'%.16f ', ctrl_z(ii)');
end
fprintf(fid,'\n');

qua_x = 0.5;
qua_y = 0.5;
qua_z = 0.5;

shape_0 = @(x) (2.0*x-1.0).*(x-1.0);
shape_1 = @(x) -4.0*x*(x-1.0);
shape_2 = @(x) x*(2.0*x-1.0);

shape_r = [shape_0(qua_x),shape_1(qua_x),shape_2(qua_x)];
shape_s = [shape_0(qua_y),shape_1(qua_y),shape_2(qua_y)];
shape_t = [shape_0(qua_z),shape_1(qua_z),shape_2(qua_z)];

R = genMap(shape_r,shape_s,shape_t);
for ii = 1 : 27
    fprintf(fid,'%.16f ', R(ii));
end
fprintf(fid,'\n');
dshape_0 = @(x) 4.0*x-3.0;
dshape_1 = @(x) -8.0*x+4.0;
dshape_2 = @(x) 4.0*x-1.0;

dshape_r = [dshape_0(qua_x),dshape_1(qua_x),dshape_2(qua_x)];
dshape_s = [dshape_0(qua_y),dshape_1(qua_y),dshape_2(qua_y)];
dshape_t = [dshape_0(qua_z),dshape_1(qua_z),dshape_2(qua_z)];

dR_dr = genMap(dshape_r,shape_s,shape_t);
dR_ds = genMap(shape_r,dshape_s,shape_t);
dR_dt = genMap(shape_r,shape_s,dshape_t);

d2R_drs = genMap(dshape_r,dshape_s,shape_t);
d2R_drt = genMap(dshape_r,shape_s,dshape_t);
d2R_dst = genMap(shape_r,dshape_s,dshape_t);

ddshape = [4.0, -8.0, 4.0];

d2R_drr = genMap(ddshape,shape_s,shape_t);
d2R_dss = genMap(shape_r,ddshape,shape_t);
d2R_dtt = genMap(shape_r,shape_s,ddshape);

xr = dot(ctrl_x,dR_dr);
xs = dot(ctrl_x,dR_ds);
xt = dot(ctrl_x,dR_dt);
yr = dot(ctrl_y,dR_dr);
ys = dot(ctrl_y,dR_ds);
yt = dot(ctrl_y,dR_dt);
zr = dot(ctrl_z,dR_dr);
zs = dot(ctrl_z,dR_ds);
zt = dot(ctrl_z,dR_dt);

xrr = dot(ctrl_x,d2R_drr);
xss = dot(ctrl_x,d2R_dss);
xtt = dot(ctrl_x,d2R_dtt);
yrr = dot(ctrl_y,d2R_drr);
yss = dot(ctrl_y,d2R_dss);
ytt = dot(ctrl_y,d2R_dtt);
zrr = dot(ctrl_z,d2R_drr);
zss = dot(ctrl_z,d2R_dss);
ztt = dot(ctrl_z,d2R_dtt);

xrs = dot(ctrl_x,d2R_drs);
xrt = dot(ctrl_x,d2R_drt);
xst = dot(ctrl_x,d2R_dst);
yrs = dot(ctrl_y,d2R_drs);
yrt = dot(ctrl_y,d2R_drt);
yst = dot(ctrl_y,d2R_dst);
zrs = dot(ctrl_z,d2R_drs);
zrt = dot(ctrl_z,d2R_drt);
zst = dot(ctrl_z,d2R_dst);

J = [xr, xs, xt; yr, ys, yt; zr, zs, zt];
D = det(J);
Jinv = inv(J);
JinvT = transpose(Jinv);
dR_dx = zeros(1,27);
dR_dy = zeros(1,27);
dR_dz = zeros(1,27);

dR_drst = [dR_dr; dR_ds; dR_dt];
dR_dxyz = JinvT * dR_drst;

dR_dx = dR_dxyz(1,:);
for ii = 1 : 27
    fprintf(fid,'%.16f ', dR_dx(ii));
end
fprintf(fid,'\n');
dR_dy = dR_dxyz(2,:);
for ii = 1 : 27
    fprintf(fid,'%.16f ', dR_dy(ii));
end
fprintf(fid,'\n');
dR_dz = dR_dxyz(3,:);
for ii = 1 : 27
    fprintf(fid,'%.16f ', dR_dz(ii));
end
fprintf(fid,'\n');

LHS = genLHS(xr,xs,xt,yr,ys,yt,zr,zs,zt);

d2R_dxx = zeros(1,27);
d2R_dyy = zeros(1,27);
d2R_dzz = zeros(1,27);
d2R_dxy = zeros(1,27);
d2R_dxz = zeros(1,27);
d2R_dyz = zeros(1,27);

for ii = 1 : 27
    temp = zeros(6,1);
    f = zeros(6,1);
    f(1) = d2R_drr(ii) - dR_dx(ii) * xrr - dR_dy(ii) * yrr - dR_dz(ii) * zrr;
    f(2) = d2R_dss(ii) - dR_dx(ii) * xss - dR_dy(ii) * yss - dR_dz(ii) * zss;
    f(3) = d2R_dtt(ii) - dR_dx(ii) * xtt - dR_dy(ii) * ytt - dR_dz(ii) * ztt;
    f(4) = d2R_drs(ii) - dR_dx(ii) * xrs - dR_dy(ii) * yrs - dR_dz(ii) * zrs;
    f(5) = d2R_drt(ii) - dR_dx(ii) * xrt - dR_dy(ii) * yrt - dR_dz(ii) * zrt;
    f(6) = d2R_dst(ii) - dR_dx(ii) * xst - dR_dy(ii) * yst - dR_dz(ii) * zst;
    temp = LHS \ f;
    d2R_dxx(ii)=temp(1);
    d2R_dyy(ii)=temp(2);
    d2R_dzz(ii)=temp(3);
    d2R_dxy(ii)=temp(4);
    d2R_dxz(ii)=temp(5);
    d2R_dyz(ii)=temp(6);
end
for ii = 1 : 27
    fprintf(fid,'%.16f ', d2R_dxx(ii));
end
fprintf(fid,'\n');
for ii = 1 : 27
    fprintf(fid,'%.16f ', d2R_dyy(ii));
end
fprintf(fid,'\n');
for ii = 1 : 27
    fprintf(fid,'%.16f ', d2R_dzz(ii));
end
fprintf(fid,'\n');
for ii = 1 : 27
    fprintf(fid,'%.16f ', d2R_dxy(ii));
end
fprintf(fid,'\n');
for ii = 1 : 27
    fprintf(fid,'%.16f ', d2R_dxz(ii));
end
fprintf(fid,'\n');
for ii = 1 : 27
    fprintf(fid,'%.16f ', d2R_dyz(ii));
end
fprintf(fid,'\n');
for ii = 1 : 3
    for jj = 1 : 3
        fprintf(fid,'%.16f ', J(ii,jj));
    end
end
fprintf(fid,'\n');
for ii = 1 : 3
    for jj = 1 : 3
        fprintf(fid,'%.16f ', Jinv(ii,jj));
    end
end
fprintf(fid,'\n');
fprintf(fid,'%.16f', D);
figure(1)
scatter3(ctrl_x,ctrl_y,ctrl_z);
for i = 1:27
   text(ctrl_x(i)+0.02,ctrl_y(i)+0.02,ctrl_z(i)+0.02,num2str(i-1));
end
grid on
function map = genMap(fr, fs, ft)
  map = zeros(1,27);
  temp = zeros(1,27);
  for kk = 1 : 3
      for jj = 1 : 3
          for ii = 1 : 3
             temp((kk-1)*9+(jj-1)*3+ii) = fr(ii)*fs(jj)*ft(kk);
          end
      end
  end
  map(1) = temp(1);
  map(2) = temp(3);
  map(3) = temp(9);
  map(4) = temp(7);
  map(5) = temp(19);
  map(6) = temp(21);
  map(7) = temp(27);
  map(8) = temp(25);
  map(9) = temp(2);
  map(10) = temp(6);
  map(11) = temp(8);
  map(12) = temp(4);
  map(13) = temp(20);
  map(14) = temp(24);
  map(15) = temp(26);
  map(16) = temp(22);
  map(17) = temp(10);
  map(18) = temp(12);
  map(19) = temp(18);
  map(20) = temp(16);
  map(21) = temp(13);
  map(22) = temp(15);
  map(23) = temp(11);
  map(24) = temp(17);
  map(25) = temp(5);
  map(26) = temp(23);
  map(27) = temp(14);
end

function mat = genLHS(xr,xs,xt,yr,ys,yt,zr,zs,zt)
  mat = zeros(6,6);

  mat(1,1) = xr * xr;
  mat(1,2) = yr * yr;
  mat(1,3) = zr * zr;
  mat(1,4) = 2.0 * xr * yr;
  mat(1,5) = 2.0 * xr * zr;
  mat(1,6) = 2.0 * yr * zr;

  mat(2,1) = xs * xs;
  mat(2,2) = ys * ys;
  mat(2,3) = zs * zs;
  mat(2,4) = 2.0 * xs * ys;
  mat(2,5) = 2.0 * xs * zs;
  mat(2,6) = 2.0 * ys * zs;

  mat(3,1) = xt * xt;
  mat(3,2) = yt * yt;
  mat(3,3) = zt * zt;
  mat(3,4) = 2.0 * xt * yt;
  mat(3,5) = 2.0 * xt * zt;
  mat(3,6) = 2.0 * yt * zt;

  mat(4,1) = xs * xr;
  mat(4,2) = ys * yr;
  mat(4,3) = zs * zr;
  mat(4,4) = xs * yr + ys * xr;
  mat(4,5) = xs * zr + zs * xr;
  mat(4,6) = ys * zr + zs * yr;

  mat(5,1) = xt * xr;
  mat(5,2) = yt * yr;
  mat(5,3) = zt * zr;
  mat(5,4) = xt * yr + yt * xr;
  mat(5,5) = xt * zr + zt * xr;
  mat(5,6) = yt * zr + zt * yr;

  mat(6,1) = xt * xs;
  mat(6,2) = yt * ys;
  mat(6,3) = zt * zs;
  mat(6,4) = xt * ys + yt * xs;
  mat(6,5) = xt * zs + zt * xs;
  mat(6,6) = yt * zs + zt * ys;
end