clear all; clc; clf;
hold on;

x = 0: 0.01 : 1;


p = 3;
%U = [-1, -2/3, -1/3, 0, 1/3, 2/3, 1, 1+1/3,1+2/3,2];

U = [0,0,0,0,1,2,3,4,5,6,6,6,6]/6;

m = length(U) - 1;
y = zeros(m-p-1,101);

for j = 0 : m-p-1
  for i = 1 : 101
    y(j+1,i) = OneBasisFun(x(i),j, p , U);
  end
end

plot(x,y);