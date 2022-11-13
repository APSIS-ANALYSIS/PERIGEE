% OneBasisFun is a function computing the value of a single basis function
% N_i,p(u).

% p is the degree
% i is the number of the basis function, range [0,1,...,m-p-1]
% u is the point we want to evaluate
% U is the knot vector

% Output is the value of the function N_i,p(u)

% This is the algorithm A2.4 in NURBS Book

% By Ju Liu, 2009 Dec.

function Nip = OneBasisFun(u,i,p,U)

m = length(U)-1;        % m+1 is the length of the knot vector; 

n = m-p-1;              % N_1 ... N_n+1, the dimensional of basis function of degree p is n+1;


if ((i==0 && u == U(1)) || (i == n && u == U(m+1)))
    Nip = 1.0;
    return;
end

if (u<U(i+1) || u >= U(i+p+2))
    Nip = 0.0;
    return;
end

% Initialize zeroth degree functions
N = zeros(p+1,1);
for j = 0 : p
    if (u >= U(i+j+1) && u < U(i+j+2))
        N(j+1) = 1.0;
    end
end

% Compute according to the trangular table
for k = 1 : p
    if N(1) == 0.0 
        saved = 0.0;
    else
        saved = ((u-U(i+1))*N(1))/(U(i+k+1)-U(i+1));
    end
    for j = 0 : p-k
        Uleft = U(i+j+2);
        Uright = U(i+j+k+2);
        if (N(j+2) == 0.0)
            N(j+1) = saved;
            saved = 0.0;
        else 
            temp = N(j+2)/(Uright-Uleft);
            N(j+1) = saved + (Uright - u)*temp;
            saved = (u - Uleft) * temp;
        end
    end
end
Nip = N(1);

    