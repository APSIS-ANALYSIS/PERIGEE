% DerOneBasisFun is a function computing the kth derivative of N_i,p for k
% =0,...,n  n<=p.

% p is the degree 
% n is the highest degree derivatives we compute
% i is the number of the basis function, ranging [0,1,...,m-p-1]
% U is the knot vector
% u is the point we evaluate

% output is ders[] where ders[k] = kth derivative of N_i,p, where k = 0,1,...,n.
% m+1 is the length of the knot vector
% This is the algorithm A2.5 in NURBS Book

% By Ju Liu, 2009 Dec.

function ders = DerOneBasisFun(p , n , i , U , u)
N = sparse(p+1,p+1);

% Local property
ders = zeros(n+1,1);
if (u<U(i+1) || u>=U(i+p+2))
    return;
end

% Initialize zeroth degree functions
for j = 0 : p
    if (u>=U(i+j+1) && u < U(i+j+2))
        N(j+1,1) = 1.0;
    end
end

% Compute full triangular table
for k = 1 : p
    if N(1,k) == 0.0
        saved = 0.0;
    else 
        saved = ((u-U(i+1))* N(1,k)) / (U(i+k+1) - U(i+1));
    end
    for j= 0 : p-k
        Uleft = U(i+j+2);
        Uright = U(i+j+k+2);
        if (N(j+2,k) == 0.0)
            N(j+1,k+1) = saved;
            saved = 0.0;
        else
            temp = N(j+2,k)/(Uright - Uleft);
            N(j+1,k+1) = saved + (Uright-u) * temp;
            saved = (u-Uleft) * temp;
        end
    end
end
% matrix N stores the entire triangular table corresponding to k = 0.

ders(1) = N(1,p+1);

% compute the derivatives up to degree n.
for k = 1 : n
    ND = zeros(k+1,1);
    for j = 0 : k
        ND(j+1) = N(j+1,p-k+1);
    end
    for jj = 1 : k
        if ND(1) == 0.0
            saved = 0.0;
        else
            saved = ND(1)/(U(i+p-k+jj+1)-U(i+1));
        end
        for j = 0 : k-jj
            Uleft = U(i+j+2);
            Uright = U(i+j+p-k+jj+2);
            if (ND(j+2) == 0.0)
                ND(j+1) = (p-k+jj)*saved;
                saved = 0.0;
            else
                temp = ND(j+2)/(Uright - Uleft);
                ND(j+1) = (p-k+jj) * (saved -temp);
                saved = temp;
            end
        end
    end
    ders(k+1) = ND(1);
end
  