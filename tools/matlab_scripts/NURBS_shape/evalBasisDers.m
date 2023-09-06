function ders = evalBasisDers(i,s,p,n,S)
 
% Computes nonzero basis functions of order p and their derivatives up to the nth
% derivative at point s located in knot span i of knot vector S
% Directly from the NURBS Book by Piegl and Tiller

ndu(1,1) = 1.0;
left = s - S( i + 1 - ( 1 : p ) );
right = S( i + ( 1 : p ) ) - s;
for j = 1:p
    %left(j) = s - S(i + 1 - j);
    %right(j) = S(i + j) - s;
    saved = 0.0;
    for r = 0 : j - 1
        ndu(j + 1,r + 1) = right(r + 1) + left(j - r);
        temp = ndu(r + 1,j) / ndu(j + 1,r + 1);
        ndu(r + 1,j + 1) = saved + right(r + 1) * temp;
        saved = left(j - r) * temp;
    end
    ndu(j + 1,j + 1) = saved;
end
ders(1,1 : p + 1) = ndu(1 : p + 1,p + 1);

for r = 0:p
    s1 = 0;
    s2 = 1; 
    a(1,1) = 1.0;
    for k = 1:n
        d = 0.0;
        rk = r - k;
        pk = p - k;
        if(r >= k)
            a(s2 + 1,1) = a(s1 + 1,1) / ndu(pk + 2,rk + 1);
            d = a(s2 + 1,1) * ndu(rk + 1,pk + 1);
        end
        if(rk >= -1)
            j1 = 1;
        else 
            j1 = -rk;
        end
        if(r - 1 <= pk)
            j2 = k - 1;
        else 
            j2 = p - r;
        end
        for j = j1:j2
            a(s2 + 1,j + 1) = (a(s1 + 1,j + 1) - a(s1 + 1,j)) / ndu(pk + 2,rk + j + 1);
            d = d + a(s2 + 1,j + 1) * ndu(rk + j + 1,pk + 1);
        end
        if(r <= pk)
            a(s2 + 1,k + 1) = -a(s1 + 1,k) / ndu(pk + 2,r + 1);
            d = d + a(s2 + 1,k + 1) * ndu(r + 1,pk + 1);
        end
        ders(k + 1,r + 1) = d;
        j = s1;
        s1 = s2;      
        s2 = j;
    end
end

r = p;
for k=1:n
     for j=0:p
         ders(k+1,j+1) = ders(k+1,j+1)*r;
     end
    r = r*(p-k);
end