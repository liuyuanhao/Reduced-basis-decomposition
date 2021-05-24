function [Y,T] = RBD(X ,eR ,dmax)
%RBD,Reduced Basis Decomposition,Yanlai Chen
%initialize set
d = 1;
Ecur = inf;
m = size(X,1);
n = size(X,2);
i = unidrnd(n);

Y = zeros(m,d);
T = zeros(d,n);

% Z = zeros(m,1);
%2.1&2.2
%Apply the modified Gram-Schmidt orthonormalization 
%to obtain the dth basis of the compressed space

while d <= dmax && Ecur > eR
%     if Z(i) == 0
%         Z(i) = 1;
%     end
    v = X(:,i);
    for j = 1:d-1
        v = v - dot(v,Y(:,j)) * Y(:,j);
    end
    
    if norm(v) < eR
        Y = Y(:,1:d-1);
        T = T(1:d-1,:);
        break;
    else
        Cd = v/norm(v);
        Y(: ,d) = Cd;
        T(d ,:) = Cd' * X;
    end
    
%2.3 
    Earray = zeros(1,n);
    for j = 1:n
        Earray(j) = norm(X(:, j) - Y(:, 1:d)*T(:, j));
    end
    Ecur = max(Earray);
    i = find(Earray == Ecur);
    i = i(1);%argmax 
    
%2.4
    if Ecur <= eR
        Y = Y(:,1:d);
        T = T(1:d ,:);
    else
        d = d + 1;
    end
end 