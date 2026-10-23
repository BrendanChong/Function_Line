clear; clc; close all;
syms xs ys xi eta;

p = phi(xs,ys,0);
gradp = gradient(p);

N = basis(xi,eta);
dN = 2/3.*jacobian(N,[xi eta]);
xi = -0.27415567780803768;
eta = 0.52359877559829893;
xs = -0.411234; ys = 0.785398;

xHat = [-1.5 -1.5; 1.5 -1.5; 1.5 1.5; -1.5 1.5];
for i = 1:size(xHat,1)
    uHat(i) = phi(xHat(i,1), xHat(i,2), 0);
end

gradBasis = double(subs([sum(dN(:,1).*uHat'); sum(dN(:,2).*uHat')]))

function N = basis(xi,eta)
    N(1) = 0.25*(1-xi)*(1-eta);
    N(2) = 0.25*(1+xi)*(1-eta);
    N(3) = 0.25*(1+xi)*(1+eta);
    N(4) = 0.25*(1-xi)*(1+eta);
end