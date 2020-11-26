function [x,v] = StdSBLSolver(y,A)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
n = size(A,2);
m = size(y,1);

a = 1e-10;
b = 1e-10;
c = 1e-10;
d = 1e-10;

iter = 0;
iter_mx = 100;
D = eye(n);
sigma2 = 1;
alpha_new = ones(n,1);
var_new = inv(A'*A/sigma2+D);
mu_old = ones(n,1);
mu_new = 1/sigma2*var_new*A'*y;
gamma_new = 1/sigma2;
while iter<iter_mx && norm(mu_new-mu_old)>1e-6
    iter = iter + 1;
    mu_old = mu_new;
    omega = diag(var_new) + mu_new .* conj(mu_new);
    alpha_new = (1+2*a) ./ (omega+2*b);
    idx1 = find(alpha_new>1e10);
    alpha_new(idx1) = 1e10;
    
    D = diag(alpha_new);
    %=============================================
    %  estimate the variance of noise
    num=(y-A*mu_old)'*(y-A*mu_old)+trace(var_new*A'*A)+2*d;
    den=m+2*c;
    sigma2=num/den;
    %==============================================
    var_new=inv(A'*A/sigma2+D);
    mu_new=1/sigma2*var_new*A'*y;
end
x = mu_new;
v = diag(var_new);
end

