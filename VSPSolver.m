function [x,v] = VSPSolver(y,A,K,beta,alpha)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
n = size(A,2);
m = size(y,1);
rho = K/n;

a = 1e-10;
b = 1e-10;
c = 1e-10;
d = 1e-10;

iter = 0;
iter_mx = 80;
alpha_new = ones(n,1);
D = diag(alpha_new);
sigma2 = 1;
var_new = inv(A'*A/sigma2+D);
mu_new = 1/sigma2*var_new*A'*y;
% gamma_new = 1/sigma2;
while iter<iter_mx %&& norm(mu_new-mu_old)>1e-6
    iter = iter + 1;
    omega = diag(var_new) + mu_new .* conj(mu_new);
    alpha_new = (1+2*a) ./ (omega+2*b);
    idx1 = find(alpha_new>1e10);
    alpha_new(idx1) = 1e10;
    
    if mod(iter,30)==0
        
        v_old = 1./alpha_new;
        
        max_in = mean(maxk(v_old,K));
        factor_in = 1/ max_in ;
        v_old = factor_in * v_old;
        
        v_old = min(v_old,1);
        %v_old(find(v_old>1e-2)) = 1;
        %                 [ v_new, ~, ~ ] = SPD_MC(v_old, K/n, L/K );
        [v_new,~] = SPD_MRF1D(v_old.', beta, alpha, 100);
        v_new = v_new.';
               
%         v_new = v_new / rho;
        
%         max_out = mean(maxk(v_new,K));
%         factor_out = 1 / rho / max_out;
        
        v_new = v_new / factor_in;
        
        alpha_new = 1./v_new;
    end
    
    D = diag(alpha_new);
    %=============================================
    %  estimate the variance of noise
    num=(y-A*mu_new)'*(y-A*mu_new)+trace(var_new*A'*A)+2*d;
    den=m+2*c;
    sigma2=num/den;
    %==============================================
    var_new=inv(A'*A/sigma2+D);
    mu_new=1/sigma2*var_new*A'*y;
end
x = mu_new;
v = diag(var_new);
end

