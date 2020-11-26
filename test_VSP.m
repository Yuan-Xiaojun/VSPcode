%% Initialization
clear;
n=200;                                          % signal dimension
m=75;                                           % number of measurements
K=30;                                           % total number of nonzero coefficients
L=1;                                            % number of nonzero blocks
% rng(42)

for simi = 1:500
    simi
    for snri = -1:6
        SNR = snri * 5;
        % generate the block-sparse signal
        x=zeros(n,1);
        r=abs(randn(L,1));
        r=r+1;
        r=round(r*K/sum(r));
        r(L)=K-sum(r(1:L-1)); % number of non-zero coefficients in each block
        g=round(r*n/K);
        g(L)=n-sum(g(1:L-1));
        g_cum=cumsum(g);
        
        for i=1:L
            % generate i-th block
            seg=1*(randn(r(i),1)+1j*randn(r(i),1)); % generate the non-zero block
            loc=randperm(g(i)-r(i)); % the starting position of non-zero block
            x_tmp=zeros(g(i), 1);
            x_tmp(loc(1):loc(1)-1+r(i)) = seg;
            x(g_cum(i)-g(i)+1:g_cum(i), 1)=x_tmp;
        end
        
        % generate the measurement matrix
        A = (randn(m,n) + 1j*randn(m,n))/sqrt(2);
        % noiseless measurements
        measure=A*x;
        sigma2 = norm(measure(:))^2/m*10^(-SNR/10);
        % Observation noise, stdnoise = std(measure)*10^(-SNR/20);
        stdnoise=sqrt(sigma2);
        noise=(randn(m,1)+1j*randn(m,1))/sqrt(2)*stdnoise;
        
        % Noisy measurements
        y=measure+noise;
        
        %% Standard SBL Solver
%         [x_new, ~] = StdSBLSolver(y,A);
        %% VSP Solver
        [x_new, ~] = VSPSolver(y,A,2*K,3,1e-4);
        %% OMP
%         [~,x_new] = SOMP2(y,A,K);
        %% PCSBL Solver
%         [x_new,~] = PCSBLSolver(y,A);
        %% LMMSE Lower Bound
%         nzindex = find(x~=0);
%         A_ = A(:,nzindex);
%         x_ = x(nzindex);
%         x_new = inv(1/var(x_)*eye(K)+1/sigma2*A_'*A_) * A_'*1/sigma2*y;
%         mse=norm(x_new-x_)^2;
%         nmse=mse/(norm(x_))^2;
%         MSE(snri+2,simi) = mse;
%         NMSE_sim(snri+2,simi) = nmse;
        
        mse=norm(x_new-x)^2;
        nmse=mse/(norm(x))^2;
        NMSE_sim(snri+2,simi) = nmse;
    end
end
NMSE = mean(NMSE_sim,2);

plot(-5:5:30,10*log10(NMSE));
xlabel('SNR (dB)');
ylabel('NMSE (dB)');