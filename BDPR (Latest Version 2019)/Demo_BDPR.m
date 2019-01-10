m = 2000; % signal canonical dimension
N = 40; % C subspace dimension
K = 40; % B subspace dimension

rng(1)
B = randn(m,K);
B = B/sqrt(m);


C = randn(m,N);
C = C/sqrt(m);

% defining the true signal subspace coefficients m*: vector of +1 and -1
%                                                h*: vector of all 1s
ms = [ones(round(N/2),1);-ones(N-round(N/2),1)];
hs = ones(K,1);

MatA1 = fft(C);
MatA2 = fft(B);

% calculating the product of bilinear signal magnitudes
yhat = abs(MatA1*ms).*abs(MatA2*hs);

% setting the BDPR_ADMM options for this demo. Please refer to the function
% BDPR_ADMM for a list of all available options and their description

options = [];
options.ADMM.maxIter = 100;
options.ADMM.plotFreq = 10;
options.ADMM.iterShowFreq = 5;
options.ADMM.convergeTol = 2e-4;
% =========================================================================
% If you want the program to run faster you can set options.ADMM.ncvxMaxRankCoef
% to values below 2 (e.g., uncomment the line below), however, there is no
% theoretical guarantee for convergence
% =========================================================================
% options.ADMM.ncvxMaxRankCoef = .1;

options.minFunc.maxFunEvals = 50;

tic;
[X1, X2] = BDPR_ADMM(MatA1, MatA2, yhat.^2, options);
toc;

