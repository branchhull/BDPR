function [X1, X2] = PSD_BH(MatA1, MatA2, delta, max_iter, iter_frequency)
% [X1 X2] = PSD_BH(A1, A2, delta, rho1, rho2, max_iter, iter_frequency)
% solves the convex program:
%
%   min_{X_1,X_2} tr(X_1) + tr(X_2)
%   subject to:   <a_{1,l}a_{1,l}^*,X_1><a_{2,l}a_{2,l}^*,X_2> >= delta_l 
%                                                           l = 1, ..., L
%                  X_1,X_2 >= 0
% In this function:
% - MatA1 is a complex matrix of size LxN
% - MatA2 is a complex matrix of size LxK
% - delta is a vector of length L of non-negative entries
% - maxiter is the maximum ADMM number of iterations
% - iter_frequency is the number of iterations that are performed before
%   each report of the program status


% making sure all elements of delta are positive
sdelta = delta >= 0;
if (prod(sdelta) == 0)
    error('Make sure all elements of delta are non-negative!');
end

L = length(delta);
N = size(MatA1,2);
K = size(MatA2,2);

MatA1 = MatA1';
MatA2 = MatA2';

AA1 = zeros(N^2,L);
AA2 = zeros(K^2,L);

for i = 1 : L
    temp = MatA1(:,i)*MatA1(:,i)';
    AA1(:,i) = temp(:);
    temp = MatA2(:,i)*MatA2(:,i)';
    AA2(:,i) = temp(:);
end

rho2 = 1;

A1 = AA1*AA1';
A2 = AA2*AA2';

rho1 = .5/(mean(real(diag(A1))) + mean(real(diag(A2))) );

% calculating A1 and A2 as specified in the paper
A1 = rho1*A1 + rho2*eye(N^2);
A2 = rho1*A2 + rho2*eye(K^2);

C1 = chol(A1);
C2 = chol(A2);

% making the Cholesky factors sparse
C1 = sparse(C1);
C2 = sparse(C2);

% Initialization of the parameters
u1 = zeros(L,1);
alpha1 = zeros(L,1);
P1 = zeros(N,N);
Z1 = zeros(N,N);

u2 = zeros(L,1);
alpha2 = zeros(L,1);
P2 = zeros(K,K);
Z2 = zeros(K,K);

for i = 1:max_iter
    ualpha1 = u1 + alpha1;
    ualpha2 = u2 + alpha2;
    [X1, X2] = solveTraceLS(C1, C2, AA1, AA2, Z1, Z2, P1, P2, ualpha1, ualpha2, rho1, rho2);
    % to make things Hermitian
    X1 = (X1 + X1')/2;
    X2 = (X2 + X2')/2;
    
    Z1 = projPSD(X1 + P1);
    Z2 = projPSD(X2 + P2);
    % to make things Hermitian
    Z1 = (Z1 + Z1')/2;
    Z2 = (Z2 + Z2')/2;
    
    aaX1 = real(AA1'*X1(:));
    aaX2 = real(AA2'*X2(:));
    [u1, u2] = projC(aaX1-alpha1, aaX2-alpha2, delta);
    alpha1 = alpha1 + u1 - aaX1;
    alpha2 = alpha2 + u2 - aaX2;
    P1 = P1 + X1 - Z1;
    P2 = P2 + X2 - Z2;
    if (round(i/iter_frequency)==(i/iter_frequency))
       disp(strcat('ADMM iteration ----> ',num2str(i))); 
%        close all;pause(.01);
%        plot(real(X1(1,:)));title('recovered signal');
    end
end


