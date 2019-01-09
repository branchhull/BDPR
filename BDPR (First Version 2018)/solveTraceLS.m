function [X1, X2] = solveTraceLS(C1, C2, AA1, AA2, Z1, Z2, P1, P2, ualpha1, ualpha2, rho1, rho2)
% function [x1 x2] = solveTraceLS(A1, A2, AA1, AA2, Z1, Z2, P1, P2, ualpha1, ualpha2, rho1, rho2)
% solves the X_j update step of the proposed ADMM scheme. In this function:
% - C1 and C2: are Cholesky matrices of size N^2xN^2 and K^2xK^2 calculated
%   through the Cholesky factorization of A_j formulation in the paper
% - AA1 and AA2 are matrices of sizes N^2xL and K^2xL where each column is 
%   vec(a_i*a_i^*)
% - Z1 and Z2, P1 and P2 are Hermitian matrices of sizes NxN and KxK
% - ualpha1 and ualpha2 are vectors of size L that are the sum of u and
%   alpha vectors
% - rho1 and rho2 are the ADMM penalty parameters 

N = size(Z1,1);
K = size(Z2,1);

rhs1 = rho1*AA1*ualpha1;
temp1 = rho2*(Z1-P1)-eye(N);
rhs1 = rhs1 + temp1(:);

rhs2 = rho1*AA2*ualpha2;
temp2 = rho2*(Z2-P2)-eye(K);
rhs2 = rhs2 + temp2(:);

x1 = C1\(C1'\rhs1);
x2 = C2\(C2'\rhs2);

X1 = reshape(x1, [N,N]);
X2 = reshape(x2, [K,K]);


