function [fval,g]=xupdate_cost(v,A,maxrank,xi,rho)
% This function computes the value and gradient of the X-Update objective 
% in the proposed ADMM scheme. The objective is
%   obj = ||V||_F^2 + rho/2*sum_l (||V'a_l||^2 - xi_l)^2
% where V is the reshaped form of v into a matrix with maxrank columns
% IMPORTANT: this code assumes v is a real signal, although A is a complex
% matrix

m = size(A,1); % m is the number of measurements
n = size(A,2); % n is the size of the signal

V = reshape(v,n,maxrank);
if length(xi) ~= m 
    error('size mismatch: length of xi must be equal to the number of rows in A!');
end

AV = A*V;
summand = (sum(abs(AV).^2,2)) - xi;
fval = norm(V,'fro')^2 + rho/2*sum(summand.^2);

g = reshape(2*V + 2*rho*real((A'*spdiags(summand,0,m,m))*AV),[],1);