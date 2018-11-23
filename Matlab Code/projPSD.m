function Xp = projPSD(X)
% function Xp = projPSD(X): takes a hermitian matrix X and returns its
% projection onto the positive semidefinite cone

% making sure the matrix is Hermitian:
if ~(ishermitian(X))
    error('Make sure the input matrix is Hermitian!');
end

n = size(X,2);

[U, D] = eig(X);
dp = max(diag(D),0);

% Dp = spdiags(dp,0,sparse(n,n));
% Xp = (U*Dp)*U';

r = min(find(dp>0));
if sum(dp>0) == 0
    Xp = zeros(size(X));
    return;
end

dpr = dp(r:n);
nr = n-r+1;
Dpr = spdiags(dpr,0,sparse(nr,nr));
Ur = U(:,r:n);
Xp = (Ur*Dpr)*Ur';

