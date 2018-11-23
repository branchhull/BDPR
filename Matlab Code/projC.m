function [u1, u2] = projC(z1, z2, del)
% function [u1 u2] = projC(z1, z2, del)
% projects a vector in the form of (z1,z2) in R^{2L} onto the set C
% defined in the paper. In this function all inputs should have identical
% sizes of L-by-1

L = size(z1,1);
if L == 1
    error('Make sure all the inputs are vertical vectors of the same size!');
end

% making sure all inputs are of the same size
AllSizes = [size(z1,1), size(z2,1), size(del,1)];
if (norm(diff(AllSizes))>0)
    error('Make sure all the inputs are vertical vectors of the same size!');
end

% making sure all elements of delta are positive

if (~all(del>=0))
    error('Make sure all the elents of delta are non-negative!');
end

% Our assumptions are based on having nonzero entries for y, so if an entry
% of y is zero we replace it with a very small number
del(del==0) = 1e-12;

% zero tolereance
mac = 1e-13;


% checking if a point meets any of the inequality constraints
C1 = del - z1.*z2;
C2 = z1;

inSet = (C1 <= mac)&(C2>=-mac);

% setting the default projection to the point itself
u1 = z1;
u2 = z2;

% forming a matrix of L-by-5 as the coefficient matrix for all L
% polynomials

pCoef = [ones(L,1), -z1, zeros(L,1), del.*z2, -del.^2];

% lopping to find the projections one after the other

for i = 51 : L
    % skipping the points with trivial projection (projection of the points
    % that are already in the set is trivial)
    if (inSet(i)>0)
        continue;
    end
    pRoots = roots(pCoef(i,:));
    pRoots = pRoots(abs(imag(pRoots))<1e-15);
    
    isFeasible = (pRoots>=-mac)&( (del(i) - z2(i)*pRoots) >=-mac );
    % if none of the roots meet the required conditions we let the
    % projection be the point itself
    if sum(isFeasible) == 0
        continue;
    end
    % if for some numerical error multiple roots meet the condition we pick
    % one of them
    if sum(isFeasible) > 1
        I = find(isFeasible);
        isFeasible(I(2:end)) = 0;
    end
    pRoots = pRoots(isFeasible);
    if abs(pRoots)<1e-12
        u1(i) = max(u1(i),0);
        u2(i) = max(u2(i),0);
    else
        
        u1(i) = pRoots;
        mu = (del(i) -z2(i)*pRoots)/(pRoots^2);
        u2(i) = z2(i) + mu*pRoots;
    end
end




