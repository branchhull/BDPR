function [X1, X2] = BDPR_ADMM(MatA1, MatA2, delta, options)
% [X1, X2] = BDPR_ADMM(MatA1, MatA2, delta, options)
% solves the convex program:
%
%   min_{X_1,X_2} tr(X_1) + tr(X_2)
%   subject to:   <a_{1,l}a_{1,l}^*,X_1><a_{2,l}a_{2,l}^*,X_2> >= delta_l
%                                                           l = 1, ..., L
%                  X_1,X_2 >= 0
% In this function:
% - MatA1 is a complex matrix of size m x n1
% - MatA2 is a complex matrix of size m x n2
% - delta is a vector of length m of non-negative entries
%   ====================================================================
% - variable options contains the ADMM and optimization parameters as
% follows:
%
% *- options.ADMM.rho: the value of penalty parameter in the proposed ADMM
%    scheme (default = 1e-4)
%
% *- options.ADMM.maxIter: is the maximum number of ADMM iterations.
%    (default = 300)
%
% *- options.ADMM.plotFreq: after every few ADMM iterates the algorithm can
%    plot the first column (or row) of the recovered X1. This parameter
%    determines such frequency (default = inf)
%
% *- options.ADMM.ncvxMaxRankCoef: as discussed in the paper, to solve the
%    X-update we use a non-convex factorization X = VV', where the second
%    dimension of V is r. In order for this framework to find the global
%    minimizer we should have r > sqrt(c*m) where c > 2. This parameter
%    represents the factor c (default = 2.5)
%
% *- options.ADMM.convergeTol: Once the change in the dual variable reaches
%    this limit, the ADMM scheme is considered converged and the iterates
%    stop (default = 5e-4)
%
% *- options.ADMM.iterShowFreq: Every few iterations, program can report
%    the current number of ADMM iterations completed (default = 10)
%   ====================================================================
% ** options.minFunc.*: To perform our nonconvex solves we use the minFunc
%    toolbox developed by:
%
%    M. Schmidt, "minFunc: unconstrained differentiable multivariate
%    optimization in Matlab". minFunc, 2012. See the following webpage to
%    download the toolbox:
%
%    https://www.cs.ubc.ca/~schmidtm/Software/minFunc.html
%
%    All the options passed to the minFunc function are transferred to our
%    options.minFunc.* variable, where * denotes a minFunc original option
%
%    Our modified default values for the minFunc options are as follows:
%    options.minFunc.maxFunEvals = 100; (max calls to the non-convex cost)
%    options.minFunc.Method = 'lbfgs'; (optimization algorithm used)
%    options.minFunc.display = 'none'; (quiet mode)
%   ====================================================================
%    BDPR_ADMM is developed by A. Ahmed, A. Aghasi and P. Hand
%
%    The current code is a large-scale solver developed in Jan 2019. Feel
%    free to use/modify this code as long as you cite the corresponding
%    publications:

%    @incollection{NIPS2018_8207,
%    title = {Blind Deconvolutional Phase Retrieval via Convex Programming},
%    author = {Ahmed, Ali and Aghasi, Alireza and Hand, Paul},
%    booktitle = {Advances in Neural Information Processing Systems 31},
%    pages = {10051--10061},
%    year = {2018},
%    publisher = {Curran Associates, Inc.},
%    url = {http://papers.nips.cc/paper/8207-blind-deconvolutional-phase-
%           retrieval-via-convex-programming.pdf}
%    }

%    For questions about the code please contact aaghasi@gsu.edu.

%  =====================================================================

%% Extracting the values from the options variable and setting defaults
% setting the default values:
rho = 1e-4;
maxIter = 300;
plotFreq = inf;
ncvxMaxRankCoef = 2.5;
convergeTol = 5e-4;
iterShowFreq = 10;


if isfield(options,'ADMM')
    if isfield(options.ADMM,'rho')
        if(options.ADMM.rho > 0)
            rho = options.ADMM.rho;
        else
            error('value of rho should be strictly positive!')
        end
    end
    %
    if isfield(options.ADMM,'maxIter')
        maxIter = options.ADMM.maxIter;
    end
    %
    if isfield(options.ADMM,'plotFreq')
        plotFreq = options.ADMM.plotFreq;
    end
    %
    if isfield(options.ADMM,'ncvxMaxRankCoef')
        ncvxMaxRankCoef = options.ADMM.ncvxMaxRankCoef;
        if(options.ADMM.ncvxMaxRankCoef <= 2)
            warning('No guarantee for global convergence: value of ncvxMaxRankCoef should be strictly greater than 2!')
        end
    end
    %
    if isfield(options.ADMM,'convergeTol')
        convergeTol = options.ADMM.convergeTol;
    end
    %
    if isfield(options.ADMM,'iterShowFreq')
        iterShowFreq = options.ADMM.iterShowFreq;
    end
end

% minFunc defaults:
if isfield(options,'minFunc')
    minFuncOpt = options.minFunc;
    if ~isfield(minFuncOpt,'maxFunEvals')
        minFuncOpt.maxFunEvals = 100;
    end
    %
    if ~isfield(minFuncOpt,'Method')
        minFuncOpt.Method = 'lbfgs';
    end
    %
    if ~isfield(minFuncOpt,'display')
        minFuncOpt.display = 'none';
    end
else
    minFuncOpt = [];
    minFuncOpt.maxFunEvals = 100;
    minFuncOpt.Method = 'lbfgs';
    minFuncOpt.display = 'none';
end


%% Main ADMM Program

% making sure all elements of delta are positive
sdelta = delta >= 0;
if (prod(sdelta) == 0)
    error('Make sure all elements of delta are non-negative!');
end

m = length(delta);
% maxrank is the second dimension of V_j, such that X_j = V_jV_j'
maxrank = ceil(sqrt(ncvxMaxRankCoef*m));

n1 = size(MatA1,2);
n2 = size(MatA2,2);


% Initialization of the ADMM parameters
u1 = ones(m,1);
alpha1 = zeros(m,1);
V1 = randn(n1,maxrank);

u2 = ones(m,1);
alpha2 = zeros(m,1);
V2 = randn(n2,maxrank);




for i = 1 : maxIter
    ualpha1 = u1 + alpha1;
    ualpha2 = u2 + alpha2;
    [v1,~] = minFunc(@xupdate_cost,V1(:),minFuncOpt,MatA1,maxrank,ualpha1,rho);
    [v2,~] = minFunc(@xupdate_cost,V2(:),minFuncOpt,MatA2,maxrank,ualpha2,rho);
    V1 = reshape(v1,[n1,maxrank]);
    V2 = reshape(v2,[n2,maxrank]);
    
    aaX1 = sum(abs(MatA1*V1).^2,2);
    aaX2 = sum(abs(MatA2*V2).^2,2);
    
    
    
    [u1, u2] = projC(aaX1-alpha1, aaX2-alpha2, delta);
    
    if norm(u1 - aaX1)/m + norm(u2 - aaX2)/m < convergeTol
        break;
    end
    
    alpha1 = alpha1 + u1 - aaX1;
    alpha2 = alpha2 + u2 - aaX2;
    
    
    if (round(i/iterShowFreq)==(i/iterShowFreq))
        disp(strcat('ADMM iteration ----> ',num2str(i)));
    end
    
    
    if (round(i/plotFreq)==(i/plotFreq))
        close all;
        X1 = V1(1,:)*V1';
        plot(real(X1(1,:)));title('plot of one column of X1: recovered m^* so far')
        pause(.1);
    end
end

X1 = V1*V1';
X2 = V2*V2';
