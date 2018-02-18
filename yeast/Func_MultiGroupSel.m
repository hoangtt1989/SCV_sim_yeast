function [Best, numITs] = Func_MultiGroupSel(X, Y, Lambda, maxIT, thresholdingWay, BInit, Nu, K, vecBased)
% Perform Group Selection for Multivariate Linear Regression
% K: the scaling constant in the TISP iterations to guaranteed convergence.
%   Note that when the convergence is guaranteed, the smaller the value of K
%   is, the faster the convergence is. 
% Nu: the ridge shrinkage parameter for the scaled model. That is, the
%   problem is defined as  |Y - X B|_F^2/(2 K) + Nu^2 |B|_F^2/2 + sparsity pen or sparsity constraint. 
% BInit: The initial quantity. If not supplied, taken to be zero. 
% vecBased: vectorize the matrix and perform thresholding (useful for
%   sparse RRR and sparse PCA (in comparison to selectable RRR and SEL-PCA)


if ~exist('vecBased', 'var') || isempty(vecBased)
    vecBased = 0; % default: matrix row norm based thresholding
end

if ~exist('K', 'var') || isempty(K)
    K = 1;
end
if ~exist('Nu', 'var') || isempty(Nu)
    Nu = 0;
end
if maxIT == -1
    maxIT = 1e+4;
end
errBnd = 1e-4;
% Compute the whole solution path
[m, n] = size(Y); d = size(X, 2);

% We assume X has been scaled to have 2-norm 1
% oriX = X;
% k0 = norm(X, 2);
% X = X / k0;

% Bs = zeros(d, n, numel(lambdas));
% for lambdaInd = 1:numel(lambdas)
%     Lambda = lambdas(lambdaInd);
        if ~exist('BInit', 'var') || isempty(BInit)
        % ZERO start
            curB = zeros(d, n); %rand(d, m); %
        else
            curB = BInit;
        end
        nIter = 1;
        while 1
            Xi = curB + X' * (Y - X * curB) / K;
            if vecBased
                newB = MultiTh_mat(Xi(:), Lambda, thresholdingWay, Nu);
                newB = reshape(newB, size(Xi));
            else
                newB = MultiTh_mat(Xi, Lambda, thresholdingWay, Nu);
            end
            if max(max(abs(curB - newB))) < errBnd | nIter >= maxIT
                break;
            else
                curB = newB;
                nIter = nIter + 1;        
            end
        end
        Best = newB;
        numITs = nIter;
% end
% Bs = Bs / k0;
% X = oriX;

