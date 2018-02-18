function [Bs_corrt] = Func_EstsCorrt(Bs, X, Y, tuningparams, nzPatts, rankBs, numBs)
% Given the estimate, compute the RRR estimate with the same sparsity pattern and rank.
Bs_corrt = zeros(size(Bs));
for paramInd = 1:numBs % 1:numel(tuningparams)
    tmpB = zeros(size(Bs(:,:,paramInd)));
    [tmpB2] = Func_RRR(X(:, nzPatts(paramInd, :)==1), Y, rankBs(paramInd));
    tmpB(nzPatts(paramInd, :)==1, :) = tmpB2;
    Bs_corrt(:, :, paramInd) = tmpB;
end

