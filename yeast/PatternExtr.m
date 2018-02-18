%%%%%%%%%%%%%%%%%%%% Sparsity and/or Low r Pattern Extraction based on An Estimator %%%%%%%%%%%%%%%%%%%%%
function [Z, newDim, TransformMat, J, r, sparPatt, Cest] = PatternExtr(X, B, pcamanner, sparPatt, rankB) 
% Arguements:   X -- model matrix, 
%               B -- possibly sparse and/or low r coefficient estimate,
% Values:       Z -- new design matrix with only newDim features as
%                   linear combinations of the orignal features: Z=X*TransformMat. 
%               newDim -- (==J or r) the number of new features, should be equal to
%               the row Jort size of B if the squeezed B has full r. Otherwise it is the r of B
%               TransformMat -- the accumulated (orthogonal) transformation matrix to X.
%               J -- the number of nonzero rows in B.
%               r -- the r of B.

if ~exist('pcamanner', 'var') || isempty(pcamanner)
    pcamanner = 0;
end
if ~exist('sparPatt', 'var') || isempty(sparPatt)
    sparPatt = find(sqrt(sum(B'.^2)) > 1e-4);
else
    sparPatt = find(sparPatt);
end

S = eye(size(X,2));
S = S(:,sparPatt); 

B = B(sparPatt,:); 

if pcamanner == 0
    % Perform Type-I feature extraction, i.e., extract the range space of
    % the estimate to make the new design. 
    [U, singvals, V] = svd(B, 'econ');
    singvals = diag(singvals);
    singvals(singvals<1e-4) = 0;
    inds = find(singvals>0);
    if ~exist('rankB', 'var') || isempty(rankB)
        rankB = numel(inds);
    else
        if rankB ~= numel(inds), error('wrong rank parameter'), end;
    end
    if rankB < min(size(B)) 
        U = U(:, inds);
        TransformMat = S * U;
    else% In this case, no need to put an extra orthornomral matrix, to be able to apply to variable selection
        TransformMat = S;
    end
    Cest = [];
else
    % PCA loading directions extraction: 
    tmpX = X(:, sparPatt);
    [tU, tS, tV] = svd(tmpX, 'econ');
    [tS, tIX] = sort(diag(tS), 'descend');
    tV = tV(:, tIX);
    if rankB > size(tV,2)
        warning('wrong rank parameter')
    end
    rankB = min(rankB, size(tV,2));
    % In whatever case, extract the PC directions to make a new design
%     if rankB < min(size(B))
        U = tV(:, 1:rankB); % 
        Cest = [tV(:, 1:rankB); zeros(size(X,2)-numel(sparPatt), rankB)]';
        TransformMat = S * U;
%     else
%         TransformMat = S;
%     end
    
end

Z = X * TransformMat;

newDim = size(Z, 2);
J = numel(sparPatt);
r = size(TransformMat,2);
