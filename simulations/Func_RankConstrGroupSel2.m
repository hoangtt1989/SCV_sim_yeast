function [Best, Shat, Vhat, nIter, remainingDims, rankEst] = Func_RankConstrGroupSel2(X, Y, Lambda, k, errBnd, maxIT, grpSelMaxIT, thresholdingWay, Nu, BInit, K, vecBased, simpFlag)
% A rewritten version of Func_RankConstrGroupSel
% Main changes: 
%   Order of two alterantive optimization steps.
%   Initalizations.
%   Add scaling constant K.

% vecBased: vectorize the matrix and perform thresholding (useful for
%   sparse RRR and sparse PCA (in comparison to selectable RRR and SEL-PCA)

output = 0;
if ~exist('simpFlag', 'var') || isempty(simpFlag)
    simpFlag = 0; % When =1, the identity is I, and thus some computations can be simplified and performed in a more efficient manner
end

if ~exist('vecBased', 'var') || isempty(vecBased)
    vecBased = 0; % default: matrix row norm based thresholding
end

if ~exist('K', 'var') || isempty(K)
    K = 1;
end
if ~exist('Nu', 'var') || isempty(Nu)
    Nu = 0;
end
% Initialization
if ~exist('BInit', 'var') || isempty(BInit) || (isempty(BInit.BInit) & isempty(BInit.SInit))
%     B_cur = zeros(size(X,2), size(Y,2)); %  the zero start can be used
    if ~simpFlag
        [B_cur, tmpr, eigvals_rrr, Shat, Vhat] = Func_RRR(X, Y, k); %Use the RRR estimate
        Vhat = Vhat';
    else
        [B_cur, Vhat, Shat] = Func_PCA(Y, k);
    end
else % Some initializations are from the arguments
    if isempty(BInit.BInit) % In this case BInit.SInit must be nonempty
        % Init based on BInit.SInit. (Compute B_cur)
        Shat = BInit.SInit;
        if ~isempty(BInit.VInit)
            Vhat = BInit.VInit;
            B_cur = Shat * Vhat';
        else
            Vhat = [];
            B_cur = zeros(size(X,2), size(Y,2));
        end
        
    else % BInit.BInit nonempty
        % Init based on BInit.BInit. If necessary, compute Shat (and Vhat).
        B_cur = BInit.BInit;
        if isempty(BInit.SInit)
            [tU,tS,tV] = svd(B_cur, 0); %  NOT svd(X * B_cur, 0) because B=SV'
            Vhat = tV(:, 1:k); 
            Shat = B_cur * Vhat;
        else 
            Shat = BInit.SInit;
            if ~isempty(BInit.VInit)
                Vhat = BInit.VInit;
            else
                Vhat = [];
            end
        end
        
    end
end
% K = norm(X, 2)^2;
% if ~exist('squeezing', 'var') || isempty(squeezing)
%     squeezing = 0;
% end
% remainingDims = [1:size(X,2)];

nIter = 1;
while 1
    % 1. Get the new V estimate
    if nIter > 1 || isempty(Vhat) % The other case: Vhat is provided
        [U_tmp, S_tmp, V_tmp] = svd(Y' * X * Shat, 'econ');
        Vhat = U_tmp * V_tmp';
    end
    % 2. Get the new S estimate
    %%% Construct new responses
    Y_pseudo = Y * Vhat;
    %%% run group selection
    [Shat, numITs] = Func_MultiGroupSel(X, Y_pseudo, Lambda, grpSelMaxIT, thresholdingWay, Shat, Nu, K, vecBased);
    
	if output, disp(['Lambda=', num2str(Lambda), ', k=', num2str(k), '; numITs ', num2str(numITs)]), end
    
%     % Squeeze the design (and maybe the response) for the purpose of screening/computational efficiency
%     if squeezing == 1
%         tmpRemainingDims = sqrt(sum(Shat'.^2)) > 1e-4;
%         remainingDims = remainingDims(tmpRemainingDims);
%         Shat = Shat(tmpRemainingDims,:);
%         X = X(:,tmpRemainingDims);
%     elseif squeezing == 2
%         tmpRemainingDims = sqrt(sum(Shat'.^2)) > 1e-4;
%         remainingDims = remainingDims(tmpRemainingDims);
%         Shat = Shat(tmpRemainingDims,:);
%         X = X(:,tmpRemainingDims);
%         Y = Y(:,tmpRemainingDims);
%     end
    % 3. Form the new B estimate
    B_new = Shat * Vhat';
    
    if max(max(abs(B_new - B_cur))) < errBnd | nIter >= maxIT
        break;
    else
        B_cur = B_new;
        nIter = nIter + 1;        
    end
end
Best = B_new;
remainingDims = find(sqrt(sum(Best'.^2)) > 1e-4);
rankEst = sum(svd(Best)>1e-4);

if output, disp(['Lambda=', num2str(Lambda), ', k=', num2str(k), '; outer ITs ', num2str(nIter)]), end

