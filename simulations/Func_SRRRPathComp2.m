function [Bs, rankBs, grpBs, gridLs, nzPatts, trnerrs, numBs, Ss, Vs, cand_ranks_grid] = Func_SRRRPathComp2(X, Y, probType, regType, ...
    rankGiven, rankGridUBnd, lambdasGiven, full2d, nPoints, maxIT, grpSelMaxIT, errBnd, thresholdingWay, Nu, rem_ctrl, output)

% A function version based on SRRRPathComp2, the revised version of SRRRPathComp.
% Main changes:
%   Call Func_RankConstrGroupSel2;
%   Make a grid for |Y - X B|_F^2/k0^2 + lambda P(B) by scaling Y as well.
%   lambdasGiven: The selection parameter grid (lambdas) can be given by lambdasGiven.
%       Note that this is for the original dataset (X, Y) instead of (X/k0, Y/k0).
% Output: gridLs is the sequence of effective lambda values for the scaled
%            model

if  ~exist('output', 'var') || isempty(output)
    output = 0;
end

if ~exist('rem_ctrl', 'var')
    rem_ctrl.zeroRem = 1;
    rem_ctrl.complexRem = 1;
    rem_ctrl.redundRem = 1;
end

% First, Compute the solution path (2-dim)
oriX = X; oriY = Y;
%     if exist('screening', 'var')
%         if screening
%             X = X(:, predScreened);
%         end
%     end
[n, m] = size(Y); d = size(X, 2);
k0 = norm(X, 2);
X = X / k0; Y = Y / k0; % We added the Y scaling to make sure the grid for lambda (tuning parameter) can be somewhat universal

if ~exist('probType', 'var') || isempty(probType) || strcmp(probType, 'SEL') % problem type: SEL-RRR/PCA
    vecBased = 0;
else % SPA RRR/PCA
    vecBased = 1;
end
% Make a grid for lambda (sparsity parameter)
if exist('lambdasGiven', 'var') && ~isempty(lambdasGiven)
    if strcmp(regType, 'pen')
        lambdas = lambdasGiven /k0^2; %lambdasGiven; % ?? lambdas * norm(original full X, 2)^2 /k0^2;
    else % constrained type problem ('constr')
        lambdas = lambdasGiven;
    end
else
    if ~exist('probType', 'var') || isempty(probType) || strcmp(probType, 'SEL') % problem type: SEL-RRR/PCA
        if ~exist('regType', 'var') || isempty(regType) || strcmp(regType, 'pen') % regularization type: penalty or constraint
            uBnd = max(sqrt(sum((X'*Y)'.^2)));
            lLogBnd = -15; %-10; % -8; %-6; %-1.5; %
            lambdas = 2.^[lLogBnd:(log2(uBnd)-lLogBnd)/nPoints:log2(uBnd)];
            %     if uBnd < 1 % 2^(-1.5) > uBnd
            %         lambdas = 0:uBnd/nPoints:(uBnd+uBnd/nPoints); % 1  %0.01:.005:2
            %         lambdas = lambdas(2:end); % remove zero
            %     end
            lambdas = lambdas(end:-1:1);
        else
%             uBnd = round(size(Y, 1) * .8); % assume support size less than 0.8 times the sample size;
            %                 uBnd = round(size(X, 2)*.8);
            %                 lBnd2 = min(uBnd, round(nPoints/2) + 1);
            %                 lambdas = [1:(lBnd2-1)];
            %                 lambdas = [lambdas, round(lBnd2:max(1, (uBnd-lBnd2)/(nPoints-lBnd2-1)):uBnd)];
            
            % new way (assuming one is only interested in small cardinality models) -- subject to change if a dense (low rank) model is the goal
%             lBnd2 = min(uBnd, nPoints);
            lBnd2 = nPoints;

            lambdas = [1:(lBnd2)];
            
            %                 disp(lambdas)
            
        end
    else
        if ~exist('regType', 'var') || isempty(regType) || strcmp(regType, 'pen') % regularization type: penalty or constraint
            
            uBnd = max(max(abs(X'*Y))) * 4;
            lLogBnd = -7; %-15; %-10; % -6; %-1.5; %
            lambdas = 2.^[lLogBnd:(log2(uBnd)-lLogBnd)/nPoints:log2(uBnd)];
            lambdas = lambdas(end:-1:1);
        else
            uBnd = round(size(Y, 1) * size(Y, 2) * .25); % assume support size less than half of the total # of observations
            
            lBnd = 1;
            lBnd2 = min(uBnd, round(nPoints/2) + 1);
            
            %             lLogBnd = 0;
            %             lambdas = round( 2.^[lLogBnd:max(0, (log2(uBnd)-lLogBnd)/nPoints):log2(uBnd)] );
            % %             lambdas = [lBnd:max(1, (uBnd-lBnd)/nPoints):uBnd];
            lambdas = [1:(lBnd2-1)];
            lambdas = [lambdas, round(lBnd2:max(1, (uBnd-lBnd2)/(nPoints-lBnd2-1)):uBnd)];
            %             lambdas = lambdas(end:-1:1);
        end
    end
end

if rankGiven > 0
    candranks = rankGiven; %         candranks = initR;
else
    if rankGridUBnd > 0
        candranks = 1:rankGridUBnd; %
    else
        candranks = 1: (min([size(X), size(Y)])/1);
    end
    %         candranks = candranks(end:-1:1);
end

Bs = zeros(d, m, numel(lambdas)*numel(candranks));
rankBs = - ones(1, size(Bs, 3));
if ~exist('probType', 'var') || isempty(probType) || strcmp(probType, 'SEL') % problem type: SEL-RRR/PCA
    grpBs = - ones(1, size(Bs, 3));
else
    grpBs = cell(size(Bs,3), 1);
    Ss = cell(size(Bs,3), 1);
    Vs = cell(size(Bs,3), 1);
end
gridLs = rankBs;
numBs = 1;
Best = Bs(:, :, 1);


if full2d == 1 % full 2D:
    %for each k, evaluate the constrained rank grp-sel
    %estimate for any given Lambda
    for k = candranks
        for lambdaInd = 1:numel(lambdas)
            Lambda = lambdas(lambdaInd);
            if rankGiven > 0 && rankGiven == min([size(X), size(Y)]) % just the group penalty; no rank constraint
                [Best, nIter] = Func_MultiGroupSel(X, Y, Lambda, grpSelMaxIT, thresholdingWay, [], [], [], vecBased);
                rankEst = sum(svd(Best)>1e-4);
                rankBs(numBs) = rankEst;
                if ~exist('probType', 'var') || isempty(probType) || strcmp(probType, 'SEL') % problem type: SEL-RRR/PCA
                    grpBs(numBs) = sum(sqrt(sum(Best.^2, 2)) > 1e-4);
                else
                    grpBs{numBs} = Func_GetIndvGrps_w(Best);
                    Ss{numBs} = Shat;
                    Vs{numBs} = Vhat;
                end
                
            else
                
                % The usual  manner (including the quantile version)
                %                     [Best, Shat, Vhat, nIter] = Func_RankConstrGroupSel(X, Y, Lambda, k, maxIT, grpSelMaxIT, thresholdingWay);
                [Best, Shat, Vhat, nIter, remainingDims, rankEst] = ...
                    Func_RankConstrGroupSel2(X, Y, Lambda, k, errBnd, maxIT, grpSelMaxIT, thresholdingWay, Nu, [], [], vecBased);
                % did not supply Nu, BInit, and K. In particular, X has been scaled down by k0 and so there's no need to
                % supply the value of K (default 1).
                
                %                     % The progressive quantile manner (for constrained problems only)
                %                             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                             %%%%%%%%% We can replace the above by a progressive screening (still constrained) %%%%%%%%%
                %                             %%% By default: the first stage has probType = 'SEL' and then 'SPA' in the second stage
                %                             firstStageScrn = 1; secondStageScrn = 0;
                %                             qVal1 = Lambda; qVal2 = Lambda/2; % reduce the # of nonzeros each S column to half of its size
                %                                 % This q-value controls the total number of selected dimensions in SEL RRR/PCA,
                %                                 % or the AVERAGE number of nonzeros in each S vector in SPA RRR/PCA
                %                                 % (and thus the total number of nonzeros is q * rank)
                %                             pcamanner = 0; % 1: Use the transpose form to handle SPA PCA
                %                             squeezing = 1; % 0  %1 % 2 %1 %
                %                                     % squeezing: 2, 1, or  0:  0: no squeezing; 1: for RRR problems; 2: for PCA problems thru self-regression
                %                                         % Adding the squeezing operations can greatly reduces the
                %                                         % computational cost in high dimensions
                %                             coolingmanner = 'sigmoid'; %'fixed'; %'sigmoid-GivenCoolingTime' %'poly-GivenCoolingTime' %   'poly-0.5' %'inverse' % 'fixed' % 'dec1' % 'halving' %
                %                                 % The manner of decreasing the 'temperature'. Comparision in terms of cooling speed: 'fixed>>'halving'>>'inverse'>>'dec1'
                %                             coolingIndex =  .1; % 2; % 0.5; % 2; %1;
                %                                                    %  Note: 1 is fast yet not accurate, 0.01 seems to be accurate but a bit slow
                %                             coolingTime = 500; %50; % NaN %
                %                                 % Given the number of cooling steps (cooling time), we sample the cooling points along the curve
                %                             keptTime = -1; % 100; %
                %                                 % We let the algorithm run for a few more iterations with tempreture fixed at the target
                %                                 % With squeezing applied, this does not change the sparsity pattern.
                %                                 % But it helps getting more accurate estimates. When squeezing is not
                %                                 % applied, theoretically the sparsity pattern can still change.
                %                             % Number of iterations (optimization steps) performed given each temperature
                %                             maxIT_pro = 1; % this is to be used with progressive dim redu (thru squeezing)
                %                             grpSelMaxIT_pro = 1; % this is to be used with progressive dim redu (thru squeezing)
                %
                %                             [Best, Shat, Vhat, remainingDims, firstStageRemainingDims, secondStageRemainingDims] ...
                %                                 = Func_HybridSELSPA_RRR2(X, Y, firstStageScrn, secondStageScrn, squeezing, squeezing,...
                %                                 coolingmanner, coolingmanner, coolingIndex, coolingIndex, coolingTime, coolingTime, ...
                %                                 keptTime, keptTime, maxIT_pro, maxIT_pro, errBnd, errBnd, grpSelMaxIT_pro, grpSelMaxIT_pro, ...
                %                                 k, k, qVal1, qVal2, Nu, Nu, pcamanner, 0);
                %                             rankEst = sum(svd(Best)>1e-4);
                %                             nIter = [];
                %                             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ends here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                cand_ranks_grid(numBs) = k;
                rankBs(numBs) = rankEst;
                if ~exist('probType', 'var') || isempty(probType) || strcmp(probType, 'SEL') % problem type: SEL-RRR/PCA
                    grpBs(numBs) = numel(remainingDims); % sum(sqrt(sum(Best.^2, 2)) > 1e-4);
                else
                    grpBs{numBs} = Func_GetIndvGrps_w(Shat);
                    Ss{numBs} = Shat;
                    Vs{numBs} = Vhat;
                end
            end
            
            if output
                if ~exist('probType', 'var') || isempty(probType) || strcmp(probType, 'SEL') % problem type: SEL-RRR/PCA
                    disp(['Lambda=', num2str(Lambda), ', k=', num2str(k), '; num(outer)ITs ', num2str(nIter), ' -- rank=', num2str(rankBs(numBs)), ', grp=', num2str(grpBs(numBs))])
                else
                    disp(['Lambda=', num2str(Lambda), ', k=', num2str(k), '; num(outer)ITs ', num2str(nIter), ' -- rank=', num2str(rankBs(numBs)), ', grp=', num2str(sum(grpBs{numBs}))])
                end
            end
            
            Bs(:, :, numBs) = Best;
            gridLs(numBs) = Lambda;
            numBs = numBs + 1;
            % Pre-Termination: If number of selected groups is larger than the rank of X times 0.8, no
            % need to decrease lambda anymore. More groups result in the same prediction (and same DF in PIC)
            % [Valid for nonconvex pens???]
            if ~exist('probType', 'var') || isempty(probType) || strcmp(probType, 'SEL') % problem type: SEL-RRR/PCA
                if grpBs(numBs-1) > min(size(X)) * 1
                    %                         if (min(min(size(X)), grpBs(numBs-1)) + size(Y, 2) - rankBs(numBs-1)) * rankBs(numBs-1) > numel(Y) * .6 % grpBs(numBs-1) > size(X, 1) * .8
                    %&& strcmp(lower(thresholdingWay), 'soft')
                    
                    break;
                end
            else
                if  sum(grpBs{numBs-1}) > size(Y,1) * size(Shat,2) * .8 %&& strcmp(lower(thresholdingWay), 'soft')
                    break;
                end;
            end
        end
        
    end
else    % adaptive 2D:
    % for each lambda, evaluate the estimate for meaningful k's
    % starting from the largest
    candranks = candranks(end:-1:1);
    for lambdaInd = 1:numel(lambdas)
        Lambda = lambdas(lambdaInd);
        k = candranks(1);
        while  sum(candranks == k)
            %                 [Best, Shat, Vhat, nIter] = Func_RankConstrGroupSel(X, Y, Lambda, k, maxIT, grpSelMaxIT, thresholdingWay);
            [Best, Shat, Vhat, nIter, remainingDims, rankEst] = ...
                Func_RankConstrGroupSel2(X, Y, Lambda, k, errBnd, maxIT, grpSelMaxIT, thresholdingWay, Nu, [], [], vecBased);
            
            k = rankEst-1; %note that rankBs(numBs) <= k is always true
            if rankEst > 0 % otherwise Best = 0. Do not record.
                Bs(:, :, numBs) = Best;
                rankBs(numBs) = rankEst;
                if ~exist('probType', 'var') || isempty(probType) || strcmp(probType, 'SEL') % problem type: SEL-RRR/PCA
                    grpBs(numBs) = numel(remainingDims); % sum(sqrt(sum(Best.^2, 2)) > 1e-4);
                else
                    grpBs{numBs} = Func_GetIndvGrps_w(Shat);
                    Ss{numBs} = Shat;
                    Vs{numBs} = Vhat;
                end
                cand_ranks_grid(numBs) = k;
                gridLs(numBs) = Lambda;
                numBs = numBs + 1;
            end
        end
    end
end

%     Bs = Bs / k0; %no need, b/c both X and Y are scaled down by k0
X = oriX; Y = oriY;

numBs = numBs - 1;
Bs = Bs(:, :, 1:numBs);
if exist('Ss', 'var') && exist('Vs', 'var') % problem type: SPA-RRR/PCA
    Ss = Ss(1:numBs); Vs = Vs(1:numBs);
end
rankBs = rankBs(1:numBs);
grpBs = grpBs(1:numBs);
gridLs = gridLs(1:numBs);

complexRem = rem_ctrl.complexRem;
redundRem = rem_ctrl.redundRem;
zeroRem = rem_ctrl.zeroRem;
% complexRem = 0; % 0; % %0: Sometimes when lambdas are given, it's better not to remove them (otherwise you'll get []
% redundRem = 0;
% zeroRem = 0;  % 1; %

if complexRem > 0
    % Make sure the estimates have the # of groups less than the rank
    % of X (or DF control??)
    if ~exist('probType', 'var') || isempty(probType) || strcmp(probType, 'SEL') % problem type: SEL-RRR/PCA
        validModelInds = find(grpBs <= min(size(X)));
        %                 validModelInds = find( (min(min(size(X)), grpBs) + size(Y, 2) - rankBs) .* rankBs < numel(Y) * .6 );
    else
        validModelInds = find( cellfun(@(x) sum(x) <= numel(Y), grpBs) );
        error('To be changed')
    end
    numBs = numel(validModelInds);
    Bs = Bs(:, :, validModelInds);
    if exist('Ss', 'var') && exist('Vs', 'var') % problem type: SPA-RRR/PCA
        Ss = Ss(validModelInds); Vs = Vs(validModelInds);
    end
    rankBs = rankBs(validModelInds);
    grpBs = grpBs(validModelInds);
    gridLs = gridLs(validModelInds);
end
if redundRem > 0
    % Remove redundant Bs
    redundantInds = find(sum(sum(abs(Bs(:,:,2:end)-Bs(:,:,1:end-1)))) <1e-4)+1;
    if ~exist('probType', 'var') || isempty(probType) || strcmp(probType, 'SEL') % problem type: SEL-RRR/PCA
        redundantInds(grpBs(redundantInds) ~= grpBs(redundantInds-1)) = [];
    end
    Bs(:, :, redundantInds) = [];
    if exist('Ss', 'var') && exist('Vs', 'var') % problem type: SPA-RRR/PCA
        Ss(redundantInds) = []; Vs(redundantInds) = [];
    end
    
    rankBs(redundantInds) = [];
    grpBs(redundantInds) = [];
    gridLs(redundantInds) = [];
    numBs = numel(rankBs);
end
if zeroRem > 0
    % Remove the zero rank est
    zerorankInds = find(rankBs == 0);
    Bs(:, :, zerorankInds) = [];
    if exist('Ss', 'var') && exist('Vs', 'var') % problem type: SPA-RRR/PCA
        Ss(zerorankInds) = []; Vs(zerorankInds) = [];
    end
    rankBs(zerorankInds) = [];
    grpBs(zerorankInds) = [];
    gridLs(zerorankInds) = [];
    numBs = numel(rankBs);
end
%     if exist('screening', 'var')
%         if screening
%             origBs = zeros(size(X, 2), size(Y, 2), size(Bs, 3));
%             origBs(predScreened, :, :) = Bs;
%             Bs = origBs;
%             clear origBs;
%         end
%     end

clear m n d;
nzPatts = -ones(numBs, size(X, 2));
trnerrs = - ones(1, numBs);
for indB = 1:numBs
    %         if ~exist('probType', 'var') || isempty(probType) || strcmp(probType, 'SEL') % problem type: SEL-RRR/PCA
    %         nzPatts(indB, :) = sqrt(sum((Bs(:,:,indB))'.^2)) > (1e-4 /k0);
    nzPatts(indB, :) = sqrt(sum((Bs(:,:,indB))'.^2, 1)) > (1e-4); % no need for back scaling because X and Y are both scaled down by k0
    %         else
    %             nzPatts(indB, :) = abs(Bs(:,:,indB)) > (1e-4);
    %         end
    trnerrs(indB) = norm(Y-X*Bs(:,:,indB),'fro')^2;
end
if ~exist('Ss', 'var') || ~exist('Vs', 'var')
    Ss = []; Vs = [];
end

end
