function [opt_cri, type] ...
    = SRR_SCV(path_ctrl, CVSplit_ctrl, biascorrt, type_ctrl)

%%%%%%%previous output (changed)
% BOpt, tuningparamOpt, ROpt, GOpt, ZOpt, UOpt, optInd, Cest, Us, Zs, Rs, Js, valerrs, valerrOpt
%%%%%%%


% Selective Projective Cross-Validation (SPCV)
% Do not perform bias-correction for the give estimates, because we only need the sparsity pattern,
% and the range space (if the estimate does not have full rank) in SPCV.
% Intercepts not included in Gaussian models.
% If pcamanner ~=0, extract the PCA loading vectors to construct the new
% design, use the PCA estimate at each training.
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%    path_ctrl - struct output from SRR_solution_path
%    CVSplit_ctrl - struct output from CVSplit
%    biascorrt - string specifying bias correction type
%    type_ctrl - struct with fields for pcamanne
%
% Outputs:
%    opt_cri

%%%%initialize from CVSplit_ctrl struct
nCV = CVSplit_ctrl.nCV;
dataIndsCV = CVSplit_ctrl.dataIndsCV;
dataIndsCVStarts = CVSplit_ctrl.dataIndsCVStarts;
dataIndsCVEnds = CVSplit_ctrl.dataIndsCVEnds;
%%%%
%%%%initialize from path_ctrl struct
Bs = path_ctrl.Bs;
tuningparams = path_ctrl.tune_params;
nzPatts = path_ctrl.nz_patts;
rankBs = path_ctrl.rank_Bs;
X = path_ctrl.X;
Y = path_ctrl.Y;
%%%%

if ~ isfield(type_ctrl, 'pcamanner')
    type_ctrl.pcamanner = 0;
end
if nargin < 3
    biascorrt = 'plugin'; %%'fractional', 'none' implemented
end


% 1. Extract the candidate sparsity and/or low rank sparsity patterns,
%    Construct the new design matrix for each estimate, and record the
%    corresponding transform matrix.
% [n, d] = size(X);
% m = size(Y, 2);
numBs = numel(tuningparams);
Zs = cell(1, numBs); % Use a cell array to store the new model matrices
Us = cell(1, numBs); % Use a cell array to store the transofrm matrices
% Cs = cell(1, numBs); % Use a cell array to store the correspoinding coefficient ests under the new design
Rs = zeros(1, numBs);% Rank of each model
Js = zeros(1, numBs);% support size of each model
Cests = cell(1, numBs); % Use a cell array to store the estimates (using all the data) for the new design. Usually empty, not so when pcammner~=0
for tInd = 1:numBs
    B = Bs(:, :, tInd);
    nzPatt = nzPatts(tInd, :);
    rankB = rankBs(tInd);
    
    %         [newdesign, newDim, U2, supp, rank] = PatternExtr(X, B);
    [newdesign, newDim, U2, supp, rank, tmppatt, Cest] = PatternExtr(X, B, type_ctrl.pcamanner, nzPatt, rankB) ;
    
    Zs{tInd} = newdesign;
    Us{tInd} = U2;
    Rs(tInd) = rank;
    Js(tInd) = supp;
    Cests{tInd} = Cest;
end

% 2. Compute the CV-error for each model. Note the CV split is the same for different candidate patterns.
oriX = X; oriY = Y;
valerrs=zeros(size(tuningparams));
% trnerr = zeros(size(tuningparams)); %store trnerrs for plug in form
q = sum(svd(oriX) > 1e-4); % rank(oriX);
for paramInd = 1:numBs
    Z = Zs{paramInd};
    valerrTmp = 0;
    for indCV = 1:nCV
        valInds = dataIndsCV( dataIndsCVStarts(indCV):dataIndsCVEnds(indCV) );
        tmpInds = 1:size(oriY, 1); tmpInds(valInds) = [];
        X = Z(tmpInds, :); Y = oriY(tmpInds, :);
        valX = Z(valInds,:); valY = oriY(valInds,:);
        if type_ctrl.pcamanner == 0
            B_est = pinv(X' * X) * (X'*Y);
        else
            B_est = Cests{paramInd};
        end
        valerrTmp = valerrTmp + norm( valY - valX * B_est, 'fro')^2;
    end
    %     valerrs(paramInd) = valerrTmp / numel(oriY); % Make sure the denominator does not affect the final valerr
    valerr_scaled = valerrTmp / numel(oriY);
    %do OLS on overall data for plugin form
    B_all = pinv(Z' * Z) * (Z' * oriY);
    trnerr = norm(oriY - Z * B_all, 'fro')^2;
    trnerr_scaled = trnerr / numel(oriY);
    
    r = Rs(paramInd);
    J = Js(paramInd);
    if type_ctrl.pcamanner == 0 % The complexity for a generaral regression setup: Y = X B + E with E iid (sub)Gaussian
        %complexity = 1 * ( A1 * min(q, J) * r + A3*(size(oriY, 2) - r)*r+ A2 * J* (1+ log(size(oriX, 2) / J)) );
        %valerrs(paramInd) = numel(oriY) * log(valerrs(paramInd)) + complexity;
        
        
        %%%%%%alpha1 alpha2 values
        if strcmp(biascorrt, 'fractional') || strcmp(biascorrt, 'none')
            alpha1 = 2; %%5 fold
            alpha2 = 2.4; %%5 fold
        elseif strcmp(biascorrt, 'plugin') || strcmp(biascorrt, 'none')
            alpha1 = 4.6; %%5 fold
            alpha2 = 3.5; %%5 fold
        else
            error('Supply valid value for biascorrt')
        end
        %%%%%%
        
        R_term = (min(q, J) - r) * r;
        IF_term = J * (1 + log(size(oriX, 2)) - log(J));
        complexity_full = ( alpha1 * (R_term + size(oriY, 2) * r) + alpha2 * IF_term ) / numel(oriY) ;
        complexity_part = ( alpha1 * R_term + alpha2 * IF_term ) / numel(oriY) ;
        complexity = complexity_part;
        
        %%%%%%%%%%HT: Rate correction from the paper
        if strcmp(biascorrt, 'fractional') || strcmp(biascorrt, 'none')
            if strcmp(biascorrt, 'fractional')
%                 alpha1 = 2;
%                 alpha2 = 2.4;
                complexity = 1 - complexity;
                if complexity_full >= 1
                    valerrs(paramInd) = inf;
                else
                    valerrs(paramInd) = numel(oriY) * valerr_scaled / complexity;
                end
            elseif strcmp(biascorrt, 'none')
                %%%%%%%%%%%%HT: this is for plain structural CV
                valerrs(paramInd) = numel(oriY) * valerr_scaled;
            end
        end
        
        %plug in form
        if strcmp(biascorrt, 'plugin')
            if complexity_full >= 1
                valerrs(paramInd) = inf;
            else
                valerrs(paramInd) = numel(oriY) * (valerr_scaled + trnerr_scaled * complexity);
            end
        end
        
    else % The complexity for PCA problems: Y = X B + E V^T.
        complexity = 1 * ( A1 * min(q, J)*r + A3*(size(oriY, 2) - 2 * r)*r+ A2 * J*(1+ log(size(oriX, 2) / J)) );
        valerrs(paramInd) = (size(oriY,1)*(size(oriY,2) - r)) * log(valerr_scaled) + complexity;
    end

end

opt_cri = valerrs;
type = strcat('SCV ', biascorrt);

% % 3. Evaluate the optimal one by minimizing the criterion.
% % figure; plot(log(1+valerrs), 'rx'); hold on, plot(log(1+valerr_scaled), 'b.')
% optInd = find(valerrs == min(valerrs), 1, 'first');
% tuningparamOpt = tuningparams(optInd);
% BOpt = Bs(:, :, optInd);
% ZOpt = Zs{optInd};
% UOpt = Us{optInd};
% valerrOpt = valerrs(optInd);
% ROpt = Rs(optInd);
% GOpt = Js(optInd);
% 
% 
% % 4. Compute the (nonpenalized) multivaraite estimate under the optimal design
% X = ZOpt; Y = oriY;
% if type_ctrl.pcamanner == 0
%     Cest = pinv(X' * X) * (X'*Y);
% else
%     Cest = Cests{optInd};
% end


end

