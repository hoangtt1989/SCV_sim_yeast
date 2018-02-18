function [ output_ctrl, conv_ctrl, type_ctrl ] = SRR_solution_path( input_ctrl, conv_ctrl, type_ctrl )

%RR_solution_path - for a given Y, X, finds the sparse and reduced rank
%regression solution path
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%    input_ctrl - a struct with X, Y matrices
%    type_ctrl - a struct with problem type information
%    conv_ctrl - a struct with convergence parameters
%
% Outputs:
%    output_ctrl
%    type_ctrl
%    conv_ctrl
% Requires:
%    Func_SRRRPathComp2
%    Func_EstsCorrt
%    FuncMultiGroupSel
%    Func_GetIndvGrps_w
%    Func_RankConstrGroupSel2

if nargin < 2
    conv_ctrl = struct;
end
if nargin < 3
    type_ctrl = struct;
end

rem_ctrl = struct;


if ~isfield(type_ctrl, 'modelType')
    type_ctrl.modelType = 'rrr';
end
if ~isfield(type_ctrl, 'probType')
    type_ctrl.probType = 'SEL';
end
if ~isfield(type_ctrl, 'regType')
    type_ctrl.regType = 'constr'; %%%code is designed for constrained form
end
if ~isfield(conv_ctrl, 'maxIT')
    conv_ctrl.maxIT = 1e+3;
end
if ~isfield(conv_ctrl, 'grpSelMaxIT')
    conv_ctrl.grpSelMaxIT = -1;
end
if ~isfield(conv_ctrl, 'errBnd')
    conv_ctrl.errBnd = 1e-4;
end
if ~isfield(conv_ctrl, 'full2d')
    conv_ctrl.full2d = 1;
end
if ~isfield(conv_ctrl, 'rank_given')
    conv_ctrl.rank_given = -1;
end
if ~isfield(conv_ctrl, 'card_upper_bnd')
    conv_ctrl.card_upper_bnd = 50; %%%the maximum cardinality of predictors for the grid (constrained form)
end
if ~isfield(conv_ctrl, 'rank_upper_bnd')
    conv_ctrl.rank_upper_bnd = 50; %%%the maximum rank number for the grid
end
if ~isfield(conv_ctrl, 'card_given')
    conv_ctrl.card_given = [];
end

%%%%%removing bad estimates
if ~isfield(conv_ctrl, 'complexRem')
    conv_ctrl.complexRem = 1;
end
if ~isfield(conv_ctrl, 'redundRem')
    conv_ctrl.redundRem = 1;
end
if ~isfield(conv_ctrl, 'zeroRem')
    conv_ctrl.zeroRem = 1;
end
rem_ctrl.complexRem = conv_ctrl.complexRem;
rem_ctrl.redundRem = conv_ctrl.redundRem;
rem_ctrl.zeroRem = conv_ctrl.zeroRem;
%%%%%

%%%initialize from input_ctrl
X = input_ctrl.X;
Y = input_ctrl.Y;
%%%

%%%checking inputs
if strcmp(type_ctrl.modelType, 'rrr')
    estCorrtWay = 'rrr'; %  'pca' %'none' %
    pcamanner = 0; %1; %
elseif strcmp(type_ctrl.modelType, 'pca-t') % This is the PCA manner 2 (recommended) by trasposing the data
    estCorrtWay = 'pca'; %'rrr' %  'none' %
    pcamanner = 1;
elseif strcmp(type_ctrl.modelType, 'pca-self') %PCA via self-regression
    estCorrtWay = 'pca'; %'rrr' %  'none' %
    pcamanner = 2; %0; %
end
%%%


%%%get the grid
% nPoints = s_true + 20; %20;%50; %100;  500; %
% rankGridUBnd = min(s_true, (2*r_true)); % r_true+1 % -1
%%%

%%%
if strcmp(type_ctrl.regType, 'pen') % regularization type: penalty or constraint
    thresholdingWay = 'hybrid'; % 'hard' %'soft' %  % hard is much better than soft which tends to select much more than necessary
    Nu = 0; %1e-4; 1e-2; %1e-4; % %  only used for hybrid
else % constraint
    thresholdingWay = 'hybrid-prop';
    Nu = 0; %1e-4; %1e-2; %   Nu: the ridge shrinkage parameter for the scaled model. That is, the
    %   problem is defined as  |Y - X B|_F^2/(2 K) + Nu^2 |B|_F^2/2 with sparsity/rank constraint.
end
%%%
% If the upper bound is s_true, from the data generation, the PC space being I must achieve the minimal test error

%%%
if pcamanner == 0 || pcamanner == 2
    Design = X; Resp = Y;
elseif pcamanner == 1
    Resp = X'; Design = eye(size(Resp,1));
else
    error('Wrong pcamanner value')
end
%%%



%%%%%%%%%%%%%%%%compute the solution path
[Bs, rankBs, grpBs, tuningparams, nzPatts, trnerrs, numBs, Ss, Vs, cand_ranks] ...
    = Func_SRRRPathComp2(Design, Resp, type_ctrl.probType, type_ctrl.regType, conv_ctrl.rank_given, conv_ctrl.rank_upper_bnd, conv_ctrl.card_given, conv_ctrl.full2d, conv_ctrl.card_upper_bnd, type_ctrl.modelType, conv_ctrl.grpSelMaxIT, conv_ctrl.errBnd, thresholdingWay, Nu, rem_ctrl);
%%%%%%%%%%%%%%%%


if pcamanner == 1
    Bs2 = zeros(size(Bs,2),  size(Bs,1), size(Bs,3));
    for tInd = 1:size(Bs,3)
        Bs2(:, :, tInd) = Bs(:, :, tInd)';
    end
    Bs = Bs2; clear Bs2;
end

% Estimate bias correction
if strcmp(type_ctrl.probType, 'SEL') % problem type: SEL-RRR/PCA
    %     estCorrtWay = 'pca' %'none' % 'rrr' %
    % Estimatimation  bias correction -- useful for, say, soft thresholding
    if strcmp(estCorrtWay, 'rrr')
        [Bs_corrt] = Func_EstsCorrt(Bs, X, Y, tuningparams, nzPatts, rankBs);
        Bs_orig = Bs;
        Bs = Bs_corrt;
    elseif strcmp(estCorrtWay, 'pca')
        [Bs_corrt] = Func_EstsCorrt_PCA(Bs, X, tuningparams, nzPatts, rankBs);
        Bs_orig = Bs;
        Bs = Bs_corrt;
    end
else
    estCorrtWay = [];
end


for indB = 1:numBs
    trnerrs(indB) = norm(Y-X*Bs(:,:,indB),'fro')^2; % we did not divide it by 2
end

if any(rankBs == 0)
    disp('Some ranks are zero') 
end

output_ctrl.Bs = Bs;
output_ctrl.rank_Bs = rankBs;
output_ctrl.card_Bs = grpBs;
output_ctrl.tune_params = tuningparams;
output_ctrl.nz_patts = nzPatts;
output_ctrl.trn_errs = trnerrs;
output_ctrl.X = X;
output_ctrl.Y = Y;
output_ctrl.Ss = Ss;
output_ctrl.Vs = Vs;
output_ctrl.cand_ranks = cand_ranks;
output_ctrl.cand_cards = tuningparams;

type_ctrl.pcamanner = pcamanner;
type_ctrl.estCorrtWay = estCorrtWay;

end

