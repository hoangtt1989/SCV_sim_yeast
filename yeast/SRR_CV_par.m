function [opt_ctrl, type] ...
    = SRR_CV_par(input_ctrl, CVSplit_ctrl, conv_ctrl, type_ctrl)

%SRR_CV - do plain cross-validation for SRR

% Inputs:
%    input_ctrl - a struct with X, Y matrices
%    CVSplit_ctrl - struct output from CVSplit
%    conv_ctrl - a struct with convergence parameters
%    type_ctrl - a struct with problem type information
%
% Outputs:
%    opt_ctrl - this is a bit different than the output from
%    SRR_find_opt_tuning because it includes the actual grids of ranks and
%    cardinalities as well as the optimal ones. After running CV then we
%    refit the model on the whole data using the positions in the candidate
%    grids.

if nargin < 3
    conv_ctrl = struct;
end
if nargin < 4
    type_ctrl = struct;
end


%%%%%%these get passed to the solution path functions - don't remove any
%%%%%%models or else we won't be able to average over all splits
conv_ctrl.complexRem = 0;
conv_ctrl.redundRem = 0;
conv_ctrl.zeroRem = 0;
%%%%%%

%%%%initialize from CVSplit_ctrl struct
nCV = CVSplit_ctrl.nCV;
dataIndsCV = CVSplit_ctrl.dataIndsCV;
dataIndsCVStarts = CVSplit_ctrl.dataIndsCVStarts;
dataIndsCVEnds = CVSplit_ctrl.dataIndsCVEnds;
%%%%
%%%%initialize from input_ctrl struct
X = input_ctrl.X;
Y = input_ctrl.Y;
%%%%

if ~isfield(conv_ctrl, 'card_upper_bnd')
    conv_ctrl.card_upper_bnd = 50; %%%the maximum cardinality of predictors for the grid (constrained form)
end
if ~isfield(conv_ctrl, 'rank_upper_bnd')
    conv_ctrl.rank_upper_bnd = 50; %%%the maximum rank number for the grid
end

num_params = conv_ctrl.card_upper_bnd * conv_ctrl.rank_upper_bnd;

% 2. Compute the CV-error for each model. Note the CV split is the same for different candidate patterns.
oriX = X; oriY = Y;

valerrs = zeros(nCV, num_params);
rank_CV = zeros(nCV, num_params);
card_CV = zeros(nCV, num_params);

cand_cards = cell(nCV, 1);
cand_ranks = cell(nCV, 1);

parpool
parfor indCV = 1:nCV
    valInds = dataIndsCV( dataIndsCVStarts(indCV):dataIndsCVEnds(indCV) );
    tmpInds = 1:size(oriY, 1); tmpInds(valInds) = [];
    trnX = oriX(tmpInds, :); trnY = oriY(tmpInds, :);
    valX = oriX(valInds,:); valY = oriY(valInds,:);
    
    trn_ctrl = struct;
    trn_ctrl.X = trnX;
    trn_ctrl.Y = trnY;
    
    trn_output_ctrl = struct;
    trn_ctrl2 = struct;
    trn_ctrl2 = trn_ctrl;
    conv_ctrl2 = struct;
    conv_ctrl2 = conv_ctrl;
    type_ctrl2 = struct;
    type_ctrl2 = type_ctrl;
    trn_output_ctrl = SRR_solution_path(trn_ctrl2, conv_ctrl2, type_ctrl2);
    cand_cards{indCV} = trn_output_ctrl.cand_cards;
    cand_ranks{indCV} = trn_output_ctrl.cand_ranks;
%     num_params = size(trn_output_ctrl.Bs, 3);
    
    valerrTmp = zeros(num_params, 1);
    
    for mods = 1:num_params
        valerrTmp(mods) = norm(valY - valX * trn_output_ctrl.Bs(:, :, mods), 'fro')^2;
    end
    
    valerrs(indCV, :) = valerrTmp;
    rank_CV(indCV, :) = trn_output_ctrl.rank_Bs;
    card_CV(indCV, :) = trn_output_ctrl.card_Bs;

end
delete(gcp('nocreate'))

opt_cri = mean(valerrs);


% Find the optimal estimate(s)
optInd = find(opt_cri == min(opt_cri), 1, 'first');
opt_ctrl.opt_ind = optInd;
opt_ctrl.opt_cri = opt_cri;
opt_ctrl.min_cri = opt_cri(optInd);

%%%%the actual grid
opt_ctrl.rank_CV = rank_CV;
opt_ctrl.card_CV = card_CV;
%%%%
%%%%the optimal values in the candidate grid
% opt_ctrl.cand_cards = trn_output_ctrl.cand_cards;
% opt_ctrl.cand_ranks = trn_output_ctrl.cand_ranks;
% opt_ctrl.cand_card_opt = trn_output_ctrl.cand_cards(optInd);
% opt_ctrl.cand_rank_opt = trn_output_ctrl.cand_ranks(optInd);
opt_ctrl.cand_cards = cand_cards{1};
opt_ctrl.cand_ranks = cand_ranks{1};
opt_ctrl.cand_card_opt = opt_ctrl.cand_cards(optInd);
opt_ctrl.cand_rank_opt = opt_ctrl.cand_ranks(optInd);
%%%%

%%%%rerun on the full data
conv_ctrl.rank_given = opt_ctrl.cand_rank_opt;
conv_ctrl.card_given = opt_ctrl.cand_card_opt;
run_output_ctrl = SRR_solution_path(input_ctrl, conv_ctrl, type_ctrl);

opt_ctrl.card_opt = run_output_ctrl.card_Bs;
opt_ctrl.rank_opt = run_output_ctrl.rank_Bs;
opt_ctrl.BOpt = run_output_ctrl.Bs;
opt_ctrl.nz_patt = run_output_ctrl.nz_patts;
%%%%

type = 'CV';


end

