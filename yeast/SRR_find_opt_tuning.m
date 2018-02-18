function [ opt_ctrl ] = SRR_find_opt_tuning( opt_cri, path_ctrl, type_ctrl )

%find_opt_tuning - after running a model selection or information criteria
%function (such as SRR_SF_PIC_fractional) use this to select the optimal model.
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%    opt_cri - the output from a model selection/information criteria
%    function
%    path_ctrl - the output from SRR_solution_path
%    type_ctrl - a struct with problem type information
%
% Outputs:
%    opt_ctrl - struct with optimal values

if nargin < 2
    type_ctrl = struct;
end

if ~isfield(type_ctrl, 'probType')
    type_ctrl.probType = 'SEL';
end

%%%initialize from tuning_ctrl struct
Bs = path_ctrl.Bs;
card_Bs = path_ctrl.card_Bs;
rank_Bs = path_ctrl.rank_Bs;
Ss = path_ctrl.Ss;
Vs = path_ctrl.Vs;
%%%

% Find the optimal estimate(s)
optInd = find(opt_cri == min(opt_cri), 1, 'first');
BOpt = Bs(:, :, optInd);

if strcmp(type_ctrl.probType, 'SEL') % problem type: SEL-RRR/PCA
    card_opt = card_Bs(optInd);
else
    card_opt = sum(sqrt(sum(BOpt.^2, 2)) > 1e-4); %sum(nzPatts(optInd, :)); % card_Bs(optInd);
    SOpt = Ss{optInd}; VOpt = Vs{optInd};
end

rank_opt = rank_Bs(optInd);

opt_ctrl.BOpt = BOpt;
opt_ctrl.card_opt = card_opt;
opt_ctrl.rank_opt = rank_opt;
opt_ctrl.opt_ind = optInd;
opt_ctrl.opt_cri = opt_cri(optInd);
opt_ctrl.nz_patt = path_ctrl.nz_patts(optInd, :);

if ~ strcmp(type_ctrl.probType, 'SEL')
    opt_ctrl.SOpt = SOpt;
    opt_ctrl.VOpt = VOpt;
end

% output_ctrl = tuning_ctrl;
% output_ctrl.card_opt = card_opt;
% output_ctrl.rank_opt = rank_opt;
% output_ctrl.opt_ind = optInd;
% 
% if ~ strcmp(type_ctrl.probType, 'SEL')
%     output_ctrl.SOpt = SOpt;
%     output_ctrl.VOpt = VOpt;
% end


end

