function [ opt_cri, type ] = SRR_EBIC( path_ctrl )
%SF_EBIC_fractional - for a given solution path compute the values of
%EBIC. path_ctrl is the output_ctrl from
%SRR_solution_path
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%    path_ctrl - a struct with solution path information
%
% Outputs:
%    opt_cri

if ~isfield(path_ctrl, 'q')
    path_ctrl.q = sum(svd(path_ctrl.X) > 1e-4);
end


%%%initialize from path_ctrl struct
X = path_ctrl.X;
Y = path_ctrl.Y;
q = path_ctrl.q;
%%%
%%%initialize from path_ctrl struct
rank_Bs = path_ctrl.rank_Bs;
card_Bs = path_ctrl.card_Bs;
trn_errs = path_ctrl.trn_errs;
%%%

penDF = (min(q, card_Bs) + size(Y, 2) - rank_Bs).*rank_Bs;
logterm = numel(Y) * log(trn_errs / numel(Y));


opt_cri = logterm + penDF * (log(numel(Y)) + 2*log(size(X, 2)));

type = 'EBIC';

end

