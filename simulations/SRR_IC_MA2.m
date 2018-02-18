function [ opt_ctrl ] = SRR_IC_MA2( path_ctrl, type )
%SRR_IC_MA - use model averaging on AIC or BIC to get a final test error or
%training error
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%    path_ctrl - the output from SRR_solution_path
%
% Outputs:
%    opt_ctrl - a struct with performance information

if nargin < 2
    type = 'AIC';
end

switch type
    case 'AIC'
        opt_cri = SRR_AIC(path_ctrl);
    case 'BIC'
        opt_cri = SRR_BIC(path_ctrl);
end

m = size(path_ctrl.Y, 2);
p = size(path_ctrl.X, 2);

num_tune = numel(path_ctrl.tune_params);

wts = exp(-.5 * opt_cri) ./ sum(exp(-.5 * opt_cri));

preds = zeros(p, m);
for i = 1:num_tune
    
    preds = preds + wts(i) * path_ctrl.Bs(:, :, i);

end


opt_ctrl.BOpt = preds;
opt_ctrl.card_opt = -1000;
opt_ctrl.rank_opt = -1000;
opt_ctrl.opt_ind = -1000;
opt_ctrl.opt_cri = -1000;
opt_ctrl.nz_patt = -1000;

end

