function [ opt_ctrl ] = SRR_opt_sel( path_ctrl, type, CVSplit_ctrl, type_ctrl )
%SRR_opt_sel - after getting a solution path, tune using a criteria then
%get the optimal B matrix. A wrapper around the various criteria functions
%and SRR_find_opt_tuning
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%    path_ctrl - the output from SRR_solution_path
%    type - a string specifying the type of criteria
%    CVSplit_ctrl - used for SRR_SCV
%    type_ctrl - struct (optional)
%
% Outputs:
%    opt_ctrl



if nargin < 4
    [n, ~] = size(path_ctrl.X);
    CVSplit_ctrl = CVSplit(n, 5); %%%5 fold CV by default
end
if nargin < 5
    type_ctrl = struct;
end

if strcmp(type, 'SF_PIC_fractional') || strcmp(type, 'AIC') || strcmp(type, 'BIC') || strcmp(type, 'EBIC') || strcmp(type, 'SCV plugin') || strcmp(type, 'SCV fractional') || strcmp(type, 'SCV none')
    
    switch type
        case 'SF_PIC_fractional'
            [opt_cri, type] = SRR_SF_PIC_fractional(path_ctrl, type_ctrl);
        case 'AIC'
            [opt_cri, type] = SRR_AIC(path_ctrl);
        case 'BIC'
            [opt_cri, type] = SRR_BIC(path_ctrl);
        case 'EBIC'
            [opt_cri, type] = SRR_EBIC(path_ctrl);
        case 'SCV plugin'
            [opt_cri, type] = SRR_SCV(path_ctrl, CVSplit_ctrl, 'plugin', type_ctrl);
        case 'SCV fractional'
            [opt_cri, type] = SRR_SCV(path_ctrl, CVSplit_ctrl, 'pic', type_ctrl);
        case 'SCV none'
            [opt_cri, type] = SRR_SCV(path_ctrl, CVSplit_ctrl, 'none', type_ctrl);
    end
    
else
    
    error('type not implemented')
    
end


opt_ctrl = SRR_find_opt_tuning(opt_cri, path_ctrl, type_ctrl);
opt_ctrl.type = type;


end

