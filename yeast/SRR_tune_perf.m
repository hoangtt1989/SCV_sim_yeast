function [ perf_ctrl ] = SRR_tune_perf( path_ctrl, test_ctrl, type, CVSplit_ctrl, type_ctrl  )
%SRR_tune_perf - after getting a solution path, tune using a criteria then
%get performance measures. A wrapper around the various criteria functions.
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%    path_ctrl - the output from SRR_solution_path
%    test_ctrl - struct with testX and testY, corrMat and true B
%    type - a string specifying the type of criteria
%    biascorrt - used for SRR_SCV
%    CVSplit_ctrl - used for SRR_SCV
%    type_ctrl - struct (optional)
%
% Outputs:
%    perf_ctrl - struct with optimal performance measures



if nargin < 4
    n = size(path_ctrl.X, 1);
    CVSplit_ctrl = CVSplit(n, 5); %%%5 fold CV by default
end
if nargin < 5
    type_ctrl = struct;
end

if strcmp(type, 'SF_PIC_fractional') || strcmp(type, 'AIC') || strcmp(type, 'BIC') || strcmp(type, 'EBIC') || strcmp(type, 'SCV_plugin') || strcmp(type, 'SCV_fractional') || strcmp(type, 'SCV_none')
    
    switch type
        case 'SF_PIC_fractional'
            opt_cri = SRR_SF_PIC_fractional(path_ctrl, type_ctrl);
        case 'AIC'
            opt_cri = SRR_AIC(path_ctrl);
        case 'BIC'
            opt_cri = SRR_BIC(path_ctrl);
        case 'EBIC'
            opt_cri = SRR_EBIC(path_ctrl);
        case 'SCV_plugin'
            opt_cri = SRR_SCV(path_ctrl, CVSplit_ctrl, 'plugin', type_ctrl);
        case 'SCV_fractional'
            opt_cri = SRR_SCV(path_ctrl, CVSplit_ctrl, 'fractional', type_ctrl);
        case 'SCV_none'
            opt_cri = SRR_SCV(path_ctrl, CVSplit_ctrl, 'none', type_ctrl);
    end
    
    opt_ctrl = SRR_find_opt_tuning(opt_cri, path_ctrl, type_ctrl);

    perf_ctrl = SRR_record_tuning(opt_ctrl, test_ctrl);
    
elseif strcmp(type, 'AIC_MA') || strcmp(type, 'BIC_MA') || strcmp(type, 'jackknife_MA')
    
    switch type
        case 'AIC_MA'
            opt_ctrl = SRR_IC_MA2(path_ctrl, 'AIC');
        case 'BIC_MA'
            opt_ctrl = SRR_IC_MA2(path_ctrl, 'BIC');
        case 'jackknife_MA'
            opt_ctrl = SRR_jackknife_MA2(path_ctrl, type_ctrl);
    end
        
    perf_ctrl = SRR_record_tuning(opt_ctrl, test_ctrl);
    
else
    
    error('type not implemented')
    
end


end

