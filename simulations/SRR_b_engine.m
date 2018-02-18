function [ perf_tables ] = SRR_b_engine( size_ctrl, sim_ctrl, sel_ctrl, type_ctrl, conv_ctrl, verbose )
%SRR_b_engine - a wrapper around SRR_sim_engine that varies b
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%    size_ctrl - a struct with information about the problem dimensions and
%    correlations - also the signal strength(s) in b_vals
%    sim_ctrl - a struct with simulation information such as nsim
%    sel_ctrl - a struct with tuning method information
%    type_ctrl - a struct with problem type information
%    conv_ctrl - a struct with convergence information - card_upper_bnd and
%    rank_upper_bnd are important here
%
% Outputs:
%    test_errs - median prediction errors for each method
%    est_ranks - median of estimated ranks
%    est_cards - median of estimated cardinalities
%    MAs - mean of miss rates
%    FAs - mean of false alarm rates


b_len = length(size_ctrl.b_vals);
sel_len = length(sel_ctrl.methods);
% %%%test error
% test_errs = zeros(b_len, sel_len);
% %%%estimated rank
% est_ranks = zeros(b_len, sel_len);
% %%%estimated cardinality
% est_cards = zeros(b_len, sel_len);
% %%%miss rate
% MAs = zeros(b_len, sel_len);
% %%%false alarm rate
% FAs = zeros(b_len, sel_len);

row_names = sel_ctrl.methods;
col_names = {'MSE_median', 'MSE_std_error', 'J_hat_median', 'r_hat_median', 'M_mean', 'FA_mean'};

perf_tables = cell(b_len, 1);
% row_names = cell(b_len, 1);
for i = 1:b_len
   
    size_ctrl.b_strength = size_ctrl.b_vals(i);
    
    %%%test error
    [test_err_raw, sd_test_errs_raw, est_ranks, est_cards, MAs, FAs] = SRR_sim_engine(size_ctrl, sim_ctrl, sel_ctrl, type_ctrl, conv_ctrl, verbose);
    test_errs = test_err_raw / size_ctrl.b_strength^2;
    sd_test_errs = sd_test_errs_raw / size_ctrl.b_strength^2;
    
    tmp_mat(:, 1) = test_errs';
    tmp_mat(:, 2) = sd_test_errs';
    tmp_mat(:, 3) = est_cards';
    tmp_mat(:, 4) = est_ranks';
    tmp_mat(:, 5) = MAs';
    tmp_mat(:, 6) = FAs';
    
    perf_tables{i} = array2table(tmp_mat, 'VariableNames', col_names, 'RowNames', row_names);

%     row_names{i} = strcat('b', num2str(size_ctrl.b_vals(i), 2));
    
end



% test_errs = array2table(test_errs, 'VariableNames', sel_ctrl.methods, 'RowNames', row_names);

end

