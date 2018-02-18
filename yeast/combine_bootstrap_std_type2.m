clear
close all

TF21_pos = csvread('TF21_pos.csv');

res = struct;
res.AIC = struct;
res.BIC = struct;
res.CV = struct;
res.EBIC = struct;
res.FSCV = struct;
res.FSCV_split = struct;
res.PIC = struct;
res.PSCV = struct;
res.PSCV_split = struct;

method_names = {'AIC', 'BIC', 'EBIC', 'FSCV', 'PIC', 'PSCV', 'CV'};

load('bootstrap_std_type2_1_40.mat')
run combine_repeated.m

load('bootstrap_std_type2_41_80.mat')
run combine_repeated.m

load('bootstrap_std_type2_81_120.mat')
run combine_repeated.m

load('bootstrap_std_type2_121_200.mat')
run combine_repeated.m

%%%%%get medians
for i = 1:length(method_names)
   
    res.(method_names{i}).test_median = median(res.(method_names{i}).test);
    res.(method_names{i}).test_mean = mean(res.(method_names{i}).test);
    res.(method_names{i}).test_std = std(res.(method_names{i}).test);
    res.(method_names{i}).card_median = median(res.(method_names{i}).card);
    res.(method_names{i}).card_std = std(res.(method_names{i}).card);
    res.(method_names{i}).rank_median = median(res.(method_names{i}).rank);
    res.(method_names{i}).rank_std = std(res.(method_names{i}).rank);
    
    disp([method_names{i}, ' median of test errors ', num2str(res.(method_names{i}).test_median)])
    disp([method_names{i}, ' mean of test errors ', num2str(res.(method_names{i}).test_mean)])
    disp([method_names{i}, ' std of test errors ', num2str(res.(method_names{i}).test_std)])
    disp([method_names{i}, ' median of cardinalities ', num2str(res.(method_names{i}).card_median)])
    disp([method_names{i}, ' median of ranks ', num2str(res.(method_names{i}).rank_median)])
    
    res.(method_names{i}).idx_all = [];
    res.(method_names{i}).idx_mat = zeros(106, 106);
    
    %%%%%get a vector of all the selected predictors
    for j = 1:length(res.(method_names{i}).test)
       res.(method_names{i}).idx_all = [res.(method_names{i}).idx_all res.(method_names{i}).idx{j}];
       res.(method_names{i}).idx_mat(j, res.(method_names{i}).idx{j}) = 1;
    end
    res.(method_names{i}).idx_pairs = res.(method_names{i}).idx_mat' * res.(method_names{i}).idx_mat;
    
%     [~, top_pairs] = sort(diag(res.(method_names{i}).idx_pairs), 1, 'descend');
    
    
    res.(method_names{i}).idx_freq = tabulate(res.(method_names{i}).idx_all);
    res.(method_names{i}).idx_freq = res.(method_names{i}).idx_freq(:, 1:2);
    
end



%%%%output data for plotting
rank_tab = table(res.PSCV.rank', res.FSCV.rank', res.CV.rank', 'VariableNames', {'PSCV', 'FSCV', 'CV'});
writetable(rank_tab, 'bootstrap_std_type2_ranks.csv')

card_tab = table(res.PSCV.card', res.FSCV.card', res.CV.card', 'VariableNames', {'PSCV', 'FSCV', 'CV'});
writetable(card_tab, 'bootstrap_std_type2_cards.csv')

filename = 'bootstrap_std_type2_all_freqs.xls';
var_names = {'idx', 'freq'};
write_names = {'PSCV', 'FSCV', 'CV'};
for i = 1:length(write_names)
    writetable(array2table(res.(write_names{i}).idx_freq, 'VariableNames', var_names), filename, 'FileType', 'spreadsheet', 'Sheet', write_names{i}, 'Range', 'A1')
end

res.PSCV.TF21_freq = res.PSCV.idx_freq(TF21_pos, :);
res.FSCV.TF21_freq = res.FSCV.idx_freq(TF21_pos, :);
res.CV.TF21_freq = res.CV.idx_freq(TF21_pos, :);

filename = 'bootstrap_std_type2_TF21_freqs.xls';
var_names = {'idx', 'freq'};
write_names = {'PSCV', 'FSCV', 'CV'};
for i = 1:length(write_names)
    writetable(array2table(res.(write_names{i}).TF21_freq, 'VariableNames', var_names), filename, 'FileType', 'spreadsheet', 'Sheet', write_names{i}, 'Range', 'A1')
end
%%%%

