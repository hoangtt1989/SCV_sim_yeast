% CV_time(curr_range) = CV_time_res(curr_range);
% PSCV_time(curr_range) = PSCV_time_res(curr_range);
% FSCV_time(curr_range) = FSCV_time_res(curr_range);
% path_time(curr_range) = path_time_res(curr_range);

for i = 1:length(method_names)
    
    test_name = strcat(method_names{i}, '_test');
    test_eval = eval(test_name);
    
    card_name = strcat(method_names{i}, '_card');
    card_eval = eval(card_name);
    
    rank_name = strcat(method_names{i}, '_rank');
    rank_eval = eval(rank_name);
    
    idx_name = strcat(method_names{i}, '_idx_cell');
    idx_eval = eval(idx_name);
   
    res.(method_names{i}).test(curr_range) = test_eval(curr_range);
    res.(method_names{i}).card(curr_range) = card_eval(curr_range);
    res.(method_names{i}).rank(curr_range) = rank_eval(curr_range);
    
    for j = curr_range
       res.(method_names{i}).idx{j} = idx_eval{j}; 
    end
    
end
