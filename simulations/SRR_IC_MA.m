function [ perf_ctrl ] = SRR_IC_MA( path_ctrl, test_ctrl, type )
%SRR_IC_MA - use model averaging on AIC or BIC to get a final test error or
%training error
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%    opt_ctrl - the output from SRR_find_opt_tuning
%    path_ctrl - the output from SRR_solution_path
%    test_ctrl - struct with testX and testY, corrMat and true B
%
% Outputs:
%    perf_ctrl - a struct with performance information

if nargin < 3
    type = 'AIC';
end

switch type
    case 'AIC'
        opt_cri = SRR_AIC(path_ctrl);
    case 'BIC'
        opt_cri = SRR_BIC(path_ctrl);
end

m = size(path_ctrl.Y, 2);

testX = test_ctrl.testX;
test_n = size(testX, 1);
if isfield(test_ctrl, 'trainX')
    train_flag = 1;
    trainX = test_ctrl.trainX;
    train_n = size(trainX, 1);
else
    train_flag = 0;
end

num_tune = numel(path_ctrl.tune_params);

wts = exp(-.5 * opt_cri) ./ sum(exp(-.5 * opt_cri));

preds = zeros(test_n, m);
if train_flag
    trains = zeros(train_n, m);
end

for i = 1:num_tune
    
    preds = preds + wts(i) * testX * path_ctrl.Bs(:, :, i);
    
    if train_flag
        trains = trains + wts(i) * trainX * path_ctrl.Bs(:, :, i);
    end
    
end

preds_err = norm(test_ctrl.testY - preds, 'fro')^2 / (test_n * m);

perf_ctrl.test_errs = preds_err;
perf_ctrl.est_rank = -1000;
perf_ctrl.est_card = -1000;
perf_ctrl.est_patt = -1000;

if train_flag
    trains_err = norm(test_ctrl.trainY - trains, 'fro')^2 / (train_n * m);
    perf_ctrl.train_errs = trains_err;
end

end

