function [ perf_ctrl ] = SRR_jackknife_MA( path_ctrl, test_ctrl, type_ctrl )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if nargin < 3
    type_ctrl = struct;
end

if ~isfield(type_ctrl, 'pcamanner')
   type_ctrl.pcamanner = 0; 
end

num_tune = numel(path_ctrl.tune_params);

testX = test_ctrl.testX;
testY = test_ctrl.testY;

if isfield(test_ctrl, 'trainX') && isfield(test_ctrl, 'trainY')
    train_flag = 1;
    trainX = test_ctrl.trainX;
    trainY = test_ctrl.trainY;
    train_n = size(trainY, 1);
else
    train_flag = 0;
end

[test_n, m] = size(testY);

Es_test = zeros(test_n, m, num_tune);
preds_test = zeros(test_n, m, num_tune);
if train_flag
    Es_train = zeros(train_n, m, num_tune);
    preds_train = zeros(train_n, m, num_tune);
end

Bs = path_ctrl.Bs;
nz_patts = path_ctrl.nz_patts;
rank_Bs = path_ctrl.rank_Bs;

%%%%%filling in Z, E matrices
for i = 1:num_tune
   
    Zs_test = PatternExtr(testX, Bs(:, :, i), type_ctrl.pcamanner, nz_patts(i, :), rank_Bs(i));
    preds_test(:, :, i) = testX * Bs(:, :, i);
    resids_test = testY - preds_test(:, :, i);
    Q = qr(Zs_test' * Zs_test);
    inv_test = Q * Q';
    proj_test = Zs_test * inv_test * Zs_test';
    D_test = diag(proj_test).^(-1);
    Es_test(:, :, i) = bsxfun(@times, D_test, resids_test);
    
    if train_flag
        Zs_train = PatternExtr(trainX, Bs(:, :, i), type_ctrl.pcamanner, nz_patts(i, :), rank_Bs(i));
        preds_train(:, :, i) = trainX * Bs(:, :, i);
        resids_train = trainY - preds_train(:, :, i);
        Q = qr(Zs_train' * Zs_train);
        inv_train = Q * Q';
        proj_train = Zs_train * inv_train * Zs_train';
        D_train = diag(proj_train).^(-1);
        Es_train(:, :, i) = bsxfun(@times, D_train, resids_train);
    end
    
end
%%%%%

S_test = zeros(num_tune, num_tune);
if train_flag
    S_train = zeros(num_tune, num_tune);
end

%%%%%filling in S matrices
for i = 1:num_tune
    
    for j = 1:num_tune
        
        S_test(i, j) = sum(sum(Es_test(:, :, i) .* Es_test(:, :, j)));
        if train_flag
            S_train(i, j) = sum(sum(Es_train(:, :, i) .* Es_train(:, :, j)));
        end

    end
    
end
%%%%%

%%%%%solving the quadratic programming problem
S_test = S_test/(test_n * m);
wts_test = quadprog(S_test * 2, zeros(num_tune, 1), [], [], ones(1, num_tune), 1, zeros(num_tune, 1), ones(num_tune, 1));
% output_ctrl.w_test = w_test;
if train_flag
    S_train = S_test/(train_n * m);
    wts_train = quadprog(S_train * 2, zeros(num_tune, 1), [], [], ones(1, num_tune), 1, zeros(num_tune, 1), ones(num_tune, 1));
%     output_ctrl.w_train = w_train;
end
%%%%%


%%%%%evaluate errors
test_errs = zeros(test_n, m);
if train_flag
    train_errs = zeros(train_n, m);
end

for i = 1:num_tune
    
    test_errs = test_errs + wts_test(i) * preds_test(:, :, i);
    
    if train_flag
        train_errs = train_errs + wts_train(i) * preds_train(:, :, i);
    end
    
end

test_err = norm(test_ctrl.testY - test_errs, 'fro')^2 / (test_n * m);

perf_ctrl.test_errs = test_err;
perf_ctrl.est_rank = -1000;
perf_ctrl.est_card = -1000;
perf_ctrl.est_patt = -1000;

if train_flag
    train_errs = norm(test_ctrl.trainY - train_errs, 'fro')^2 / (train_n * m);
    perf_ctrl.train_errs = train_errs;
end
    
end

