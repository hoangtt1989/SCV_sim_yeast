clear
close all

rng('default')


resp = csvread('Y_mat.csv');
pred = csvread('X_mat.csv');


n = size(pred, 1);
m = size(resp, 2);

rank_up = 15;
card_up = 106;


nrep = 200;

pred_train_cell = cell(nrep, 1);
pred_train_orig = cell(nrep, 1);
pred_test_cell = cell(nrep, 1);
pred_test_orig = cell(nrep, 1);

resp_train_cell = cell(nrep, 1);
resp_train_orig = cell(nrep, 1);
resp_test_cell = cell(nrep, 1);
resp_test_orig = cell(nrep, 1);

CVSplit_cell = cell(nrep, 1);

cov_eps = 1e-5;

for i = 1:nrep
    
    
    CVSplit_cell{i} = CVSplit(n, 5);
    
    train_idx = randsample(n, n, true);
%     test_idx = setdiff(1:n, train_idx);
    
    %%%%center and decorrelate response
    resp_train = resp(train_idx, :);
    resp_train_orig{i} = resp_train;

    resp_train_mn = mean(resp_train);
    resp_train_std = std(resp_train);
    
    resp_train_cent = bsxfun(@minus, resp_train, resp_train_mn);
%     icov = sqrtm(inv(resp_train_cent' * resp_train_cent / (length(train_idx)-1) + cov_eps * eye(m)));
%     resp_train_dec = resp_train_cent * icov;
    resp_train_dec = bsxfun(@rdivide, resp_train_cent, resp_train_std);
    resp_train_cell{i} = resp_train_dec;
    %%%%
    
    %%%%center and scale predictors
    pred_train = pred(train_idx, :);
    
    pred_train_orig{i} = pred_train;
    
    pred_train_mn = mean(pred_train);
    
    pred_train_sd = std(pred_train);
    pred_train_sd(pred_train_sd == 0) = 1;
    
    pred_train_std = bsxfun(@minus, pred_train, pred_train_mn);
    pred_train_std = bsxfun(@rdivide, pred_train_std, pred_train_sd);
    pred_train_cell{i} = pred_train_std;
    %%%%
    
end

path_time_res = zeros(nrep, 1);
PSCV_time_res = zeros(nrep, 1);
CV_time_res = zeros(nrep, 1);

PIC_test = zeros(nrep, 1);
PIC_card = zeros(nrep, 1);
PIC_rank = zeros(nrep, 1);
PIC_idx_cell = cell(nrep, 1);

EBIC_test = zeros(nrep, 1);
EBIC_card = zeros(nrep, 1);
EBIC_rank = zeros(nrep, 1);
EBIC_idx_cell = cell(nrep, 1);

BIC_test = zeros(nrep, 1);
BIC_card = zeros(nrep, 1);
BIC_rank = zeros(nrep, 1);
BIC_idx_cell = cell(nrep, 1);

AIC_test = zeros(nrep, 1);
AIC_card = zeros(nrep, 1);
AIC_rank = zeros(nrep, 1);
AIC_idx_cell = cell(nrep, 1);

PSCV_test = zeros(nrep, 1);
PSCV_card = zeros(nrep, 1);
PSCV_rank = zeros(nrep, 1);
PSCV_idx_cell = cell(nrep, 1);

FSCV_test = zeros(nrep, 1);
FSCV_card = zeros(nrep, 1);
FSCV_rank = zeros(nrep, 1);
FSCV_idx_cell = cell(nrep, 1);

SCV_test = zeros(nrep, 1);
SCV_card = zeros(nrep, 1);
SCV_rank = zeros(nrep, 1);
SCV_idx_cell = cell(nrep, 1);

CV_test = zeros(nrep, 1);
CV_card = zeros(nrep, 1);
CV_rank = zeros(nrep, 1);
CV_idx_cell = cell(nrep, 1);

CV_split_test = zeros(nrep, 1);
CV_split_card = zeros(nrep, 1);
CV_split_rank = zeros(nrep, 1);
CV_split_idx_cell = cell(nrep, 1);

FSCV_split_test = zeros(nrep, 1);
FSCV_split_card = zeros(nrep, 1);
FSCV_split_rank = zeros(nrep, 1);
FSCV_split_idx_cell = cell(nrep, 1);

PSCV_split_test = zeros(nrep, 1);
PSCV_split_card = zeros(nrep, 1);
PSCV_split_rank = zeros(nrep, 1);
PSCV_split_idx_cell = cell(nrep, 1);

curr_range = 1:200;
curr_range

tic
parpool
parfor rep_it = curr_range
    
    disp(rep_it)
    
    input_ctrl = struct;
    input_ctrl.X = pred_train_cell{rep_it};
    input_ctrl.Y = resp_train_cell{rep_it};
    
    test_ctrl = struct;
    test_ctrl.testX = pred_train_cell{rep_it};
    test_ctrl.testY = resp_train_cell{rep_it};
    test_ctrl.trainX = input_ctrl.X;
    test_ctrl.trainY = input_ctrl.Y;
    
    [n, ~] = size(input_ctrl.Y);
    
    path_ctrl = struct;
    conv_ctrl = struct;
    conv_ctrl.card_upper_bnd = card_up;
    conv_ctrl.rank_upper_bnd = rank_up;
    conv_ctrl.zeroRem = 0;
%     type_ctrl = struct;
    
    tic
    [path_ctrl, conv_ctrl, type_ctrl] = SRR_solution_path(input_ctrl, conv_ctrl);
    path_time = toc;
    path_time_res(rep_it) = path_time;
    
    
    CVSplit_ctrl = CVSplit_cell{rep_it};
    
    %%%%PIC
    perf_ctrl_PIC = SRR_tune_perf(path_ctrl, test_ctrl, 'SF_PIC_fractional');
    PIC_idx = find(perf_ctrl_PIC.est_patt ~= 0);
%     PIC_res{rep_it} = perf_ctrl_PIC;
    PIC_test(rep_it) = perf_ctrl_PIC.test_errs;
    PIC_rank(rep_it) = perf_ctrl_PIC.est_rank;
    PIC_card(rep_it) = perf_ctrl_PIC.est_card;
    PIC_idx_cell{rep_it} = PIC_idx;
    
    %%%%EBIC
    perf_ctrl_EBIC = SRR_tune_perf(path_ctrl, test_ctrl, 'EBIC');
    EBIC_idx = find(perf_ctrl_EBIC.est_patt ~= 0);
%     EBIC_res{rep_it} = perf_ctrl_EBIC;
    EBIC_test(rep_it) = perf_ctrl_EBIC.test_errs;
    EBIC_rank(rep_it) = perf_ctrl_EBIC.est_rank;
    EBIC_card(rep_it) = perf_ctrl_EBIC.est_card;
    EBIC_idx_cell{rep_it} = EBIC_idx;
    
    %%%%BIC
    perf_ctrl_BIC = SRR_tune_perf(path_ctrl, test_ctrl, 'BIC');
    BIC_idx = find(perf_ctrl_BIC.est_patt ~= 0);
%     BIC_res{rep_it} = perf_ctrl_BIC;
    BIC_test(rep_it) = perf_ctrl_BIC.test_errs;
    BIC_rank(rep_it) = perf_ctrl_BIC.est_rank;
    BIC_card(rep_it) = perf_ctrl_BIC.est_card;
    BIC_idx_cell{rep_it} = BIC_idx;
    
    %%%%AIC
    perf_ctrl_AIC = SRR_tune_perf(path_ctrl, test_ctrl, 'AIC');
    AIC_idx = find(perf_ctrl_AIC.est_patt ~= 0);
%     AIC_res{rep_it} = perf_ctrl_AIC;
    AIC_test(rep_it) = perf_ctrl_AIC.test_errs;
    AIC_rank(rep_it) = perf_ctrl_AIC.est_rank;
    AIC_card(rep_it) = perf_ctrl_AIC.est_card;
    AIC_idx_cell{rep_it} = AIC_idx;
    
    %%%%SCV plugin
%     tic
    perf_ctrl_PSCV = SRR_tune_perf(path_ctrl, test_ctrl, 'SCV_plugin', CVSplit_ctrl);
%     PSCV_time = toc;
    PSCV_idx = find(perf_ctrl_PSCV.est_patt ~= 0);
%     PSCV_res{rep_it} = perf_ctrl_PSCV;
    PSCV_test(rep_it) = perf_ctrl_PSCV.test_errs;
    PSCV_rank(rep_it) = perf_ctrl_PSCV.est_rank;
    PSCV_card(rep_it) = perf_ctrl_PSCV.est_card;
    PSCV_idx_cell{rep_it} = PSCV_idx;
    
    %%%%SCV fractional
%     tic
    perf_ctrl_FSCV = SRR_tune_perf(path_ctrl, test_ctrl, 'SCV_fractional', CVSplit_ctrl);
%     FSCV_time = toc;
    %%%look at pattern
    FSCV_idx = find(perf_ctrl_FSCV.est_patt ~= 0);
%     FSCV_res{rep_it} = perf_ctrl_FSCV;
    FSCV_test(rep_it) = perf_ctrl_FSCV.test_errs;
    FSCV_rank(rep_it) = perf_ctrl_FSCV.est_rank;
    FSCV_card(rep_it) = perf_ctrl_FSCV.est_card;
    FSCV_idx_cell{rep_it} = FSCV_idx;
    %%%%
    
    %%%%SCV no correction
    perf_ctrl_SCV = SRR_tune_perf(path_ctrl, test_ctrl, 'SCV_none', CVSplit_ctrl);
    SCV_idx = find(perf_ctrl_SCV.est_patt ~= 0);
%     SCV_res{rep_it} = perf_ctrl_SCV;
    SCV_test(rep_it) = perf_ctrl_SCV.test_errs;
    SCV_rank(rep_it) = perf_ctrl_SCV.est_rank;
    SCV_card(rep_it) = perf_ctrl_SCV.est_card;
    SCV_idx_cell{rep_it} = SCV_idx;
    %%%%
    
    %%%%CV regular
    opt_CV = SRR_CV(input_ctrl, CVSplit_ctrl, conv_ctrl, type_ctrl);
    perf_ctrl_CV = SRR_record_tuning(opt_CV, test_ctrl);
    CV_idx = find(perf_ctrl_CV.est_patt ~= 0);
    CV_test(rep_it) = perf_ctrl_CV.test_errs;
    CV_rank(rep_it) = perf_ctrl_CV.est_rank;
    CV_card(rep_it) = perf_ctrl_CV.est_card;
    CV_idx_cell{rep_it} = CV_idx;
    %%%%
    
    
    type_ctrl.white = 0;
    input_ctrl = struct;
    input_ctrl.oriX = pred_train_orig{rep_it};
    input_ctrl.oriY = resp_train_orig{rep_it};
    input_ctrl.transX = pred_train_cell{rep_it};
    input_ctrl.transY = resp_train_cell{rep_it};

    
end
toc

delete(gcp('nocreate'))

filename = ['bootstrap_std_type2' '_' num2str(curr_range(1)) '_' num2str(curr_range(end)) '.mat'];
clear pred_train_cell pred_test_cell pred_train_orig pred_test_orig resp_train_cell resp_test_cell resp_train_orig resp_test_orig
save(filename)
