function [ perf_ctrl ] = SRR_record_tuning( opt_ctrl, test_ctrl )
%SRR_record_tuning - after running getting the optimal model, record some
%performance measures
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%    opt_ctrl - the output from SRR_find_opt_tuning
%    path_ctrl - the output from SRR_solution_path
%    test_ctrl - struct with testX and testY, corrMat and true B
%
% Outputs:
%    record_ctrl - struct with measures


%%%initialize from opt_ctrl struct
BOpt = opt_ctrl.BOpt;
%%%

% Record the statistics useful for multiple run evaluation of the algorithms
%%%initialize from test_ctrl struct
if isfield(test_ctrl, 'testX') && isfield(test_ctrl, 'testY')
    tstY = test_ctrl.testY;
    tstX = test_ctrl.testX;
    %%%
    modelErrs = norm(tstY-tstX*BOpt,'fro')^2/numel(tstY)/1;
    perf_ctrl.test_errs = modelErrs;
end
if isfield(test_ctrl, 'trainX') && isfield(test_ctrl, 'trainY')
    trn_errs = norm(test_ctrl.trainY - test_ctrl.trainX * BOpt, 'fro')^2 / numel(test_ctrl.trainY) / 1;
    perf_ctrl.train_errs = trn_errs;
end

if isfield(test_ctrl, 'B_true')
    B_true = test_ctrl.B_true;
    estErrs = norm(B_true-BOpt,'fro')^2/numel(B_true)/1;
    perf_ctrl.est_errs = estErrs;
    if isfield(test_ctrl, 'corrMat')
        corrMat = test_ctrl.corrMat;
        new_tmp = (BOpt - B_true);
        new_modelErrs = trace(new_tmp'*corrMat*new_tmp);
        perf_ctrl.test_errs_corr = new_modelErrs;
    end
end

%%%%these are taken from opt_ctrl
perf_ctrl.est_rank = opt_ctrl.rank_opt;
perf_ctrl.est_card = opt_ctrl.card_opt;
perf_ctrl.est_patt = opt_ctrl.nz_patt;
perf_ctrl.BOpt = BOpt;
%%%%

end

