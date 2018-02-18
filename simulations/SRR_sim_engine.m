function [ test_errs, sd_test_errs, est_ranks, est_cards, MAs, FAs ] = SRR_sim_engine( size_ctrl, sim_ctrl, sel_ctrl, type_ctrl, conv_ctrl, verbose  )
%SRR_sim_engine - for some number of simulations, get the test errors using
%different tuning methods
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%    size_ctrl - a struct with information about the problem dimensions and
%    correlations
%    sim_ctrl - a struct with simulation information such as nsim
%    sel_ctrl - a struct with tuning method information
%    type_ctrl - a struct with problem type information
%    conv_ctrl - a struct with convergence information - card_upper_bnd and
%    rank_upper_bnd are important here
%
% Outputs:
%    test_errs - median prediction errors for each method

rng('default')

if nargin < 4
    type_ctrl = struct;
end
if nargin < 5
    conv_ctrl = struct;
end
if nargin < 6
    verbose = 1;
end

if ~isfield(size_ctrl, 'n_true')
    error('Must supply size_ctrl.n_true')
end
if ~isfield(size_ctrl, 'p_true')
    error('Must supply size_ctrl.p_true')
end
if ~isfield(size_ctrl, 'J_true')
    error('Must supply size_ctrl.J_true')
end
if ~isfield(size_ctrl, 'm_true')
    error('Must supply size_ctrl.m_true')
end
if ~isfield(size_ctrl, 'r_true')
    error('Must supply size_ctrl.r_true')
end
if ~isfield(size_ctrl, 'dataName')
    if verbose
        disp('size_ctrl.dataName not supplied using default of mynorm-medcorr')
    end
    size_ctrl.dataName = 'mynorm-medcorr';
end
if ~isfield(size_ctrl, 'b_strength')
    error('Must supply size_ctrl.b_strength')
end
if ~isfield(size_ctrl, 'sigma')
    if verbose
        disp('size_ctrl.sigma not supplied using default of 1')
    end
    size_ctrl.sigma = 1;
end
if ~isfield(size_ctrl, 'valN')
    if verbose
        disp('size_ctrl.valN not supplied using default of 1e+4')
    end
    size_ctrl.valN = 1e+4;
end
if ~isfield(size_ctrl, 'tstN')
    if verbose
        disp('size_ctrl.tstN not supplised using default of 1e+4')
    end
    size_ctrl.tstN = 1e+4;
end
if ~isfield(sim_ctrl, 'nsim')
    if verbose
        disp('sim_ctrl.nsim not supplied using default of 100')
    end
    sim_ctrl.nsim = 100;
end
if ~isfield(sim_ctrl, 'genBOnce')
    if verbose
        disp('sim_ctrl.genBOnce not supplied using default of 1')
    end
    sim_ctrl.genBOnce = 1;
end
    

%%%initialize from size_ctrl struct
n_true = size_ctrl.n_true;
m_true = size_ctrl.m_true;
p_true = size_ctrl.p_true;
J_true = size_ctrl.J_true;
r_true = size_ctrl.r_true;
sigma = size_ctrl.sigma;
b_strength = size_ctrl.b_strength;
dataName = size_ctrl.dataName;
valN = size_ctrl.valN;
tstN = size_ctrl.tstN;
%%%
%%%initialize from sim_ctrl struct
nsim = sim_ctrl.nsim;
genBOnce = sim_ctrl.genBOnce;
%%%

%%%for model selection/averaging
% tune_methods = sel_ctrl.methods;
sel_len = length(sel_ctrl.methods);
test_errs = zeros(nsim, sel_len);
est_ranks = zeros(nsim, sel_len);
est_cards = zeros(nsim, sel_len);
MAs = zeros(nsim, sel_len);
FAs = zeros(nsim, sel_len);
%%%

% CVSplit_ctrls = [];
    
for timeInd = 1:nsim

    %%%%%generate data
    if genBOnce == 1 && timeInd == 1
        B_raw = randn(J_true, r_true) * randn(r_true, m_true);
        B_raw = [B_raw; zeros(p_true-J_true, m_true)];
    end
    
    [X, Y, B_true, valX, valY, tstX, tstY, E, corrMat, rho] = ...
        Func_GeneSRRRData(n_true, m_true, p_true, J_true, r_true, sigma, b_strength, dataName, valN, tstN, B_raw);
    %%%%%
    
    q = sum(svd(X) > 1e-4);
    
    %%%%CV splits
    CVSplit_ctrl5 = CVSplit(n_true, 5);
    CVSplit_ctrl2 = CVSplit(n_true, 2);
    CVSplit_ctrl10 = CVSplit(n_true, 10);
    %%%%
    
    %%%%initializing structs
    test_ctrl.testX = tstX;
    test_ctrl.testY = tstY;
    test_ctrl.corrMat = corrMat;
    test_ctrl.B_true = B_true;
    
    input_ctrl.X = X;
    input_ctrl.Y = Y;
    %%%%
    
    %%%%get true sparsity pattern
    true_nz = find(sqrt(sum(B_true'.^2)) > 1e-4);
    %%%
    
    %%%%solution path
    [path_ctrl, conv_ctrl, type_ctrl] = SRR_solution_path(input_ctrl, conv_ctrl, type_ctrl);
    path_ctrl.q = q;
    %%%%
    
    
    for tunes = 1:sel_len
        
        if strcmp(sel_ctrl.methods{tunes}, 'SCV_fractional') || strcmp(sel_ctrl.methods{tunes}, 'SCV_plugin') || strcmp(sel_ctrl.methods{tunes}, 'SCV_none5')
            CVSplit_ctrl = CVSplit_ctrl5;
        elseif strcmp(sel_ctrl.methods{tunes}, 'SCV_none2')
            CVSplit_ctrl = CVSplit_ctrl2;
        elseif strcmp(sel_ctrl.methods{tunes}, 'SCV_none10')
            CVSplit_ctrl = CVSplit_ctrl10;
        else
            CVSplit_ctrl = CVSplit_ctrl5;
        end
        
        perf_ctrl = struct;
        
        perf_ctrl = SRR_tune_perf(path_ctrl, test_ctrl, sel_ctrl.methods{tunes}, CVSplit_ctrl, type_ctrl);
        
        est_nz = find(perf_ctrl.est_patt ~= 0);
        
        %%%test error
        test_errs(timeInd, tunes) = perf_ctrl.test_errs_corr;
        %%%estimated rank
        est_ranks(timeInd, tunes) = perf_ctrl.est_rank;
        %%%estimated cardinality
        est_cards(timeInd, tunes) = perf_ctrl.est_card;
        %%%miss rate
        MAs(timeInd, tunes) = numel(setdiff(true_nz, est_nz)) / numel(true_nz);
        %%%false alarm rate
        FAs(timeInd, tunes) = numel(setdiff(est_nz, true_nz)) / (p_true - J_true);
%         FAs(timeInd, tunes) = numel(setdiff(est_nz, true_nz)) / numel(est_nz);
%         FAs(timeInd, tunes) = numel(setdiff(est_nz, true_nz)) / numel(true_nz);
        
    end
    
    
end

sd_test_errs = std(test_errs / m_true);
test_errs = median(test_errs) / m_true;
est_ranks = median(est_ranks);
est_cards = median(est_cards);
MAs = mean(MAs);
FAs = mean(FAs);

end

