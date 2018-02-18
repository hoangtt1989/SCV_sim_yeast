function [ opt_ctrl ] = SRR_jackknife_MA2( path_ctrl, type_ctrl )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    type_ctrl = struct;
end

if ~isfield(type_ctrl, 'pcamanner')
   type_ctrl.pcamanner = 0; 
end

num_tune = numel(path_ctrl.tune_params);

X = path_ctrl.X;
Y = path_ctrl.Y;
[n, m] = size(Y);
p = size(X, 2);

Es = zeros(n, m, num_tune);
% preds = zeros(n, m, num_tune);

Bs = path_ctrl.Bs;
nz_patts = path_ctrl.nz_patts;
rank_Bs = path_ctrl.rank_Bs;

%%%%%filling in Z, E matrices
for i = 1:num_tune
   
    Zs = PatternExtr(X, Bs(:, :, i), type_ctrl.pcamanner, nz_patts(i, :), rank_Bs(i));
%     preds = X * Bs(:, :, i);
    resids = Y - X * Bs(:, :, i);
    Q = qr(Zs' * Zs);
    inv_inner = Q * Q';
    proj_Zs = Zs * inv_inner * Zs';
    D_test = diag(proj_Zs).^(-1);
    Es(:, :, i) = bsxfun(@times, D_test, resids);
    
end
%%%%%

S_mat = zeros(num_tune, num_tune);

%%%%%filling in S matrices
for i = 1:num_tune
    
    for j = 1:num_tune
        
        S_mat(i, j) = sum(sum(Es(:, :, i) .* Es(:, :, j)));

    end
    
end
%%%%%

%%%%%solving the quadratic programming problem
S_mat = S_mat/(n * m);
wts = quadprog(S_mat * 2, zeros(num_tune, 1), [], [], ones(1, num_tune), 1, zeros(num_tune, 1), ones(num_tune, 1));
%%%%%


%%%%%average B matrix
BOpt = zeros(p, m);
for i = 1:num_tune
   
    BOpt = BOpt + wts(i) * Bs(:, :, i);
    
end
%%%%%


opt_ctrl.BOpt = BOpt;
opt_ctrl.card_opt = -1000;
opt_ctrl.rank_opt = -1000;
opt_ctrl.opt_ind = -1000;
opt_ctrl.opt_cri = -1000;
opt_ctrl.nz_patt = -1000;
    
end

