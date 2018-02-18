function [ optCri, type, type_ctrl ] = SRR_SF_PIC_fractional( path_ctrl, type_ctrl )

%SF_PIC_fractional - for a given solution path compute the values of
%scale-free PIC fractional form. path_ctrl is the output_ctrl from
%SRR_solution_path
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%    path_ctrl - a struct with solution path information
%    type_ctrl - a struct with problem type information
%
% Outputs:
%    opt_cri
%    type_ctrl


if nargin < 2
    type_ctrl = struct;
end

if ~isfield(path_ctrl, 'q')
    q = sum(svd(path_ctrl.X) > 1e-4);
else
    q = path_ctrl.q;
end
if ~isfield(type_ctrl, 'pcamanner')
    type_ctrl.pcamanner = 0;
end


%%%initialize from path_ctrl struct
X = path_ctrl.X;
Y = path_ctrl.Y;
% q = path_ctrl.q;
%%%
%%%initialize from path_ctrl struct
rank_Bs = path_ctrl.rank_Bs;
card_Bs = path_ctrl.card_Bs;
trn_errs = path_ctrl.trn_errs;
%%%


% SF-PIC: fractional form
A1 = 2.0;
A2 = 1.8;
A3 = 2.0;

if type_ctrl.pcamanner == 0 % For general problems with iid E:
    pen = A1 * (min(q, card_Bs)).*rank_Bs + A3*(size(Y, 2) - 1 * rank_Bs).*rank_Bs + ...
        A2 * card_Bs.*(1+ log(size(X, 2)./card_Bs));
    pen = 1* pen;
    optCri = trn_errs ./ (numel(Y) - pen);
    optCri(optCri<0) = inf;
elseif type_ctrl.pcamanner == 1 % For PCA-transpose problems
    error('Implement it NOW!')
elseif type_ctrl.pcamanner == 2 % For PCA-selregression problems with noise E V^T
    pen = A1 * (min(q, card_Bs)).*rank_Bs + A3*(size(Y, 2) - 2 * rank_Bs).*rank_Bs + ...
        A2 * card_Bs.*(1+ log(size(X, 2)./card_Bs));
    %         pen = 1* sigma^2 * pen;
    pen = 1 * pen;
    %         optCri = trn_errs/ (numel(Y) - pen); optCri(optCri<0) = inf;
    optCri = trn_errs ./ ( size(Y,1).*(size(Y,2) - rank_Bs) - pen ); optCri(optCri<0) = inf;
end

type = 'SF_PIC_fractional';



end

