function [Best, ks, eigvals, B0Opt, B1Opt] = Func_RRR(X, Y, ks)
%Reduced rank regression 

tmpUniEst = X' * Y;
tmpSigma = X' * X;

%B_unpen = tmpSigma\tmpUniEst; % should be generalized inverse
% % pinv(tmpSigma) sometimes chokes SVD because of the scale
% scale_tmp = 100; % norm(X, 'fro')^2; % 100; %10*max(max(abs(tmpSigma)));
% tmpSigmaGInv = pinv(tmpSigma/scale_tmp)/scale_tmp; 
% Instead, use the following way:    
tmpXGInv = pinv(X);
tmpSigmaGInv = tmpXGInv * tmpXGInv';

B_unpen = tmpSigmaGInv * tmpUniEst;

eigTarget = tmpUniEst' * B_unpen;
[eigV0, eigD0] = eig((eigTarget+eigTarget')/2);
[eigDdiag, IX] = sort(diag(eigD0), 'descend');
eigD = eigD0(IX, IX); 
eigV = eigV0(:, IX);

[m, n] = size(Y); d = size(X, 2);
r = min(n, sum(svd(X)>1e-6)); % r = min(n, d); 
% ks may be changed!
ks(ks >= r) = r;
    
B0 = B_unpen * eigV;
B1 = eigV';

if numel(ks) == 1
    k = ks;
    
    Best = B0(:,1:k) * B1(1:k, :);
    eigvals = eigDdiag;
    B0Opt = B0(:,1:k);
    B1Opt = B1(1:k, :);
else
    Best = zeros(d, n, numel(ks));
    for kInd=1:numel(ks)
        k = ks(kInd);
        Best(:, :, kInd) = B0(:,1:k) * B1(1:k, :);
    end

    eigvals = eigDdiag;
    B0Opt = B0; 
    B1Opt = B1;   
end

