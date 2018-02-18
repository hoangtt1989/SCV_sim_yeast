function [X, Y, B_true, valX, valY, tstX, tstY, E, corrMat, rho] = ...
    Func_GeneSRRRData(n, m, d, s, r, sigma, bStrength, dataName, valN, tstN, B_raw, q)

if ~exist('B_raw', 'var') || isempty(B_raw) 
    B_raw = randn(s, r) * randn(r, m);
    B_raw = [B_raw; zeros(d-s, m)];
end
B_true = bStrength * B_raw; % randn(d, r) * randn(r, n);
totN = n + valN + tstN;
if strcmp(dataName(1:2), 'my')   
    if ~exist('q', 'var')
        X = randn(totN, d);
    else
        X = randn(totN, q) * randn(q, d);
    end
    if strcmp(dataName, 'mynorm-iid')
        rho=0;  
        corrMat = rho .^ abs([1:d]' * ones(1, d) - ones(d, 1) * [1:d]);
        X = X * chol(corrMat);
    elseif strcmp(dataName, 'mynorm-highcorr')
        rho=.9;  
        corrMat = rho .^ abs([1:d]' * ones(1, d) - ones(d, 1) * [1:d]);
        X = X * chol(corrMat);
    elseif strcmp(dataName, 'mynorm-medcorr')
        rho=.5;  
        corrMat = rho .^ abs([1:d]' * ones(1, d) - ones(d, 1) * [1:d]);
        X = X * chol(corrMat);
    elseif strcmp(dataName, 'mynorm-mildcorr')
        rho=.1;  
        corrMat = rho .^ abs([1:d]' * ones(1, d) - ones(d, 1) * [1:d]);
        X = X * chol(corrMat);        
    end        
    Y = X * B_true + sigma * randn(totN, m);
    
    tstX = X((n+valN+1):end, :); valX = X((n+1):(n+valN), :); X = X(1:n, :); 
    tstY = Y((n+valN+1):end, :); valY = Y((n+1):(n+valN), :); Y = Y(1:n, :); 

    centerX = 0; scaleX = 0;
    E = Y - X * B_true;
end


end
