function X_th = MultiTh_mat(X, lambda, thresholdingWay, Nu)
% Note: X can be a vector as well! 

if ~exist('Nu', 'var') || isempty(Nu)
    Nu = 0;
end
if ~exist('thresholdingWay', 'var')
    thresholdingWay = 'hard';
end

X_th = zeros(size(X));
grpnorms = sqrt(sum(X.^2, 2));

% disp(grpnorms');  

switch lower(thresholdingWay)
    case {'soft'}    
        ttInds = find(grpnorms>lambda);
        X_th(ttInds, :) = diag((grpnorms(ttInds) - lambda)./grpnorms(ttInds)) * X(ttInds, :); 
 
    case {'hard'}
        ttInds = find(grpnorms>lambda);
        X_th(ttInds, :) = X(ttInds, :);  
        
    case {'prop'}
        quan = 1 - lambda / numel(grpnorms);
        quan(quan < 0) =0; 
        thval = quantile(grpnorms, quan);                     
        thval(quan==0) = 0;

        ttInds = find(grpnorms > thval);        
        X_th(ttInds, :) = X(ttInds, :);  % X(ttInds, :) / (1+1e-5)

	case {'hybrid'}
        ttInds = find(grpnorms>lambda);
        X_th(ttInds, :) = X(ttInds, :) ./(1+Nu);  
        
    case {'hybrid-prop'}
        quan = 1 - lambda / numel(grpnorms);
        quan(quan < 0) =0; 
        thval = quantile(grpnorms, quan);                     
        thval(quan==0) = 0;

        ttInds = find(grpnorms > thval);        
        X_th(ttInds, :) = X(ttInds, :) ./ (1+Nu);
                    
end




