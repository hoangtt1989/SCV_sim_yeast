
function [pSpa, spaErr, pNz, sucIden] = SparStat(resBetas, beta, T)
%%%%% Parameters: 
%%%%% resBetas - each column is an estimate; beta - true beta (column vector); T - sparsity matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    d = size(beta, 1);
    
    diffinds = nchoosek(1:d,2); 
    T = zeros(size(diffinds, 1), d);
    for tInd = 1:size(T, 1)
        T(tInd, diffinds(tInd, 2)) = -1; T(tInd, diffinds(tInd, 1)) = 1;
    end
    T = [eye(d); 1/sqrt(2) * T];
end

times = size(resBetas, 2);
if prod(size(resBetas)) == 0
    pSpa = -1; spaErr = -1;
    return;
end
tmpSpa = T*resBetas;
tmpSpa(abs(tmpSpa)<1e-3) = 0;
tmpSpa = (tmpSpa ~= 0);
tmpRefSpa = (T*beta ~= 0);
tmpRefInds = find(tmpRefSpa == 0);

pSpa = mean( tmpSpa(tmpRefInds,:) == repmat(tmpRefSpa(tmpRefInds), [1, times]), 1 ) ;
pNz =  mean( tmpSpa(tmpRefSpa,:) == repmat(tmpRefSpa(tmpRefSpa), [1, times]), 1 ) ;
spaErr = mean( tmpSpa ~= repmat(tmpRefSpa, [1, times]), 1 );

sucIden = mean(mean(tmpSpa(tmpRefSpa,:))==1);
