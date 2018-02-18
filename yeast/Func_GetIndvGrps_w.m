function grpNums = Func_GetIndvGrps_w(S)
grpNums = zeros(size(S,2), 1);
for ind = 1:size(S, 2)
    grpNums(ind) = sum( abs(S(:,ind)) > 1e-4 );
end
