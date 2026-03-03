function X = enforce_opr_bounds(X, varNames, lb, ub, OPR_min, OPR_max)
% Resample LPR/HPR for any rows violating OPR bounds

iL = find(varNames=="LPR",1);
iH = find(varNames=="HPR",1);
if isempty(iL) || isempty(iH)
    error("LPR/HPR not found");
end

for i = 1:size(X,1)
    opr = X(i,iL) * X(i,iH);
    tries = 0;
    while (opr < OPR_min || opr > OPR_max) && tries < 200
        % resample LPR/HPR uniformly within bounds
        X(i,iL) = lb(iL) + rand*(ub(iL)-lb(iL));
        X(i,iH) = lb(iH) + rand*(ub(iH)-lb(iH));
        opr = X(i,iL) * X(i,iH);
        tries = tries + 1;
    end
end
end
