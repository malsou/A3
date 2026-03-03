function X = enforce_opr_bounds(X, varNames, lb, ub, OPR_min, OPR_max)
% Resample LPR/HPR for rows violating OPR bounds.
% Throws an error if a row cannot be repaired within maxTries.

iL = find(varNames=="LPR",1);
iH = find(varNames=="HPR",1);
if isempty(iL) || isempty(iH)
    error("LPR/HPR not found");
end

maxTries = 200;
for i = 1:size(X,1)
    opr = X(i,iL) * X(i,iH);
    tries = 0;
    while (opr < OPR_min || opr > OPR_max) && tries < maxTries
        % resample LPR/HPR uniformly within bounds
        X(i,iL) = lb(iL) + rand*(ub(iL)-lb(iL));
        X(i,iH) = lb(iH) + rand*(ub(iH)-lb(iH));
        opr = X(i,iL) * X(i,iH);
        tries = tries + 1;
    end

    if opr < OPR_min || opr > OPR_max
        error("Failed to enforce OPR bounds for row %d after %d retries (OPR=%.4f).", i, maxTries, opr);
    end
end
end