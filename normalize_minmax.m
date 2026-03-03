function Xn = normalize_minmax(X, lb, ub)
Xn = (X - lb) ./ (ub - lb);
end
