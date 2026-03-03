function X = denormalize_minmax(Xn, lb, ub)
X = lb + Xn .* (ub - lb);
end
