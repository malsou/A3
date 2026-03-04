function mor = morris_analyze(Xn, Y, k, r, delta)
%MORRIS_ANALYZE Morris elementary effects for one or more responses.
% Returns (k x m):
%   mu     = mean(di)            (signed mean effect)
%   mu_star= mean(abs(di))       (mean absolute effect)
%   sigma2 = var(di,1)           (population variance, normalized by N)
%   sigma  = sqrt(sigma2)
%
% Categorical factors are expected to be handled outside this function as
% separate scenarios. This function only analyzes continuous dimensions in Xn.

[N, k2] = size(Xn);
if k2 ~= k
    error("k mismatch");
end
if N ~= r*(k+1)
    error("Row count mismatch: expected %d rows from r*(k+1), got %d.", r*(k+1), N);
end
m = size(Y,2);

EE = NaN(r, k, m); % trajectory, variable, response

idx = 1;
for t = 1:r
    Xt = Xn(idx:idx+k, :);
    Yt = Y(idx:idx+k, :);

    for step = 1:k
        dx = Xt(step+1,:) - Xt(step,:);
        changed = find(abs(dx) > 1e-12);
        if numel(changed) ~= 1
            error("Invalid Morris step in trajectory %d, step %d (changed vars=%d).", t, step, numel(changed));
        end

        j = changed(1);
        stepSize = dx(j);
        if abs(abs(stepSize) - delta) > 1e-8
            warning("Step size %.6f differs from expected delta %.6f (trajectory %d, step %d).", stepSize, delta, t, step);
        end

        ee = (Yt(step+1,:) - Yt(step,:)) / stepSize;
        tinyImag = abs(imag(ee)) <= 1e-12;
        ee(~tinyImag) = NaN;
        ee(tinyImag) = real(ee(tinyImag));
        ee(~isfinite(ee)) = NaN;  % robust handling of failed model calls
        EE(t,j,:) = ee;
    end
    idx = idx + (k+1);
end

mu = NaN(k,m);
mu_star = NaN(k,m);
sigma2 = NaN(k,m);

for j = 1:k
    for q = 1:m
        d = squeeze(EE(:,j,q));
        d = d(isfinite(d));
        if isempty(d)
            mu(j,q) = 0;
            mu_star(j,q) = 0;
            sigma2(j,q) = 0;
        else
            mu(j,q) = mean(d);
            mu_star(j,q) = mean(abs(d));
            % Use var(di,1): population variance (normalization by N)
            sigma2(j,q) = var(d, 1);
        end
    end
end

mor.EE = EE;
mor.mu = mu;
mor.mu_star = mu_star;
mor.sigma2 = sigma2;
mor.sigma = sqrt(max(sigma2,0));
end
