function X = morris_design(k, p, r)
% Generate Morris trajectories in normalized space [0,1]^k
% Output X: N x k, with N = r*(k+1)
%
% Uses delta = p/(2*(p-1)) on a p-level grid in [0,1].
% Step direction is randomized and flipped if it would leave [0,1].

levels = linspace(0, 1, p);
delta  = p/(2*(p-1));  % classic Morris step (e.g., p=6 -> 0.6)

X = zeros(r*(k+1), k);
row = 1;

for t = 1:r
    % Base point on grid (any level allowed)
    x0 = levels(randi(numel(levels), [1, k]));
    order = randperm(k);

    Xt = zeros(k+1, k);
    Xt(1,:) = x0;

    for i = 1:k
        j = order(i);

        Xt(i+1,:) = Xt(i,:);

        % Random direction (+/-). Flip if out of bounds.
        dir = 1;
        if rand > 0.5, dir = -1; end

        cand = Xt(i,j) + dir*delta;
        if cand < 0 || cand > 1
            dir = -dir;
            cand = Xt(i,j) + dir*delta;
        end

        % If still out of bounds (should be rare), clamp.
        cand = min(max(cand, 0), 1);

        % Snap to nearest grid level to avoid floating error
        [~, idx] = min(abs(levels - cand));
        Xt(i+1,j) = levels(idx);
    end

    X(row:row+k, :) = Xt;
    row = row + (k+1);
end
end