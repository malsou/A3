function X = lhs_simple(n, k)
% Latin hypercube in [0,1] (basic)
X = zeros(n,k);
for j = 1:k
    edges = linspace(0,1,n+1);
    u = rand(n,1);
    X(:,j) = edges(1:n)' + u.*(edges(2:n+1)'-edges(1:n)');
    X(:,j) = X(randperm(n),j);
end
end
