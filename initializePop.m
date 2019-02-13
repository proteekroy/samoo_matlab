function pop = initializePop(N, k, Xmin, Xmax)

X = zeros(N, k);
for i=1:k
    X(1:N, i) = repmat(Xmin(i), N, 1) + (Xmax(i) - Xmin(i)) * rand(N,1);
end
pop = zeros(N, k+1+1); % [x, f, TotConsViol]
pop(:, 1:k) = X;