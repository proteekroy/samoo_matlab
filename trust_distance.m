
% opt.epsilonTrust parameter controls smoothness of trusted region
% the larger epsilon, the isolated region it is

function div = trust_distance(opt, population, nondominaed_point)

    N = size(population, 1);
    data = vertcat(population, nondominaed_point);
    DistMatrix = squareform(pdist(data, 'euclidean'));
    K = exp(-opt.epsilonTrust*DistMatrix.^2);
    D = sum(K,2);
    P = K./repmat(D, 1, size(K,2));
    [U,S, ~]=svd(P);  % eigendecomposition of symmetric matrix
    S = diag(S);
    [S,I] = sort(S,'descend'); % sort according to magnitude of eigenvalues
    U = U(:, I);
    Map = U.*repmat(S', size(U,1), 1);
    diff_dist = squareform(pdist(Map,'euclidean'));
    div = min(diff_dist(:,N+1:end),[],2);
    div = div(1:N);
    
    %---------normalize diversity metric-----------------------------------
    %div = (div - min(div))./(max(div)-min(div));

end