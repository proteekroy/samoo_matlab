% Copyright [2017] [Proteek Chandan Roy]
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

%
% Proteek Chandan Roy, 2017
% Contact: royprote@msu.edu, proteek_buet@yahoo.com


% opt.epsilonTrust parameter controls smoothness of trusted region
% the larger epsilon, the isolated region it is

function div = trust_region_distance(epsilonTrust, population, nondominaed_point)

    N = size(population, 1);
    data = vertcat(population, nondominaed_point);
    DistMatrix = squareform(pdist(data, 'euclidean'));
    K = exp(-epsilonTrust*DistMatrix.^2);
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