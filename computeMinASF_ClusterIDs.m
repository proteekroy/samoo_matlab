function [pop, cluster_sizes] = computeMinASF_ClusterIDs(opt, Np, N, pop, k, M, P, normF, ncon)

ASF = zeros( Np, N); % Initialize the ASF matrix
for i=1:Np
    ASF = compute_ASF( ASF, P(i,:), i, N, normF);
end

if ncon == 0
    [minASF, ASF_idx] = min( ASF ); % MinASFs - for clustering
    pop(:, k+M+2) = minASF';
    pop(:, k+M+3) = ASF_idx'; % Cluster IDs
else
    feasID = find(pop(:, k+M+1) == 0);
    infeasID = find(pop(:, k+M+1) ~= 0);
    Nfeas = size(feasID,1);
    
    if Nfeas > 0
        % Assign MinASF value for the feasible solutions
        [minASF, ASF_idx] = min( ASF );
        pop(feasID, k+M+2) = minASF(feasID)';
        pop(:, k+M+3) = ASF_idx'; % Cluster ID is assigned for ALL solutions !!!!!
        
        %find maximum ASF of feasible solutions
        maxASF = max(pop(feasID, k+M+2));%proteek
        
        % Assign minimum CV value for infeasible solutions
        pop(infeasID, k+M+2) = pop(infeasID, k+M+1); % Minimum CV
        %pop(infeasID, k+M+2) = maxASF + pop(infeasID, k+M+1);%proteek % max(ASF)+CV
        
    else % If all solutions are infeasible
        [~, ASF_idx] = min( ASF ); % This is performed only to determine the cluster IDs
        pop(:, k+M+2) = pop(:, k+M+1); % Minimum CV
        pop(:, k+M+3) = ASF_idx';
        
    end
end


% CLUSTER (in obj. space)
cluster_sizes = zeros(1,Np);
for i=1:Np
    for j = 1:N
        if pop(j, end) == i
            cluster_sizes( i ) = cluster_sizes( i ) + 1;
        end
    end
end

end

