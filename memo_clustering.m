function [opt, pop3, cluster_sizes] = memo_clustering(opt, pop3)

    k = opt.V;
    M = opt.M;
    ncon = opt.C;
    p = opt.numdir-1;
    %pop = zeros( size(pop1) );
    extremeARC = zeros(opt.M, opt.V+opt.M+1);
    pop3 = [pop3; extremeARC];

    %=====================NORMALIZATION====================================
    [normF, extremeARC] = memo_normalize(pop3, k, M, extremeARC, ncon);
    %======================================================================

    if ncon > 0
        feasID = find(pop3(:, k+M+1) == 0);
        Nfeas = size(feasID,1);
        if Nfeas > 0 % If there is at least one feasible solution, move P to its location
            [~, popBestFeas] = min( sum( normF(feasID,:), 2 ) );
            P = shiftRefSet(NaN, normF(feasID(popBestFeas), :), M, p);
        else % If there are no feasible solutions, move P to the location of the solution that has the minimum CV value.
            infeasID = find(pop3(:, k+M+1) ~= 0);
            [~, popBestInfeas] = min( pop3(infeasID, k+M+1) );
            P = shiftRefSet(NaN, normF(infeasID(popBestInfeas), :), M, p);
        end
    else
        [~, popBest] = min( sum( normF, 2 ) );
        P = shiftRefSet(NaN, normF(popBest, :), M, p);
    end
    P=P';
    Np = size(P,1);
    N = size(pop3,1);

    [pop3, cluster_sizes] = computeMinASF_ClusterIDs(opt, Np, N, pop3, k, M, P, normF, ncon);
    opt.nsga2.popASF = pop3(1:end-opt.M,k+M+2);
    opt.nsga2.popCluster = pop3(1:end-opt.M,k+M+3);
    pop3 = pop3(1:end-opt.M,:);%get rid of extremeARC

end