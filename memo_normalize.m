function [normF, extremeARC] = memo_normalize(pop, k, M, extremeARC, ncon)

if ncon == 0 % Unconstrained Problem
    
    N = size(pop,1);
    % Step-1: F - Fmin
    min_F = min( pop(:, k+1:k+M) );
    normF = pop(:, k+1:k+M) - repmat(min_F, N, 1);
    % Step-2: Find the extremes
    epsn = 1e-6;
    rho = 1e-4;
    % extremeARC = zeros(M, k+M); % Archieve for the extreme values
    for i=1:M
        weight = epsn * ones(1,M);
        weight(i) = 1.0;
        rho_matrix = normF ./ repmat(weight, N, 1); % for constraint handling, normF should be composed of only feasible solutions!!!
        rhosum = sum(rho_matrix, 2);
        [~, extremeID] = min(max(rho_matrix, [], 2) + rho * rhosum);
        extremeARC(i, 1:k+M+1) = pop(extremeID, 1:k+M+1);
    end
    % Step-3: Adjust if extremes are absurd. Identify ND points of pop
    if ( rank( extremeARC(:, k+1:k+M) ) == M ) % Gauss elimination works
        alpha = ones(1,M) ./ ( pinv( extremeARC(:, k+1:k+M) - repmat(min_F, M, 1)) * ones(M,1) )';
        index = alpha < epsn;
        if any( index > 0 )
            disp('Gauss elimination works, but alpha < epsn')
            NDpop = find_pareto(pop(:, k+1:k+M));  % constrained nondomination check!!!
            nadir = max( pop(NDpop, k+1:k+M) ); % find max of ND points
            alpha(index) = nadir(index) - min_F(index);
        end
    else      % Gaussian elimination fails
        %disp('Gaussian elimination fails')
        NDpop = find_pareto(pop(:, k+1:k+M)); % constrained nondomination check!!!
        nadir = max( pop(NDpop, k+1:k+M) ); % find max of ND points
        alpha = nadir - min_F;
    end
    alpha = alpha + 1e-10;
    normF = ( pop(:, k+1:k+M) - repmat(min_F, N, 1)) ./ repmat(alpha, N, 1);
%     feasPopID = [];
    
else % Constrained Problem
    
    feasPopID = find(pop(:,k+M+1) == 0);
    Nfeas = size(feasPopID, 1);
    
    if Nfeas > 0 % If there are at least one feasible solution
        
        N = size(pop,1);
        % Step-1: F - Fmin
        min_F = min( pop(feasPopID, k+1:k+M) );
        if Nfeas == 1, min_F = pop(feasPopID, k+1:k+M); end
        normF = pop(:, k+1:k+M) - repmat(min_F, N, 1);
        % Step-2: Find the extremes
        epsn = 1e-6;
        rho = 1e-4;
        % extremeARC = zeros(M, k+M); % Archieve for the extreme values
        for i=1:M
            weight = epsn * ones(1,M);
            weight(i) = 1.0;
            rho_matrix = normF ./ repmat(weight, N, 1); % for constraint handling, normF should be composed of only feasible solutions!!!
            rhosum = sum(rho_matrix, 2);
%             [~, extremeID] = min(max(rho_matrix, [], 2) + rho * rhosum);
%             extremeARC(i, 1:k+M+1) = pop( extremeID, 1:k+M+1 );
            temp = max(rho_matrix, [], 2) + rho * rhosum;
            [~, extremeID] = min(temp(feasPopID));
            extremeARC(i, 1:k+M+1) = pop( feasPopID(extremeID), 1:k+M+1 );
        end
        % Step-3: Adjust if extremes are absurd. Identify ND points of pop
        if ( rank( extremeARC(:, k+1:k+M) ) == M ) % Gauss elimination works
            alpha = ones(1,M) ./ ( pinv( extremeARC(:, k+1:k+M) - repmat(min_F, M, 1)) * ones(M,1) )';
            index = alpha < epsn;
            if any( index > 0 )
                disp('Gauss elimination works, but alpha < epsn')
                NDpop = find_pareto(pop(feasPopID, k+1:k+M));  % constrained nondomination check!!!
                nadir = max( pop(feasPopID(NDpop), k+1:k+M) ); % find max of ND points
                if Nfeas == 1, nadir = pop(feasPopID(NDpop), k+1:k+M); end
                alpha(index) = nadir(index) - min_F(index);
            end
        else      % Gaussian elimination fails
            %disp('Gaussian elimination fails')
            NDpop = find_pareto(pop(feasPopID, k+1:k+M)); % constrained nondomination check!!!
            nadir = max( pop(feasPopID(NDpop), k+1:k+M) ); % find max of ND points
            alpha = nadir - min_F;
        end
%         normF = ( pop(feasPopID, k+1:k+M) - repmat(min_F, Nfeas, 1)) ./ repmat(alpha, Nfeas, 1);
        alpha = alpha + 1e-10;
        normF = normF ./ repmat(alpha, N, 1);
        
    else % If there are no feasible solutions
        
        N = size(pop,1);
        % Step-1: F - Fmin
        min_F = min( pop(:, k+1:k+M) );
        normF = pop(:, k+1:k+M) - repmat(min_F, N, 1);
        % Step-2: Find the extremes
        epsn = 1e-6;
        rho = 1e-4;
        % extremeARC = zeros(M, k+M); % Archieve for the extreme values
        for i=1:M
            weight = epsn * ones(1,M);
            weight(i) = 1.0;
            rho_matrix = normF ./ repmat(weight, N, 1); % for constraint handling, normF should be composed of only feasible solutions!!!
            rhosum = sum(rho_matrix, 2);
            [~, extremeID] = min(max(rho_matrix, [], 2) + rho * rhosum);
            extremeARC(i, 1:k+M) = pop(extremeID, 1:k+M);
        end
        % Step-3: Adjust if extremes are absurd. Identify ND points of pop
        if ( rank( extremeARC(:, k+1:k+M) ) == M ) % Gauss elimination works
            alpha = ones(1,M) ./ ( pinv( extremeARC(:, k+1:k+M) - repmat(min_F, M, 1)) * ones(M,1) )';
            index = alpha < epsn;
            if any( index > 0 )
                disp('Gauss elimination works, but alpha < epsn')
                NDpop = find_pareto(pop(:, k+1:k+M));  % constrained nondomination check!!!
                nadir = max( pop(NDpop, k+1:k+M) ); % find max of ND points
                alpha(index) = nadir(index) - min_F(index);
            end
        else      % Gaussian elimination fails
            %disp('Gaussian elimination fails')
            NDpop = find_pareto(pop(:, k+1:k+M)); % constrained nondomination check!!!
            nadir = max( pop(NDpop, k+1:k+M) ); % find max of ND points
            alpha = nadir - min_F;
        end
        alpha = alpha + 1e-10;
        normF = ( pop(:, k+1:k+M) - repmat(min_F, N, 1)) ./ repmat(alpha, N, 1);
        
    end
    
end


