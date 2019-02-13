function M = BuildStochasticMatrix(opt, population)

 
    switch(opt.UtilityOption)
        case 1 
            n = size(population,1);%size of population;
            M = sparse(n, n);%large sqaure matrix
            u = computeUtility1(population);%suppose this is L1-norm
                        
            for i = 1:n
                for j = 1:n
                    if(i~=j)
                        if u(i)>u(j)
                            M(j, i) = u(i) - u(j);%direction from (i) to (j)
                        elseif u(i)<u(j)
                            M(i, j) = u(j) - u(i);
                        end
                    end
                end
            end
            
            %normalize to make stochastic
            for i = 1:n
                M(:,i) = M(:,i)./sum(M(:,i));
            end
              
        case 2
            n = size(population,1);%size of population
            M = sparse(n, n);%declare a sparse matrix
            Mdl = KDTreeSearcher(population,'Distance','euclidean');
            u = computeUtility1(population);%suppose this is L1-norm
            d = size(population,2);
            [Idx,~] = knnsearch(Mdl,population,'K',d+1);%k-d tree
            for i = 1:n
                for j = 2:d+1 %1st member is itself
                    s = Idx(j);
                    if(i~=s)
                        if u(i)>u(s)
                            %direction from (i) to (j), value iteration is reverse, it takes indegree
                            %rather than outdegree
                            M(s, i) = u(i) - u(s);
                        elseif u(i)<u(s)
                            M(i, s) = u(s) - u(i);
                        end
                    end
                end
            end
            
            %normalize M to make column stochastic
            for i = 1:n
                M(:,i) = M(:,i)./sum(M(:,i));
            end
            
        otherwise
            
            
    end

end

function u = computeUtility1(population)
    u = sum(population,2); %should be normalized   ----   %%sum(population.*0.5,2);
end