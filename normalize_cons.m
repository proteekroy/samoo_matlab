function [normalized_cons] = normalize_cons(popCons, popCV)

%     n = size(popCons, 1);
%     normalized_cons = popCons./(repmat(mean(popCons),n,1));
    
    %%{
    normalized_cons = zeros(size(popCons));
    
    for i=1:size(popCons,2)
        index = popCons(:, i)<= 0;
        feasibleCons = popCons(index,i);
        infeasibleCons = popCons(~index,i);
        
        if ~isempty(feasibleCons)
            max_val = max(feasibleCons);
            min_val = min(feasibleCons);
            if size(feasibleCons,1) > 2               
                n = size(feasibleCons,1);
                %normalized_cons(index, i) = -1e-15 - (feasibleCons-repmat(min_val,n,1))./((repmat(max_val,n,1)+1e-16)-repmat(min_val,n,1));
                normalized_cons(index, i) = -(feasibleCons./(repmat(mean(feasibleCons),n,1)));
            end
        end


        if ~isempty(infeasibleCons)
            max_val = max(infeasibleCons);
            min_val = min(infeasibleCons);
            if size(infeasibleCons,1) > 2              
                n = size(infeasibleCons,1);
                %normalized_cons(~index, i) = 1e-15 + (infeasibleCons-repmat(min_val,n,1))./((repmat(max_val,n,1)+1e-16)-repmat(min_val,n,1));
                normalized_cons(~index, i) = infeasibleCons./(repmat(mean(infeasibleCons),n,1));
            end
        end
        
    end
    %}
        
end