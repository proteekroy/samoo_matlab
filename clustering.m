

function [selected_pop] = clustering(opt, pop_cluster, asf, pop)


    %---------------COLLECT CLUSTER MEMBERS AND SORT THEM------------------
    index = cell(1, opt.numdir);
    
    for i = 1:opt.numdir %number of clusters
        index{i} = find(pop_cluster == i); %all solutions in cluster i  
        obj = asf(index{i},:); %objectives of cluster i
        [~,I] = sort(obj(:,1)); %should be non-dominated sort for multiple objectives
        index{i} = index{i}(I,:);%store in a sorted order        
    end

    
    %---------SELECTION BASED ON CLUSTERING AND MIN ASF VALUE--------------
    
    selected_pop = zeros(opt.numdir, opt.V);
    k = 1;
    while k<= opt.numdir
        for j=1:opt.numdir
            if ~isempty(index{j})
                selected_pop(k, :) = pop(index{j}(1),1:opt.V);
                k = k + 1;
            end
            if k>opt.numdir
                break;
            end
        end
    end
    %selected_pop(k:opt.numdir,:) = [];
    

end