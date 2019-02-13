

function [clusterNo] = cds_cluster(popCons)

    g = popCons;
    for i=1:size(pop_cons,2)
        g(g(:,i)<0,i)=0;
    end

    index = g>0;
    
    [~,ia,~] = unique(index,'rows');
    
    %cluster = cell(1,size(ia,1));
    clusterNo = zeros(1,size(popCons,1));
    
    for i=1:size(ia,1)%for all unique points
        for j=1:size(popCons,1)%for all solutions
            if isequal(index(i,:),index(j,:))%check if index values are equal
                %cluster{i} = horzcat(cluster{i},j);
                clusterNo(j) = i;
            end
        end
    end
    
    
    
end