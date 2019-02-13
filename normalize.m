function [normalized_obj] = normalize(~, popObj, min_val, max_val)

    %{
    [~,ia,~] = unique(opt.archiveObj,'rows');%find unique objective values
    obj = opt.archiveObj(ia,:);%corresponding objective values
    index = ones(size(obj,1),1);
    if opt.C>0
        cv = opt.archiveCV(ia,:);%corresponding constraint violation 
        index = cv<=0;%
    end
    %}

    %min_val = opt.min_val;%[-300 0];
    %max_val = opt.max_val;%[1 1];%[-30 80];
    %{
    index = paretofront(opt.archiveObj);%find pareto fronts
    temp_pop = opt.archiveObj(index==1,:);

    if(size(temp_pop,1)==1)
        max_val = temp_pop;
        min_val = temp_pop;
    else
        max_val = max(temp_pop);
        min_val = min(temp_pop);
    end
    %}
    %max_val = max(opt.archiveObj);
    %min_val = max(opt.archiveObj);  
    
    if min(abs(max_val-min_val))>1e-16 %if values are not close
        normalized_obj = (popObj-repmat(min_val,size(popObj,1),1))./((repmat(max_val,size(popObj,1),1)+0.000001)-repmat(min_val,size(popObj,1),1));      
    else
        normalized_obj = popObj;
    end
    
    
end