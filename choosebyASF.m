function index = choosebyASF(opt, popObj, popCV)
    
    asfAll = cell(1, opt.numdir);
    asfCVAll = cell(1, opt.numdir);
    
    min_val = min(popObj);
    max_val = max(popObj);
    
    if min(abs(max_val-min_val))>1e-6
        normalized_obj = (popObj-repmat(min_val,size(popObj,1),1))./(repmat(max_val,size(popObj,1),1)-repmat(min_val,size(popObj,1),1));
    else
        normalized_obj = popObj;
    end
    
    for i=1:size(opt.dirs,1)
        w = opt.dirs(i,:);
        asf = calculate_Asf(normalized_obj, w);
        asfAll{i} = asf; 
    end
    
    %opt.archiveASFCVAll = opt.archiveASFAll;
    feasbile_index = popCV<=0;
    infeasible_index = find(popCV>0);%infeasible indices
    asfCVAll = asfAll;
    
    for i=1:size(opt.dirs,1)
        if opt.C>0
            feasibleASF =  asfAll{i}(feasbile_index,:);
            
            if ~isempty(feasibleASF)
                fmax = max(feasibleASF);%maximum feasible ASF
            else
                fmax = 0;
            end
            
            if ~isempty(infeasible_index)
                asfCVAll{i}(infeasible_index) = fmax + popCV(infeasible_index);
            end
        end
    end
    
    index = zeros(1, opt.numdir);
    
    for i=1:size(opt.dirs,1)
        [~, index(i)] = min(asfCVAll{i});
    end
    
    %index = unique(index);
    %pop = pop(index, :);
end

function asf = calculate_Asf(obj, dir)   
    asf = max( (obj-repmat(dir, size(obj,1),1)),[],2);
end
