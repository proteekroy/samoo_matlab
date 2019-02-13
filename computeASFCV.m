%This function calculates ASF, ASFCV
function [ASF, ASFCV] = computeASFCV(opt, objective, CV)
    

    ASF = cell(1, opt.numdir);
    normalizedObj =  nsga3_normalization(opt, objective);
    
    for i=1:size(opt.dirs,1)
        w = opt.dirs(i,:);
        asf = calculate_Asf(normalizedObj, w);
        ASF{i} = asf; 
    end
    
    ASFCV = ASF;
    feasbile_index = CV<=0;
    
    for i=1:size(opt.dirs,1)
        if opt.C>0
            feasibleASF =  ASF{i}(feasbile_index,:);
            
            if ~isempty(feasibleASF)
                fmax = max(feasibleASF);%maximum feasible ASF
            else
                fmax = 0;
            end
            infeasible_index = find(CV>0);%infeasible indices
            if ~isempty(infeasible_index)
                ASFCV{i}(infeasible_index) = fmax + CV(infeasible_index);
            end
        end
    end
    
end