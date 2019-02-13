function [opt, tempTrustRadiusDeltaK] = adjust_trust_region_asf(opt, predictedpop, predictedObj, predictedCV, hifiObj, hifiCV, archive, archiveObj, archiveCV)
    
    %here Objective value is ASF
    n = size(predictedpop, 1);
    N = size(archive, 1);
    
    allpop = normalize(opt, archive, opt.bound(1,:), opt.bound(2,:));
    norm_pop = normalize(opt, predictedpop, opt.bound(1,:), opt.bound(2,:)); 
    div = pdist2(norm_pop, allpop);%distance in normalized space   
    
    
    neighbor = cell(1, n);
    neighbor_dv = cell(1, n);
    archiveDeltaK = cell(1, N);
    
    for i=1:n
        neighbor{i} = find(div(i,:)<=opt.TrustRadiusDeltaK);
        neighbor_dv{i} = div(i, neighbor{i});
        [~, tempIndex] = sort(neighbor_dv{i});
        neighbor{i} = neighbor{i}(tempIndex);
    end
    
    tempTrustRadiusDeltaK = [];
    
    for i=1:n
           
        neighbor_size = size(neighbor{i}, 2);
        
        if neighbor_size>0
            
            deltaKneighbors = [];
            for j=1:min(5, neighbor_size)
                
                neighbor_index = neighbor{i}(j);
                
                %-----------------FIND r-----------------------------------
                if archiveCV(neighbor_index)<=0 && hifiCV(i)<=0 %Both are feasible
                    r = (archiveObj(neighbor_index) - hifiObj(i))/(archiveObj(neighbor_index) - predictedObj(i));
                elseif archiveCV(neighbor_index)>0 && hifiCV(i)<=0 %DOE point is infeasible
                    r = opt.TrustR2+0.1;%better point is found
                elseif archiveCV(neighbor_index)<=0 && hifiCV(i)>0%New point is infeasible
                    r = opt.TrustR1-0.1;%worse point is found
                elseif archiveCV(neighbor_index)>0 && hifiCV(i)>0%Both are infeasible
                    r = (archiveCV(neighbor_index) - hifiCV(i))/(archiveCV(neighbor_index) - predictedCV(i));
                end
                
                %s = div(i, neighbor_index);%normalized distance
                
                %---------------------Update Delta K-----------------------
                if r < opt.TrustR1
                    %tempDelta = opt.TrustC1*s;
                    tempDelta = opt.TrustC1*opt.TrustRadiusDeltaK(neighbor_index);
                elseif r > opt.TrustR2
                    tempDelta = min(opt.TrustC2*opt.TrustRadiusDeltaK(neighbor_index), opt.TrustDeltaStar);
                else
                    %tempDelta = opt.TrustC1*opt.TrustRadiusDeltaK(neighbor_index);
                    %tempDelta = s;
                    tempDelta = opt.TrustRadiusDeltaK(neighbor_index);
                end
                archiveDeltaK{neighbor_index} = vertcat(archiveDeltaK{neighbor_index}, tempDelta);
                deltaKneighbors = vertcat(deltaKneighbors, tempDelta);
            end 
        end
        
        if neighbor_size>0
            tempTrustRadiusDeltaK = horzcat(tempTrustRadiusDeltaK, max(deltaKneighbors));
        else
            tempTrustRadiusDeltaK = horzcat(tempTrustRadiusDeltaK, opt.TrustInit);
        end  
        
    end
    
    for i=1:N
        if size(archiveDeltaK{i},1)>0
            opt.TrustRadiusDeltaK(i) = min(archiveDeltaK{i});
        end
    end
    
    
end



