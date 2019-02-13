function [opt, tempTrustRadiusDeltaK] = adjust_trust_region(opt, predictedpop, predictedObj, predictedCV, hifiObj, hifiCV, archive, archiveObj, archiveCV)
    
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
    
    %opt.TrustRadiusDeltaK
    %suppose, DOE point = p, new predicted point = q
    
    data_points = archiveObj;
    
%     total_archive_hifiObj = vertcat(archiveObj, hifiObj);
%     ub1 = max(total_archive_hifiObj);
%     hv_archive_hifiObj = hypeIndicatorExact(total_archive_hifiObj, ub1, 1 );
%     
%     total_archive_PredictedObj = vertcat(archiveObj, predictedObj);
%     ub2 = max(total_archive_PredictedObj);
%     hv_archive_PredictedObj = hypeIndicatorExact(total_archive_PredictedObj, ub2, 1 );
%     
%     
%     hv_archive = hv_archive_hifiObj(1:N);
%     hv_q = hv_archive_hifiObj(N+1:end);
%     
%     hv_predicted_q = hv_archive_PredictedObj(N+1:end);
    
    
    ub = max(data_points,[], 2);
    ub = ub';
    hv_data = lebesgue_measure(data_points, ub);
    %hv_data = hypeIndicatorExact(data_points, ub, 1 );
    tempTrustRadiusDeltaK = [];
    
    for i=1:n
           
        neighbor_size = size(neighbor{i}, 2);
        
        if neighbor_size>0
            
            deltaKneighbors = [];
            for j=1:min(5, neighbor_size)
                
                %disp(['Neighbor size = ' num2str(neighbor_size)]);
                neighbor_index = neighbor{i}(j);
                
                %-----------------FIND r-----------------------------------
                if archiveCV(neighbor_index)<=0 && hifiCV(i)<=0 %Both are feasible
                    data_with_q = vertcat(data_points, hifiObj(i, :));              
                    hv_with_q = lebesgue_measure(data_with_q, ub);
                    %hv_with_q = hypeIndicatorExact(data_with_q, ub, 1);
                    data_with_q_pred = vertcat(data_points, predictedObj(i, :));
                    hv_q_pred = lebesgue_measure(data_with_q_pred, ub);
                    %hv_q_pred = hypeIndicatorExact(data_with_q_pred, ub, 1);
                    r = (hv_data - hv_with_q)/(hv_data - hv_q_pred);
                    %r = (hv_archive(neighbor_index) - hv_q(i))/(hv_archive(neighbor_index) - hv_predicted_q(i));
                    
                elseif archiveCV(neighbor_index)>0 && hifiCV(i)<=0 %DOE point is infeasible
                    r = opt.TrustR2+0.1;%better point is found
                elseif archiveCV(neighbor_index)<=0 && hifiCV(i)>0%New point is infeasible
                    r = opt.TrustR1-0.1;%worse point is found
                elseif archiveCV(neighbor_index)>0 && hifiCV(i)>0%Both are infeasible
                    r = (archiveCV(neighbor_index) - hifiCV(i))/(archiveCV(neighbor_index) - predictedCV(i));
                end
                
                %s = div(i, neighbor_index);%normalized distance
                
                %---------------------Find Delta K-------------------------
                if r < opt.TrustR1
                    %tempDelta = opt.TrustC1*s;
                    tempDelta = opt.TrustC1*opt.TrustRadiusDeltaK(neighbor_index);
                elseif r > opt.TrustR2
                    tempDelta = min(opt.TrustC2*opt.TrustRadiusDeltaK(neighbor_index), opt.TrustDeltaStar);
                else
                   %tempDelta = s;
                   tempDelta = opt.TrustRadiusDeltaK(neighbor_index);
                end
                archiveDeltaK{neighbor_index} = vertcat(archiveDeltaK{neighbor_index}, tempDelta);
                deltaKneighbors = vertcat(deltaKneighbors, tempDelta);
            end
            
        end
           
        %assign delta-K for new points
        if neighbor_size>0
            tempTrustRadiusDeltaK = horzcat(tempTrustRadiusDeltaK, max(deltaKneighbors));
        else
            tempTrustRadiusDeltaK = horzcat(tempTrustRadiusDeltaK, opt.TrustInit);
        end        
        
    end
    
    %adjust delta-K for old points
    for i=1:N
        if size(archiveDeltaK{i},1)>0
            opt.TrustRadiusDeltaK(i) = min(archiveDeltaK{i});
        end
    end
    
    
end



