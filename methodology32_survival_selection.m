% Copyright [2018] [Proteek Chandan Roy]
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

%
% Proteek Chandan Roy, 2018
% Contact: royprote@msu.edu, proteek_buet@yahoo.com


function [selected_pop, selected_pop_obj, selected_pop_cv] = methodology32_survival_selection(opt, pop, popObj, popCV)

    selected_index = zeros(1, opt.N);
    
    if opt.trust_region_option_active==1
        
        switch(opt.trust_region_update_option)
            
            case 1 % continuous decreasing
                allpop = normalize(opt, opt.archive, opt.bound(1,:), opt.bound(2,:));
                norm_pop = normalize(opt, pop, opt.bound(1,:), opt.bound(2,:)); 
                div = pdist2(norm_pop, allpop);%distance in normalized space 
                div = min(div,[],2);
                I1 = find(div<=opt.TrustDistVar);
                %I2 = find(div>opt.TrustDistVar);
            case 2 % non-linear/diffusion trust region
                leader = normalize(opt, opt.archive(opt.ParetoIndex,:), opt.bound(1,:), opt.bound(2,:));
                if isempty(leader)
                    leader = normalize(opt, opt.archive, opt.bound(1,:), opt.bound(2,:));
                end
                norm_pop = normalize(opt, pop, opt.bound(1,:), opt.bound(2,:));
                div = trust_distance(opt, norm_pop, leader);
                div = min(div,[],2);
                I1 = find(div<=opt.delta);
                %I2 = find(div>opt.delta);                
            case 3 % adaptive
                n  = size(opt.archive, 1);
                allpop = normalize(opt, opt.archive, opt.bound(1,:), opt.bound(2,:));
                norm_pop = normalize(opt, pop, opt.bound(1,:), opt.bound(2,:)); 
                div = pdist2(norm_pop, allpop);%distance in normalized space 
                index = any(   div  <=  repmat(opt.TrustRadiusDeltaK(1:n), size(norm_pop,1),1),      2);
                I1 = find(index==1);
                %I2 = find(index==0);
                
            otherwise
                input('trust region option invalid');
        end
        
    else %no trust region option
        I1 = (1:1:size(pop, 1))';%every solution is valid
        %I2 = [];
    end 
    
    if size(I1,1)>opt.N % more solutions are inside trust region 
        
        index = find(popCV<=0);
        max_feasible_ASF = max(popObj(index,:));%for each direction, max feasible ASF 
        if isempty(index)
            max_feasible_ASF = zeros(1, opt.numdir);
        end
        infeasible_index = find(popCV>0);
        if ~isempty(infeasible_index) 
            popObj(infeasible_index,:) = popCV(infeasible_index,:) + repmat(max_feasible_ASF, size(infeasible_index,1), 1);
        end
        [~, sorted_index] = sort(popObj);


        count = 1;
        done = zeros(1, size(pop, 1));
        for i=1:size(popObj,1)
            for j=1:size(popObj,2) 
                s = sorted_index(i,j);
                if done(s)==0
                    done(s) = 1;
                    selected_index(count) = s;
                    count = count + 1;
                    if count > opt.N
                        break;
                    end
                end
            end
            if count > opt.N
                break;
            end
        end
    else % not enough solutions inside trust region
        selected_index = I1;
        temp_N = opt.N - size(I1,1);
        if temp_N>0
            div = min(div,[],2);
            [~, I3] = sort(div);
            I3 = setdiff(I3, I1);
            temp_selectedPopIndex = I3(1:temp_N);
            selected_index = vertcat(selected_index, temp_selectedPopIndex);
        end
        
    end
    selected_pop = pop(selected_index, 1:opt.V); 
    selected_pop_obj = popObj(selected_index, :);
    selected_pop_cv = popCV(selected_index, :);
    
%     %------------FIND CLUSTER NUMBER FOR THE SOLUTIONS---------------------
%     
%     archive = normalize(opt, opt.archive, opt.bound(1,:), opt.bound(2,:));
%     normalized_pop = normalize(opt, pop, opt.bound(1,:), opt.bound(2,:));
%     
%     closest_point_index = dsearchn(archive,normalized_pop);
%     pop_cluster = opt.archiveCluster(closest_point_index);
%     
%     
%     %---------------COLLECT CLUSTER MEMBERS AND SORT THEM------------------
%     index = cell(opt.numdir,1);
%     
%     for i = 1:opt.numdir %number of clusters
%         
%         index{i} = find(pop_cluster == i); %all solutions in cluster i
%         if ~isempty(index{i})
%             obj = popObj(index{i},:); %objectives of cluster i
%             [~,I] = sort(obj(:,1)); %should be non-dominated sort for multiple objectives
%             index{i} = index{i}(I,:);%store in a sorted order    
%         end
%     end
% 
%     
%     %---------SELECTION BASED ON CLUSTERING AND MIN ASF VALUE--------------
%     k = 1;
%     while(k<=opt.N)
%         
%         for j=1:opt.numdir
%             if ~isempty(index{j})
%                 selected_pop(k, :) = pop(index{j}(1),1:opt.V); 
%                 selected_pop_obj(k,:) = popObj(index{j}(1),:);
%                 selected_pop_cv(k,:) = popCV(index{j}(1),:);
%                 index{j}(1) = [];%delete from the cluster as taken
%                 k = k + 1;
%             end
%             if(k > opt.N)
%                 break;
%             end
%         end
%     end

    
    
    
%     %---------------COLLECT CLUSTER MEMBERS AND SORT THEM------------------
%     index = cell(opt.numdir,1);
%     
%     for i = 1:opt.numdir %number of clusters
%         [~, index{i}] = sort(popObj(:,i),'ascend');
%     end
%     
%     %---------SELECTION BASED ON CLUSTERING AND MIN ASF VALUE--------------
%     k = 1;
%     selected_index = [];
%     while(k<=opt.N)
%         
%         for j=1:opt.numdir
%             if ~ismember(index{j}(1), selected_index)
%                 selected_index = horzcat(selected_index, index{j}(1));
%                 selected_pop(k, :) = pop(index{j}(1),1:opt.V); 
%                 selected_pop_obj(k,:) = popObj(index{j}(1),:);
%                 selected_pop_cv(k,:) = popCV(index{j}(1),:);
%                 k = k + 1;
%             end
%             if k>opt.N
%                 break;
%             end
%             index{j}(1) = [];%delete from the cluster as taken
%         end
%         
%     end
% 
% end


