% Copyright [2017] [Proteek Chandan Roy]
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
% Proteek Chandan Roy, 2017
% Contact: royprote@msu.edu, proteek_buet@yahoo.com

function [selected_pop, selected_pop_obj, selected_pop_cv] = trusted_survival_selection(opt, pop, popObj, popCV, num)

    %--------------RETURN VALUES-------------------------------------------
    selected_size = 0;
    selected_pop = zeros(opt.N, opt.V);
    selected_pop_obj = zeros(opt.N, size(popObj,2));
    selected_pop_cv = zeros(opt.N, size(popCV,2));

    
    
    %------------FIND WHO IS WITHIN TRUST REGION---------------------------
    
    %leader = opt.archive(opt.ParetoIndex,:);
    leaderAll = cell2mat(opt.LeaderIndex);
    leader = opt.archive(leaderAll,:);
    if isempty(leader)
        [~, I] = min(opt.archiveCV);
        leader =  opt.archive(I,:);
    end
    
    %normalize leaeders and population in x-space with bound
    leader = normalize(opt, leader, opt.bound(1,:), opt.bound(2,:));
    norm_pop = normalize(opt, pop, opt.bound(1,:), opt.bound(2,:)); 
    div = trust_distance(opt, norm_pop, leader);
    div = min(div,[],2);
    %div = div/opt.V;
    I1 = find(div<opt.delta);
    I2 = find(div>opt.delta);
    
    if size(I1,1)>opt.N
        pop = pop(I1,:);
        popObj = popObj(I1, :);
        popCV = popCV(I1,:);
    else
        selected_pop(1:size(I1,1), :) = pop(I1,1:opt.V); 
        selected_pop_obj(1:size(I1,1),:) = popObj(I1,:);
        selected_pop_cv(1:size(I1,1),:) = popCV(I1,:);
        selected_size = size(I1,1);
        pop = pop(I2,:);
        popObj = popObj(I2,:);
        popCV = popCV(I2,:);
    end
    
    
    
    %-------------FIND CLUSTER NUMBERS-------------------------------------
    
    archive = normalize(opt, opt.archive, opt.bound(1,:), opt.bound(2,:));
    normalized_pop = normalize(opt, pop, opt.bound(1,:), opt.bound(2,:));
    
    closest_point_index = dsearchn(archive,normalized_pop);
    pop_cluster = opt.archiveCluster(closest_point_index);
    
    
    %---------------COLLECT CLUSTER MEMBERS AND SORT THEM------------------
    index = cell(opt.numdir,1);
    
    for i = 1:opt.numdir %number of clusters
        
        index{i} = find(pop_cluster == i); %all solutions in cluster i
        if ~isempty(index{i})
            obj = popObj(index{i},:); %objectives of cluster i
            [~,I] = sort(obj(:,1),'ascend'); %should be non-dominated sort for multiple objectives
            
            
            index{i} = index{i}(I,:);%store in a sorted order    
        end
    end

    
    %---------SELECTION BASED ON CLUSTERING AND MIN ASF VALUE--------------
    k = selected_size+1;
    while(k<=num)
        
        for j=1:opt.numdir
            if ~isempty(index{j})
                selected_pop(k, :) = pop(index{j}(1),1:opt.V); 
                selected_pop_obj(k,:) = popObj(index{j}(1),:);
                selected_pop_cv(k,:) = popCV(index{j}(1),:);
                index{j}(1) = [];%delete from the cluster as taken
                k = k + 1;
            end
            if(k > num)
                break;
            end
        end
    end

end


