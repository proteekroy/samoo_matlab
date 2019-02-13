% Copyright [2016] [Proteek Chandan Roy]
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
% Proteek Chandan Roy, 2016
% Contact: royprote@msu.edu, proteek_buet@yahoo.com



function [selected_pop, selected_pop_obj, selected_pop_cv] = memo_survival_selection(opt, pop, popObj, popCV, num)



    %--------------RETURN VALUES-------------------------------------------
    selected_pop = zeros(opt.N, opt.V);
    selected_pop_obj = zeros(opt.N, size(popObj,2));
    selected_pop_cv = zeros(opt.N, size(popCV,2));

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
            [~,I] = sort(obj(:,1)); %should be non-dominated sort for multiple objectives
            index{i} = index{i}(I,:);%store in a sorted order    
        end
    end

    
    %---------SELECTION BASED ON CLUSTERING AND MIN ASF VALUE--------------
    k = 1;
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


