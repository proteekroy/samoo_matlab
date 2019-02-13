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
% Contact: royprote@egr.msu.edu, proteek_buet@yahoo.com


function [selectedPopIndex] = apply_memo_selection(opt, popObj, pop_cluster, N)
    

    %---------------COLLECT CLUSTER MEMBERS AND SORT THEM------------------
    index = cell(opt.numdir,1);
    
    for i = 1:opt.numdir %number of clusters
        
        index{i} = find(pop_cluster == i); %all solutions in cluster i
        if ~isempty(index{i})
            obj = popObj(index{i}); %objectives of cluster i
            [~,I] = sort(obj); %should be non-dominated sort for multiple objectives
            index{i} = index{i}(I,:);%store in a sorted order    
        end
    end

    selectedPopIndex = zeros(N, 1);
    %---------SELECTION BASED ON CLUSTERING AND MIN ASF VALUE--------------
    k = 1;
    while(k<=N)
        
        for j=1:opt.numdir
            if ~isempty(index{j})
                selectedPopIndex(k, :) = index{j}(1);
                index{j}(1) = [];%delete from the cluster as taken
                k = k + 1;
            end
            if(k > N)
                break;
            end
        end
    end
        
        
        
    
end