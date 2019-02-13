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



function selected_pop = memo_tournament_selection2(opt, pop, popObj, popCV)

    N = opt.N;
    k = opt.V;
    selected_pop = zeros(N, k);
    
    %-------------FIND CLUSTER NUMBERS-------------------------------------
    pop_cluster = niching(opt,pop);
%     archive = normalize(opt, opt.archive, opt.bound(1,:), opt.bound(2,:));
%     normalized_pop = normalize(opt, pop, opt.bound(1,:), opt.bound(2,:));
%     
%     closest_point_index = dsearchn(archive,normalized_pop);
%     pop_cluster = opt.archiveCluster(closest_point_index);
    
    %closest_point_index = dsearchn(opt.archive(:,1:opt.V),pop(:,1:opt.V));%for each element in pop, find closest solution in x-space
    %pop_cluster = opt.archiveCluster(closest_point_index); %cluster values of those solutions
    
    
    
    
    %---------------COLLECT CLUSTER MEMBERS AND SORT THEM------------------
    %index = cell(1, opt.numdir);
    k = 1;
    for i = 1:opt.numdir %number of clusters
        
        member = find(pop_cluster == i); %all solutions in cluster i
        if ~isempty(member)
            index{k} = member;
            k = k + 1;
        end       
    end
    k = k - 1;
    
    
    %---------SELECTION BASED ON CLUSTERING AND MIN ASF VALUE--------------
    l = 1;
    while(l <= N)
        
        for j = 1:k % for each cluster
            
            %-----------PICK FIRST PARENT----------------------------------
            if isempty(index{j})%cluster cannot be empty
                disp('memo_tournament_selection is wrong');
                continue; 
            end
            
            i1 = randi(size(index{j},1));
            p1 = index{j}(i1);
            set = setdiff(index{j}, p1);%all points except p1
            if ~isempty(set)
                i2 = randi(size(set,1));
                p2 = set(i2);
                p = compare(p1, p2, popCV, popObj);
            else
                p = p1;
            end
            selected_pop(l, :) = pop(p,1:opt.V);
            
            %----------PICK SECOND PARENT----------------------------------
            cluster_set = setdiff(1:k, j);
            if isempty(cluster_set)%there is only one cluster
                j2 = j;
            else
                i = randi(size(cluster_set,1));
                j2 = cluster_set(i);%cluster number
            end
            
            i1 = randi(size(index{j2},1));%random solution from cluster
            p1 = index{j2}(i1);
            set = setdiff(index{j2}, p1);%all points except p1
            if ~isempty(set)
                i2 = randi(size(set,1));
                p2 = set(i2);
                p = compare(p1, p2, popCV, popObj);
            else
                p = p1;
            end
            
            selected_pop(l+1, :) = pop(p,1:opt.V); 
            l = l + 2;
            if l>N
                break
            end
            
        end
        
    end
    
end

%---------------- BINARY TOURNAMENT COMPARISON-----------------------------
function [p] = compare(p1, p2, popCV, popObj)

    
    if (popCV(p1)<=0 && popCV(p2)<=0) || (popCV(p1)-popCV(p2))<1e-16 %both are feasible or same CV
        obj1 = popObj(p1,:);
        obj2 = popObj(p2,:);
        d = lex_dominate(obj1, obj2);

        if d == 1 % p1 dominates p2
            p = p1;
        elseif d == 3 % p2 dominates p1
            p =  p2 ;
        else
            p = p1;
        end
    else
        if(popCV(p1)<popCV(p2))%p1 less constraint violation
            p = p1;
        else
            p = p2; 
        end
    end

end


%------------------------------END OF -FILE--------------------------------
