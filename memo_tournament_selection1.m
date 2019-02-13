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



function selected_pop = memo_tournament_selection1(opt, pop, popObj, popCV)

    N = size(pop,1);
    k = opt.V;
    selected_pop = zeros(N, k);
    
    %-------------FIND CLUSTER NUMBERS-------------------------------------
    archive = normalize(opt, opt.archive, opt.bound(1,:), opt.bound(2,:));
    normalized_pop = normalize(opt, pop, opt.bound(1,:), opt.bound(2,:));
    
    closest_point_index = dsearchn(archive,normalized_pop);
    pop_cluster = opt.archiveCluster(closest_point_index);
    
    %closest_point_index = dsearchn(opt.archive(:,1:opt.V),pop(:,1:opt.V));%for each element in pop, find closest solution in x-space
    %pop_cluster = opt.archiveCluster(closest_point_index); %cluster values of those solutions
    
    %-------------SELECT MATING SOLUTIONS----------------------------------
    
    for p1=1:N %each solution will get chance for at least one tournament
        
        c = pop_cluster(p1);%cluster of the solution
        
        index = find(pop_cluster == c); %all solutions in cluster i
        index = setdiff(index, p1);%all points except p1
        
        
        if isempty(index)
            selected_pop(p1,1:k) = pop(p1,1:k);%This is the only element in that cluster
        else
            p = randi(size(index,1));%pick random solution for tournament
            p2 = index(p);

            
            if (popCV(p1)<=0 && popCV(p2)<=0) || (popCV(p1)-popCV(p2))<1e-16 %both are feasible or same CV
                obj1 = popObj(p1,:);
                obj2 = popObj(p2,:);
                d = lex_dominate(obj1, obj2);

                if d == 1 %p1 dominates p2
                    selected_pop(p1, :) = pop(p1,1:opt.V);
                elseif d == 3 % p2 dominates p1
                    selected_pop(p1, :) = pop(p2,1:opt.V); 
                else
                    selected_pop(p1, :) = pop(p1,1:opt.V);%initially p1 was randomly choosen
                end
            else
                if(popCV(p1)<popCV(p2))%p1 less constraint violation
                    selected_pop(p1, :) = pop(p1,1:opt.V);
                else
                    selected_pop(p1, :) = pop(p2,1:opt.V); 
                end
            end
        end

    end    


end

%------------------------------END OF -FILE--------------------------------
