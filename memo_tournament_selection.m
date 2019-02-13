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



function selected_pop = memo_tournament_selection(opt, pop, popASF, popCV, popCluster)

    N = opt.N;
    selected_pop = zeros(N, opt.V);
    
    %-------------FIND CLUSTER NUMBERS-------------------------------------
    pop_cluster = popCluster;
    asf = popASF;
    
    l = 1;
    for i = 1:opt.numdir %number of clusters
        
        member = find(pop_cluster == i); %all solutions in cluster i
        if length(member) == 1
            selected_pop(l, :) = pop( member, 1:opt.V);
            selected_pop(l+1, :) = pop( member, 1:opt.V);
            l = l+2;
        elseif length(member) > 1
            tour1 = randperm(length(member));
            tour2 = randperm(length(member));
            
            for j=1:length(member)
                
                if popCV(member(tour1(j)))==0 && popCV(member(tour2(j)))==0 % both feasible
                    if asf(member(tour1(j)))<=asf(member(tour2(j)))
                        selected_pop(l, :) = pop(member(tour1(j)),1:opt.V);
                    else
                        selected_pop(l, :) = pop(member(tour2(j)),1:opt.V);
                    end
                else
                    if popCV(member(tour1(j)))<=popCV(member(tour2(j)))
                        selected_pop(l, :) = pop(member(tour1(j)),1:opt.V);
                    else
                        selected_pop(l, :) = pop(member(tour2(j)),1:opt.V);
                    end
                    
                end
                l = l+1;
                if l>N
                    break;
                end
            end
        end
        if l>N
            break;
        end
    end
    
    
end

%------------------------------END OF -FILE--------------------------------
