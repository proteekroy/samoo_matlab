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



function pop_sel = kktpm_tournament_selection(opt, pop)


    N = size(pop,1);
    k = opt.V;
    pop_sel = zeros(N, opt.V);
    
    %find cluster info for rGA population
    closest_point_index = dsearchn(opt.activearchive(:,1:opt.V),pop(:,1:opt.V));%for each element in pop, find closest solution in x-space
    pop_cluster = opt.activearchiveCluster(closest_point_index); %cluster values of those solutions
    
    
    for i=1:N %for each solution
        
        c = pop_cluster(i);%cluster of the solution
        index = [];
        s=0;
        
        for j=1:N %collect all solution of same cluster
            if pop_cluster(j)==c && i~=j
                s = s+1;
                index(s) = j;
            end
        end
        if ~isempty(index)
            p = index(randi([1 size(index,2)]));%pick random solution for tournament
            if pop(i,k+1)<=pop(p,k+1)%less kktpm error is better
                pop_sel(i,1:k) = pop(i,1:k);
            else
                pop_sel(i,1:k) = pop(p,1:k);
            end
        else
            pop_sel(i,1:k) = pop(i,1:k);
        end
    end
    


end

%------------------------------END OF -FILE--------------------------------
