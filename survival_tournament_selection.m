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



function [selected_pop, selected_pop_obj, selected_pop_cv] = survival_tournament_selection(opt, pop, popObj, popCV)


N = opt.N;
N2 = N*2;

%----TOURNAMENT CANDIDATES-------------------------------------------------
tour = randperm(N2);

%----START TOURNAMENT SELECTION--------------------------------------------

selected_pop = zeros(N, opt.V); % Only the design variables of the selected members
selected_pop_obj = zeros(N, size(popObj,2));
selected_pop_cv = zeros(N, size(popCV,2));

item = 1;

for i = 1:2:N2
    p1 = tour(i);
    p2 = tour(i+1);
    if(p1==p2)
       p2 = p1+1;
       if p2>N2
           p2 = 1;%ensure two solutions are different
       end
    end
    if(p1>N2 || p2>N2)
        disp('Inside');
    end
    
    if (popCV(p1)<=0 && popCV(p2)<=0) || (popCV(p1)-popCV(p2))<1e-16 %both are feasible or same CV
        
        obj1 = popObj(p1,:);
        obj2 = popObj(p2,:);
        d = lex_dominate(obj1, obj2);
        
        if d == 1 %p1 dominates p2
            selected_pop(item,:) = pop(p1,1:opt.V);
            selected_pop_obj(item,:) = popObj(p1,:);
            selected_pop_cv(item,:) = popCV(p1,:);
        elseif d == 3 %p2 dominates p1
            selected_pop(item, :) = pop(p2,1:opt.V); 
            selected_pop_obj(item,:) = popObj(p2,:);
            selected_pop_cv(item,:) = popCV(p2,:);
        else
            selected_pop(item, :) = pop(p1,1:opt.V);%initially p1 was randomly choosen
            selected_pop_obj(item,:) = popObj(p1,:);
            selected_pop_cv(item,:) = popCV(p1,:);
        end
    else
        if(popCV(p1)<popCV(p2))%p1 less constraint violation
            selected_pop(item, :) = pop(p1,1:opt.V);
            selected_pop_obj(item,:) = popObj(p1,:);
            selected_pop_cv(item,:) = popCV(p1,:);
        else
            selected_pop(item, :) = pop(p2,1:opt.V); 
            selected_pop_obj(item,:) = popObj(p2,:);
            selected_pop_cv(item,:) = popCV(p2,:);
        end
    end
    item = item + 1;
end


end

%------------------------------END OF -FILE--------------------------------