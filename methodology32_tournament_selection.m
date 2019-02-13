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



function selected_pop = methodology32_tournament_selection(opt, pop, popObj, popCV)



%--------------------TAKE MIN ASF AS OBJECTIVE VALUE-----------------------

popObj = min(popObj, [], 2);


N = opt.N;

%----TOURNAMENT CANDIDATES-------------------------------------------------

tour1 = randperm(N);
tour2 = randperm(N);


%----START TOURNAMENT SELECTION--------------------------------------------


selected_pop = zeros(N, opt.V); % Only the design variables of the selected members

for i = 1:N
    p1 = tour1(i);
    p2 = tour2(i);
    if(p1==p2)
       p2 = p1+1;
       if p2>N
           p2 = 1;%ensure two solutions are different
       end
    end
    
    cv1 = min(popCV(p1,:));
    cv2 = min(popCV(p2,:));
    
    if (cv1<=0 && cv2<=0) || (cv1-cv2)<1e-16 %both are feasible or same CV
        
        obj1 = popObj(p1,:);
        obj2 = popObj(p2,:);
        d = lex_dominate(obj1, obj2);
        
        if d == 1 %p1 dominates p2
            selected_pop(i, :) = pop(p1,1:opt.V);
        elseif d == 3 % p2 dominates p1
            selected_pop(i, :) = pop(p2,1:opt.V); 
        else
            if min(popObj(p1,:)) <= min(popObj(p2,:))
                selected_pop(i, :) = pop(p1,1:opt.V);%initially p1 was randomly choosen
            else
                selected_pop(i, :) = pop(p2,1:opt.V);
            end
        end
    else
        if(cv1<=cv2)%p1 less constraint violation
            selected_pop(i, :) = pop(p1,1:opt.V);
        else
            selected_pop(i, :) = pop(p2,1:opt.V); 
        end
    end
    
end


end