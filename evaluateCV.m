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

function [cv, cv_modified]= evaluateCV(pop_cons)
    

    g = pop_cons;
    
    for i=1:size(pop_cons,2)
        g(g(:,i)<0,i)=0;
    end
    cv = sum(g, 2);
    cv_modified = cv;

    feasible_index = find(cv==0);
    g2 = sum(pop_cons, 2);
    cv_modified(feasible_index) = g2(feasible_index); 

    
end