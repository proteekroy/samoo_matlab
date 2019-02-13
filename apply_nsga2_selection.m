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
% Contact: royprote@egr.msu.edu, proteek_buet@yahoo.com


function [selectedPopIndex, CD] = apply_nsga2_selection(opt, pop, popObj, popCV, N)
    

    [n, ~] = size(popObj);
    R = zeros(size(popObj,1),1);   
     
    selectedPopIndex = [];
    index = popCV<=0;
    FeasiblePopIndex = find(index == 1);    
    InfeasiblePopIndex = find(index == 0);

    %---------------Find Non-dominated Sorting of Feasible Solutions-------
    F = cell(n,1);
    
    if ~isempty(FeasiblePopIndex)
                
        [RR,~] = bos(popObj(FeasiblePopIndex,:));
        
        for i=1:size(FeasiblePopIndex,1)
            F{RR(i)} = horzcat(F{RR(i)}, FeasiblePopIndex(i));
        end
    end 
    
    %---------------Store Ranking of Feasible Solutions--------------------
    if ~isempty(FeasiblePopIndex)
        R(FeasiblePopIndex) = RR;
    end
    %--------------Rank the Infeasible Solutions---------------------------
    
    if ~isempty(InfeasiblePopIndex)
        
        CV = popCV(InfeasiblePopIndex);

        [~,index] = sort(CV,'ascend');
        c = max(R) + 1;

        for i = 1: size(index,1)
            if i>1 && (CV(index(i))==CV(index(i-1)))%If both CV are same, they are in same front
                R(InfeasiblePopIndex(index(i))) = R(InfeasiblePopIndex(index(i-1)));
                b = R(InfeasiblePopIndex(index(i)));
                F{b} = horzcat(F{b}, InfeasiblePopIndex(index(i)));
            else
                R(InfeasiblePopIndex(index(i))) = c ;
                F{c} = horzcat(F{c}, InfeasiblePopIndex(index(i)));
                c = c + 1;
            end
        end
    
    end
    

    
    %----------------Select High Rank Solutions----------------------------
    count = zeros(n,1);
    for i=1:n
        count(i) = size(F{i},2);
    end

    cd = cumsum(count);
    p1 = find(N<=cd);
    lastfront = p1(1);
    
    for i=1:lastfront-1
        selectedPopIndex = horzcat(selectedPopIndex, F{i});
    end
    
    %------------CROWDING DISTANCE PART------------------------------------
    
    CD = zeros(size(popObj,1),1);
    for i=1:max(R)
        front = F{i};
        front_cd = crowdingDistance(opt, front, popObj(front,:));
        CD(front) = front_cd;
    end
    
   
    if size(selectedPopIndex,2)<N
        index = F{lastfront};
        CDlastfront = CD(index);

        [~,I] = sort(CDlastfront,'descend');

        j = 1;
        for i = size(selectedPopIndex,2)+1: N
            selectedPopIndex = horzcat(selectedPopIndex, index(I(j)));
            j = j + 1;
        end
    end
    
end