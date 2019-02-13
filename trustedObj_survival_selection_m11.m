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

%THIS FUNCTION APPLIES SURVIVAL TOURNAMENT SELECTION

function [selected_pop, selected_pop_obj, selected_pop_cv] = trustedObj_survival_selection_m11(opt, pop, popObj, popCV)

    n  = size(opt.archive, 1);
    allpop = normalize(opt, opt.archive, opt.bound(1,:), opt.bound(2,:));
    norm_pop = normalize(opt, pop, opt.bound(1,:), opt.bound(2,:)); 
    div = pdist2(norm_pop, allpop);%distance in normalized space 
    trust_index = any(   div  <=  repmat(opt.TrustRadiusDeltaK(1:n), size(norm_pop,1),1),      2);
    mindiv = min(div,[],2);
    %I1 = find(index==1);
    %I2 = find(index==0);
    
    %consider only active set
    
    
%     allpop = normalize(opt, opt.archive(opt.activeArchiveIndex,:), opt.bound(1,:), opt.bound(2,:));
%     norm_pop = normalize(opt, pop, opt.bound(1,:), opt.bound(2,:)); 
%     div = pdist2(norm_pop, allpop);%distance in normalized space 
%     trust_index = any(   div  <=  repmat(opt.TrustRadiusDeltaK(1,opt.activeArchiveIndex), size(norm_pop,1),1),      2);
%     mindiv = min(div,[],2);
    
    

    index = popCV>0;
    popObj(index) = max(popObj)+ popCV(index);
    popCV = zeros(size(popObj,1),1);

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
        selectedIndividualIndex = -1; 

        if trust_index(p1)==1 && trust_index(p2)==1 %div(p1)<=delta && div(p2)<=delta% they are within trust region

            if (popCV(p1)<=0 && popCV(p2)<=0) || (popCV(p1)-popCV(p2))<1e-16 %both are feasible or same CV

                obj1 = popObj(p1,:);
                obj2 = popObj(p2,:);
                d = lex_dominate(obj1, obj2);

                if d == 1 %p1 dominates p2
                    selectedIndividualIndex = p1;
                elseif d == 3 %p2 dominates p1
                    selectedIndividualIndex = p2;
                else
                    selectedIndividualIndex = p1;
                end
            else
                if(popCV(p1)<popCV(p2))%p1 less constraint violation
                    selectedIndividualIndex = p1;
                else
                    selectedIndividualIndex = p2;
                end
            end
        else
            if mindiv(p1)<=mindiv(p2) %choose the one with nearest 
                selectedIndividualIndex = p1;
            else
                selectedIndividualIndex = p2;
            end
        end
        selected_pop(item, :) = pop(selectedIndividualIndex,1:opt.V); 
        selected_pop_obj(item,:) = popObj(selectedIndividualIndex,:);
        selected_pop_cv(item,:) = popCV(selectedIndividualIndex,:);
        item = item + 1;
    end


end