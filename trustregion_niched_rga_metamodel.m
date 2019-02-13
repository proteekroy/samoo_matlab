
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
% Contact: royprote@msu.edu, proteek_buet@yahoo.com

function pop = trustregion_niched_rga_metamodel(opt)


    %----------INCLUDE FEASIBLE PARETO SOLUTION FOR LOW-FI EVALUATION------
    pop = lhsamp_model(opt.N, opt);%LHC initialize
    [popObj, popCV, popCons] = evaluate_pop_metamodel(opt, pop);%objective value prediction
    %popCV = evaluateCV(opt, popCons);
    
    %---------SIMPLE REAL PARAMETER GENETIC ALGORITHM----------------------
    opt.gen = 1;
    totalpopObj = zeros(2*opt.N,1);
    totalpop = zeros(2*opt.N,opt.V);
    totalpopCV = zeros(2*opt.N,1);
    
    while opt.gen <= opt.G

        %------------Mating Selection and Child Creation-------------------
        popChild = memo_tournament_selection2(opt, pop, popObj, popCV);
        [popChild, opt.nrealcross] = memo_sbx(opt, popChild, opt.pcross, opt.nrealcross, opt.eta_c, opt.bound(1,:), opt.bound(2,:));%simulated binary crossover
        [popChild, opt.nrealmut] = pol_mut( popChild, opt.pmut, opt.nrealmut, opt.eta_m, opt.bound(1,:), opt.bound(2,:));%polynomial mutation
        
        %--------------Evaluate--------------------------------------------
        [popChildObj, popChildCons] = evaluate_pop_metamodel(opt, popChild);%objective value prediction
        popChildCV = evaluateCV(opt, popChildCons);

        %--------------Merge-----------------------------------------------
        totalpopObj(1:opt.N,:) = popChildObj;
        totalpopObj(opt.N+1:end,:) = popObj;
        totalpop(1:opt.N,:) = popChild;
        totalpop(opt.N+1:end,:) = pop;
        totalpopCV(1:opt.N,:) = popChildCV;
        totalpopCV(opt.N+1:end,:) = popCV;
        
        %---------------Survival Selection---------------------------------
        [pop, popObj, popCV]  = trusted_survival_selection(opt, totalpop, totalpopObj, totalpopCV, opt.N);
        opt.gen = opt.gen+1;

    end


    %-----------UNIQUE SOLUTION--------------------------------------------
    [~,ia,~] = unique(pop,'rows');%find unique points
    pop = pop(ia,:);
    popObj = popObj(ia,:);
    
    
    
    %----------------------------------------------------------------------
    %------RETURN AT MOST 1 SOLUTION FROM 1 CLUSTER------------------------
    %----------------------------------------------------------------------
    
    
    %-------------FIND CLUSTER NUMBERS-------------------------------------
    
    
%     archive = normalize(opt, opt.archive, opt.bound(1,:), opt.bound(2,:));
%     normalized_pop = normalize(opt, pop, opt.bound(1,:), opt.bound(2,:));
%     
%     closest_point_index = dsearchn(archive,normalized_pop);%for each element in pop, find closest solution in x-space
%     pop_cluster = opt.archiveCluster(closest_point_index); %cluster values of those solutions
    pop_cluster = niching(opt, pop);
    pop = clustering(opt, pop_cluster, popObj, pop);
    %[~,I] = sort(popObj);  
    %pop = pop(I(1:opt.numdir),:);
    %---------------REMOVE DUPLICATE---------------------------------------
    %{
    remove = [];
    for i = 1:size(pop,1)
        for j = 1:opt.numdir 
            if ~isempty(opt.bestPop{j})
                eud = norm(opt.bestPop{j}-pop(i,:));
                if(eud<1e-6) %exclude those which are near than euclidean distance 10^-6 in x-space
                    remove = horzcat(remove,i);
                    break;
                end
            end
        end    
    end
    if ~isempty(remove)
        pop(remove,:) = [];
    end
    %}
end

%------------------------------END OF -FILE--------------------------------

