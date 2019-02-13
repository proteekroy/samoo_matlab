
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

function pop = methodology1_rga_metamodel(opt)

    
    pop = lhsamp_model(opt.N, opt);%LHC initialize
    [popObj, popCons] = evaluate_pop_metamodel(opt, pop);%objective value prediction
    popCV = evaluateCV(opt, popCons);%  
    
    %tic
    %---------SIMPLE REAL PARAMETER GENETIC ALGORITHM----------------------
    opt.gen = 1;
    totalpopObj = zeros(2*opt.N,1);
    totalpop = zeros(2*opt.N,opt.V);
    totalpopCV = zeros(2*opt.N,1);
    totalpopCons = zeros(2*opt.N,1);
    
    %tic
    while opt.gen <= opt.G

        %------------Mating Selection and Child Creation-------------------
        popChild = constrained_tournament_selection(opt, pop, popObj, popCV);
        [popChild, opt.nrealcross] = sbx( popChild, opt.pcross, opt.nrealcross, opt.eta_c, opt.bound(1,:), opt.bound(2,:));%simulated binary crossover
        [popChild, opt.nrealmut] = pol_mut( popChild, opt.pmut, opt.nrealmut, opt.eta_m, opt.bound(1,:), opt.bound(2,:));%polynomial mutation

        %--------------Evaluate--------------------------------------------
        [popChildObj, popChildCV, popChildCons] = evaluate_pop_metamodel(opt, popChild);%objective value prediction
        %popChildCV = evaluateCV(opt, popChildCons);

        %--------------Merge-----------------------------------------------
        totalpopObj(1:opt.N,:) = popChildObj;
        totalpopObj(opt.N+1:end,:) = popObj;
        totalpop(1:opt.N,:) = popChild;
        totalpop(opt.N+1:end,:) = pop;
        totalpopCV(1:opt.N,:) = popChildCV;
        totalpopCV(opt.N+1:end,:) = popCV;
        totalpopCons(1:opt.N,:) = popChildCons;
        totalpopCons(opt.N+1:end,:) = popChildCons;
        
        
        
        %---------------Survival Selection---------------------------------
        [pop, popObj, popCV]  = survival_tournament_selection_m1(opt, totalpop, totalpopObj, totalpopCV,totalpopCons);

        opt.gen = opt.gen+1;

    end

    %t = toc;
    %disp(t);
    %-----------UNIQUE SOLUTION--------------------------------------------
    [~,ia,~] = unique(pop,'rows');%find unique points
    pop = pop(ia,:);
    popObj = popObj(ia,:);
    popCV = popCV(ia,:);
    
    %------------RETURN ONE SOLUTION---------------------------------------
    index = find(popCV<=0);
    if isempty(index)%all are infeasible
        [~,index] = min(popCV);
    else %more than one feasible
        [~,index2] = sort(popObj(index,:));
        index = index(index2(1));
    end
    
    pop = pop(index,1:opt.V);

    
end


%------------------------------END OF -FILE--------------------------------
