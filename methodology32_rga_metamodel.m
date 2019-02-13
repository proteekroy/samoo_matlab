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

function pop = methodology32_rga_metamodel(opt)

    
    pop = lhsamp_model(opt.N, opt);%LHC initialize
    [popObj, popCV, popCons] = evaluate_model_space(opt, pop);%objective value prediction
    
    %---------SIMPLE REAL PARAMETER GENETIC ALGORITHM----------------------
    opt.gen = 1;
    totalpopObj = zeros(2*opt.N, opt.numdir);
    totalpop = zeros(2*opt.N, opt.V);
    totalpopCV = zeros(2*opt.N, 1);
    
    while opt.gen <= opt.G
        %disp(opt.gen)
        %------------Mating Selection and Child Creation-------------------
        popChild = methodology32_tournament_selection(opt, pop, popObj, popCV);
        [popChild, opt.nrealcross] = sbx( popChild, opt.pcross, opt.nrealcross, opt.eta_c, opt.bound(1,:), opt.bound(2,:));%simulated binary crossover
        [popChild, opt.nrealmut] = pol_mut( popChild, opt.pmut, opt.nrealmut, opt.eta_m, opt.bound(1,:), opt.bound(2,:));%polynomial mutation

        %--------------Evaluate--------------------------------------------
        [popChildObj, popChildCV, popChildCons] = evaluate_model_space(opt, popChild);%objective value prediction
        %popChildCV = evaluateCV(opt, popChildCons);

        %--------------Merge-----------------------------------------------
        totalpopObj(1:opt.N,:) = popChildObj;
        totalpopObj(opt.N+1:end,:) = popObj;
        totalpop(1:opt.N,:) = popChild;
        totalpop(opt.N+1:end,:) = pop;
        totalpopCV(1:opt.N) = popChildCV;
        totalpopCV(opt.N+1:end) = popCV;
        
        %---------------Survival Selection---------------------------------
        [pop, popObj, popCV]  = methodology32_survival_selection(opt, totalpop, totalpopObj, totalpopCV);
        opt.gen = opt.gen+1;

    end
    %-----------UNIQUE SOLUTION--------------------------------------------
    %     [~,ia,~] = unique(pop,'rows');%find unique points
    %     pop = pop(ia,:);
    %     popObj = min(popObj(ia,:),[],2);
    %     popCV = min(popCV(ia,:),[],2);
    pop = pop(randperm(size(pop,1)),:);
    pop =  pop(1:opt.numdir,:);    

    
end


%------------------------------END OF -FILE--------------------------------
