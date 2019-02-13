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

% Proteek Chandan Roy, 2018
% Contact: royprote@egr.msu.edu, proteek_buet@yahoo.com

function [opt, pop] = nsga3_main(opt)

    
    
    %-------------------PLOT INITIAL SOLUTIONS----------------------------- 
    %plot_population(opt, opt.archiveObj);
    
    %--------------- OPTIMIZATION -----------------------------------------
    opt.nsga2.funcEval = size(opt.archive,1);
    opt.nsga2.popCV = zeros(size(opt.nsga2.pop,1), 1);
           
    opt.nsga2.gen = 1;
    %------------Initialize------------------------------------------------
        
    opt.nsga2.pop = lhsamp_model(opt.nsga2.initpopsize, opt);%LHS Sampling
    [opt.nsga2.popObj, opt.nsga2.popCV, opt.nsga2.popCons] = evaluate_pop_metamodel(opt, opt.nsga2.pop);

    %-------------------FIND CLUSTER---------------------------------------
    %[~,opt.nsga2.CD] = apply_nsga3_selection(opt, opt.nsga2.pop, opt.nsga2.popObj, opt.nsga2.popCV, opt.nsga2.N);
        
    while opt.nsga2.gen < opt.nsga2.G
            %-------------Mating Selection---------------------------------
            opt.nsga2.popChild = constrained_tournament_selection_nsga3(opt, opt.nsga2.pop, opt.nsga2.popObj, opt.nsga2.popCV);
            %opt.nsga2.popChild = memo_tournament_selection(opt, opt.nsga2.pop, opt.nsga2.popASF, opt.nsga2.popCV, opt.nsga2.popCluster);

            %-------------------Crossover----------------------------------
            [opt.nsga2.popChild, opt.nsga2.nrealcross] = sbx2(opt.nsga2.popChild, opt.nsga2.pcross, opt.nsga2.nrealcross, opt.nsga2.eta_c, opt.bound(1,:), opt.bound(2,:), opt.nsga2.Epsilon);

            %--------------------Mutation----------------------------------
            [opt.nsga2.popChild, opt.nsga2.nrealmut] = pol_mut2(opt.nsga2.popChild, opt.nsga2.pmut, opt.nsga2.nrealmut,  opt.nsga2.eta_m,  opt.bound(1,:), opt.bound(2,:) );

            %---------------EVALUATION-------------------------------------
            [opt.nsga2.popChildObj, opt.nsga2.popChildCV, opt.nsga2.popChildCons] = evaluate_pop_metamodel(opt, opt.nsga2.popChild);

            %---------------MERGE------------------------------------------
            opt.nsga2.totalpopObj = vertcat(opt.nsga2.popChildObj, opt.nsga2.popObj);
            opt.nsga2.totalpop = vertcat(opt.nsga2.popChild, opt.nsga2.pop);
            opt.nsga2.totalpopCV = vertcat(opt.nsga2.popChildCV, opt.nsga2.popCV);
            opt.nsga2.totalpopCons = vertcat(opt.nsga2.popChildCons, opt.nsga2.popCons);

            %---------------SELECTION--------------------------------------
            opt = trusted_nsga3_selection(opt);
            %---------------COUNTER----------------------------------------
            opt.nsga2.gen = opt.nsga2.gen+1;
            %disp(opt.nsga2.gen);
    end
        
    %------------------RETURN VALUE----------------------------------------
    %[opt.pop, opt.popObj] = calculate_feasible_paretofront(opt, opt.pop, opt.popObj, opt.popCV);
    
    selectedPopIndex = choosebyASF(opt, opt.nsga2.popObj, opt.nsga2.popCV);
    
    %[selectedPopIndex,~] = apply_nsga2_selection(opt, opt.nsga2.pop, opt.nsga2.popObj, opt.nsga2.popCV, opt.numdir);
    %[~,opt.nsga2.CD] = apply_nsga3_selection(opt, opt.nsga2.totalpop, opt.nsga2.totalpopObj, opt.nsga2.totalpopCV, opt.nsga2.N);
     
    
    opt.nsga2.pop =  opt.nsga2.pop(selectedPopIndex,:);
    opt.nsga2.popObj = opt.nsga2.popObj(selectedPopIndex,:);
    opt.nsga2.popCV = opt.nsga2.popCV(selectedPopIndex,:);
    opt.nsga2.popCons = opt.nsga2.popCons(selectedPopIndex,:);
    %opt.nsga2.CD = opt.nsga2.CD(selectedPopIndex,:);
    
    
    pop = opt.nsga2.pop;

end%end of function
%------------------------------END OF -FILE--------------------------------

