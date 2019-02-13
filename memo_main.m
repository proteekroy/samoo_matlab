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

function opt = memo_main(opt)

    %------------Initialize------------------------------------------------
    opt.nsga2.hifipop = opt.archive;
    opt.nsga2.hifipopObj = opt.archiveObj;
    opt.nsga2.hifipopCons = opt.archiveCons;
    opt = buildModel(opt);
    opt.nsga2.pop = lhsamp_model(opt.nsga2.initpopsize, opt);%LHS Sampling
    [opt.nsga2.popObj, opt.nsga2.popCV, opt.nsga2.popCons] = evaluate_pop_metamodel(opt, opt.nsga2.pop);
    [opt.nsga2.asf, opt.nsga2.pop_cluster] = find_memo_cluster(opt, opt.nsga2.popObj);
    
    %-------------------PLOT INITIAL SOLUTIONS----------------------------- 
    %plot_cluster(opt, opt.archiveObj, opt.archiveCluster);
    %plot_population(opt, opt.archiveObj);
    
    %--------------- OPTIMIZATION -----------------------------------------
    opt.nsga2.funcEval = opt.nsga2.initpopsize;
    progressiveInterval = opt.nsga2.buildModelInterval;
    iter = 0;
    while opt.funcEval <= opt.nsga2.totalFuncEval

        if(opt.nsga2.funcEval >= progressiveInterval)
            
            if(size(opt.nsga2.pop,1)+opt.funcEval>opt.nsga2.totalFuncEval)%discard excess solution in last iteration
                temppop = opt.nsga2.pop(randsample(opt.nsga2.N,opt.nsga2.totalFuncEval-opt.funcEval),:);
            else
                temppop = opt.nsga2.pop;
            end
            [temppopObj, temppopCons] = evaluate_pop(opt, temppop);
            tempCV = evaluateCV(opt, temppopCons);
            opt.nsga2.hifipop = vertcat(opt.nsga2.hifipop, temppop);
            opt.nsga2.hifipopObj = vertcat(opt.nsga2.hifipopObj, temppopObj);
            opt.nsga2.hifipopCons = vertcat(opt.nsga2.hifipopCons, temppopCons);          
            opt = buildModel(opt);
            [opt] = store_results(opt, temppop, temppopObj, temppopCons);
            opt.funcEval = size(opt.archive, 1);%number of function evaluations
            igd = igd_calculation(opt.paretofront, opt.archiveObj(opt.FeasibleIndex,:));
            gd = igd_calculation(opt.archiveObj(opt.ParetoIndex,:), opt.paretofrontGD);
            
            disp(['Method: ' num2str(opt.methodology) ', Eval: ' num2str(opt.funcEval) ',  IGD: ' num2str(igd) ...
            ',  GD: ' num2str(gd)  ', TD(var): ' num2str(opt.TrustDistVar) ', TD(obj): ' num2str(opt.TrustDistObj) ', Delta: ' num2str(opt.delta)]);
        
            if opt.funcEval>=opt.nsga2.totalFuncEval
                plot_population(opt, temppopObj(tempCV<=0,:));
                break;
            end
            
            if opt.trust_region_option_nsga2==1
                if (0.95*opt.TrustDistVar)<0.03
                    opt.TrustDistVar = 0.03;
                else
                    opt.TrustDistVar =  0.95*opt.TrustDistVar;
                end
                if opt.TrustDistObj*0.95 < 0.03
                    opt.TrustDistObj = 0.03;
                else
                    opt.TrustDistObj = opt.TrustDistObj*0.95;
                end
            elseif opt.trust_region_option_nsga2==2
                if opt.delta*0.95<0.005
                     opt.delta = 0.005;
                else
                    opt.delta = 0.95*opt.delta;
                end
            end
            progressiveInterval  =  opt.nsga2.funcEval + opt.nsga2.buildModelInterval;
            plot_cluster(opt, opt.archiveObj, opt.archiveCluster);
            plot_population(opt, temppopObj(tempCV<=0,:));
            
        end
        
        %-------------Mating Selection-------------------------------------
        opt.nsga2.popChild = memo_tournament_selection(opt, opt.nsga2.pop, opt.nsga2.popObj, opt.nsga2.popCV);%regular tournament selection, since we have obj and cv values
        
        %-------------------Crossover--------------------------------------
        [opt.nsga2.popChild, opt.nsga2.nrealcross] = sbx2(opt.nsga2.popChild, opt.nsga2.pcross, opt.nsga2.nrealcross, opt.nsga2.eta_c, opt.bound(1,:), opt.bound(2,:), opt.nsga2.Epsilon);
        
        %--------------------Mutation--------------------------------------
        [opt.nsga2.popChild, opt.nsga2.nrealmut] = pol_mut2(opt.nsga2.popChild, opt.nsga2.pmut, opt.nsga2.nrealmut,  opt.nsga2.eta_m,  opt.bound(1,:), opt.bound(2,:) );
        
        %---------------EVALUATION-----------------------------------------
        [opt.nsga2.popChildObj, opt.nsga2.popChildCV, opt.nsga2.popChildCons] = evaluate_pop_metamodel(opt, opt.nsga2.popChild);%get obj and constraints
        %[opt.nsga2.popChildObj, opt.nsga2.popChildCons] = evaluate_pop(opt, opt.nsga2.popChild);
        %opt.nsga2.popChildCV = evaluateCV(opt, opt.nsga2.popChildCons);
        
        %---------------MERGE----------------------------------------------
        opt.nsga2.totalpopObj = vertcat(opt.nsga2.popChildObj, opt.nsga2.popObj);
        opt.nsga2.totalpop = vertcat(opt.nsga2.popChild, opt.nsga2.pop);
        opt.nsga2.totalpopCV = vertcat(opt.nsga2.popChildCV, opt.nsga2.popCV);
        opt.nsga2.totalpopCons = vertcat(opt.nsga2.popChildCons, opt.nsga2.popCons);
        
        %---------------SELECTION------------------------------------------
        opt = trusted_memo_selection(opt);
        opt.nsga2.funcEval = opt.nsga2.funcEval + opt.nsga2.N;
               
        %-------------------PLOT NEW SOLUTIONS----------------------------- 
        plot_population(opt, opt.nsga2.popObj);
        iter = iter+1;
    end
        
    %------------------RETURN VALUE----------------------------------------
    %[opt.pop, opt.popObj] = calculate_feasible_paretofront(opt, opt.pop, opt.popObj, opt.popCV);

end%end of function
%------------------------------END OF -FILE--------------------------------

