
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

function opt = methodology6_kktpm(opt)%generative, one by one


    %opt.modality = 3;
    %-------------------NORMALIZE------------------------------------------      
    
    opt.normalizedObj = normalize(opt, opt.archiveObj, opt.min_val, opt.max_val);
    %opt.normalizedObj =  nsga3_normalization(opt, opt.archiveObj);
    [opt.archiveASF, opt.archiveCluster] = evaluateASF(opt);%evaluate ASF for all solutions, previous and new
    
    %----------------FIND ACTIVE SET-----------------------------------
    opt = find_active_set(opt);         
    
    %--------------UNIQUES POPULATION--------------------------------------         
    %opt = unique_population(opt);                
    [~,ia,~] = unique(opt.archive,'rows');
    temp_archive = opt.archive(ia,:);
    temp_archiveKKTPM = opt.activeArchiveKKTPM(ia, :);
    
    %---------------------MODEL KKTPM--------------------------------------
    %if opt.M<3
    [opt.dmodel_kktpm, ~] = dacefit(temp_archive, temp_archiveKKTPM, @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model kktpm surface
    %else
    %[opt.dmodel_kktpm, ~] = dacefit(opt.activeArchive, opt.activeArchiveKKTPM, @regpoly1, @corrgauss, opt.theta, opt.lob, opt.upb);%model kktpm surface
    %end
    %pop = memo_niched_rga_metamodel(opt);%one or two good solution
    %opt.net = train_neural_network(opt.archive, opt.archiveKKTPM);
    pop = kktpm_niched_rga_metamodel(opt);%same as previous code      
                
    if(size(pop,1)+opt.funcEval>opt.totalFuncEval)%discard excess solution in last iteration
        pop = pop(1:opt.totalFuncEval-opt.funcEval,:);
    end
                
                              
    %-----------------STORE RESULTS & PLOT-----------------------------
    [popObj, popCons] = evaluate_pop(opt, pop);%high fidelity computation    
    opt = store_results(opt, pop, popObj, popCons); %HI-FI+ASF+KKTPM+CLUSTER
    
    if opt.plotOption==1
        plot_points(opt);       
    end
                
end


%------------------------------END OF -FILE--------------------------------