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


function opt = methodology6_nsga2(opt)%target multiple optimal solution


    %-------------------NORMALIZE------------------------------------------      
    
    opt.normalizedObj = normalize(opt, opt.archiveObj, opt.min_val, opt.max_val);
    [opt.selectionFunctionValue, ~] = selectionFunction(opt);%evaluate Selection value for all solutions, previous and new
    
                
    %--------------UNIQUES POPULATION--------------------------------------
                    
    opt = unique_population(opt);
                    
                    
    %----------MODEL COMBINED (ASF+CV)-------------------------------------
    
    
    %[opt.dmodel_selection, ~] = dacefit(opt.archive, opt.selectionFunctionValue, @regpoly1, @corrgauss, opt.theta, opt.lob, opt.upb);%model asf               
    %[opt.dmodel_asf, ~] = dacefit(opt.activeArchive, opt.activeArchiveASF, @regpoly1, @corrgauss, opt.theta, opt.lob, opt.upb);%model asf               
    %[opt.dmodel_asf, ~] = dacefit(opt.archive, opt.archiveASF, @regpoly2, @corrcubic, opt.theta, opt.lob, opt.upb);%model asf                
    opt.net = train_neural_network(opt.archive, opt.selectionFunctionValue);      
    %draw_selection_function_kriging(opt.dmodel_selection);
    %draw_selection_function_neural_net(opt.net);
    pop = memo_niched_rga_metamodel(opt);%one or two good solution
                
                
    if(size(pop,1)+opt.funcEval>opt.totalFuncEval)%discard excess solution in last iteration
        pop = pop(1:opt.totalFuncEval-opt.funcEval,:);
        
    end
                
                              
    %-----------------STORE RESULTS & PLOT-----------------------------
    [popObj, popCons] = evaluate_pop(opt, pop);%high fidelity computation    
    opt = store_results(opt, pop, popObj, popCons); %HI-FI+ASF+KKTPM+CLUSTER
    plot_population(opt, popObj);
    
end


%------------------------------END OF -FILE--------------------------------

