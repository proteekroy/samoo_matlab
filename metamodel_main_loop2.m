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

function [data, opt] = metamodel_main_loop2(opt)


%-------------------INITIAL HI-FI EVALUATION-------------------------------

    [pop, popObj, popCons] = initialization(opt);
    
%-----------------EVALUATE INITIAL SOLUTIONS-------------------------------
    
    if opt.initOption == 1 %if not evaluated while initialization
        [popObj, popCons] = evaluate_pop(opt, pop); 
    end
    
    opt = store_results(opt, pop, popObj, popCons);%ASF+CV+KKTPM+CLUSTER+PARETO+WRITE FILE+NORMALIZE
    %plot_population(opt, popObj);
    
    
%-------------------MAIN LOOP----------------------------------------------
    opt.funcEval = size(pop,1);
    opt.iter = 1;
    opt.selection_function_option = opt.methodology;
    
    igd_array = zeros(opt.totalFuncEval, 1);
    while(opt.funcEval<opt.totalFuncEval)
        
        switch(opt.methodology)
            case 1%Model M1
                opt = methodology11(opt);
            case 2%Model M2
                disp('Methodology 2 not implemented');
            case {3,4}%Model M3/M4    
                opt = methodology3(opt);
            case 5%Model M5
                opt = methodology5(opt);
            case 6%Model KKTPM
                opt.selection_function_option = 2;
                opt = methodology6_kktpm(opt);                
            case 7 %MODEL MEMO
                opt.selection_function_option = 3;
                opt = methodology6_memo(opt);
                %{
                if opt.funcEval>opt.switch_method
                    opt.InitialActiveSetSize = 50;
                    opt.methodology = 1;
                end
                %}
            case 8 %MODEL TrustRegion-RGA
                if opt.funcEval<opt.switch_method
                    opt.hybrid = 8;
                    opt.selection_function_option = 3;
                    %opt = methodology6_TrustRegion(opt);
                    opt = methodology_hybrid1_TrustRegion(opt);                
                else %if opt.funcEval>1000
                    opt.hybrid = 3;
                    opt.activeSetSize = 50;
                    opt = methodology3(opt);
                end
            
            case 9 %MODEL NSGA-II
                opt.selection_function_option = 4;
                opt = methodology6_nsga2(opt);
            otherwise
                disp('Methodology Out Of Bound Exception');
        end
        
        igd = IGD(opt.paretofront', opt.archiveObj(opt.FeasibleIndex,:)');
        igd_array(opt.iter) = igd;
        gd = IGD(opt.archiveObj(opt.ParetoIndex,:)', opt.paretofrontGD');
        %real_igd = IGD(opt.paretofrontGD', opt.archiveObj');
        opt.funcEval = size(opt.archive,1);%number of function evaluations
        %disp(['Function Evaluation: ' num2str(opt.funcEval) ',  IGD: ' num2str(igd) ',  GD: ' num2str(gd) ', Delta: ' num2str(opt.delta), ', real IGD: ' num2str(real_igd)]);
        disp(['Function Evaluation: ' num2str(opt.funcEval) ',  IGD: ' num2str(igd) ',  GD: ' num2str(gd) ', Trust Distance: ' num2str(opt.TrustDistVar)]);
        
        
        opt.activeSetSize = min(opt.funcEval, opt.InitialActiveSetSize);% opt.funcEval;%min(opt.funcEval, 1000);
        
        
        %==========UPDATE DELTA============================================
        if opt.funcEval> opt.changeDelta
            opt.changeDelta = opt.changeDelta+200;                 

            if opt.delta*0.5<0.005
                 opt.delta = 0.05;
            else
                opt.delta = opt.delta/2;%0.002;
            end
        end
        opt.iter = opt.iter + 1;
        %opt.TrustDistObj = (1-opt.TrustDistLearningRate)*opt.TrustDistObj;
        %==================PLOT====================================
        %if mod(opt.iter,10)==0
        %    plot_points(opt);
            %plot_trust_region(opt);
        %end
        
    end %----------------------END OF OPTIMIZATION-------------------------
    
    %if opt.methodology==8
    %    disp(igd_array);
    %end
%-----------------WRITE IGD TO FILE----------------------------------------

    %dlmwrite('igd.txt', pop, 'delimiter','\n','precision','%.10f','-append');
%--------FIND UNIQUE POINTS FOR ALL APPROACHES-----------------------------

    [opt] = unique_population(opt);
    
%---------------------------RETURN VALUES----------------------------------
    
    data.var = opt.archive;%return value, solution
    data.obj = opt.archiveObj;%return value, high fidelity objectives   
    data.cv = opt.archiveCV;%return constraint violation function
    data.cons = opt.archiveCons;%return exact constraint violation
    
    
end



%------------------------------END OF -FILE--------------------------------