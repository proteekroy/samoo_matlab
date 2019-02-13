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
% Contact: proteek_buet@yahoo.com



function [opt] = metamodel_main_loop(opt)


%-------------------INITIAL HI-FI EVALUATION-------------------------------

    [pop, popObj, popCons] = initialization(opt); 
    opt = store_results(opt, pop, popObj, popCons);%ASF+CV+KKTPM+CLUSTER+PARETO+WRITE FILE+NORMALIZE
    
%-------------------MAIN LOOP----------------------------------------------
    opt.funcEval = size(pop,1);
    opt.iter = 1;
    opt.selection_function_option = opt.methodology;
    igd_array = zeros(opt.totalFuncEval, 1);
    
    
    while(opt.funcEval<opt.totalFuncEval)
        
        %clf(opt.fig)
        switch(opt.methodology)
            case 11%Model M11
                opt = methodology11(opt);
            case 12%Model M12
                opt = methodology12(opt);
            case 13%hybrid of Model M12 and M11
                opt = methodology1112(opt);
            case 122%Model M122
                opt = methodology122(opt);
            case 21%Model M21
                opt = methodology21(opt);
            case 22%Model M22
                opt = methodology22(opt);
            case 32%Model M22
                opt = methodology32(opt);
            case {3, 4}%Model M3-1 or M4-1    
                opt = methodology3(opt);
                %plot_points(opt);
                %plot_population(opt, opt.archiveObj(opt.FeasibleIndex));
            case 5%Model M5
                opt = methodology5(opt);
            case 6%Model KKTPM
                opt.selection_function_option = 2;
                opt = methodology6_kktpm(opt);                
            case 7%MODEL MEMO
                opt.selection_function_option = 3;
                opt = methodology6_memo(opt);
            otherwise
                disp('Methodology Not Implemeted');
        end
        
        %igd = igd_calculation(opt.paretofront, opt.archiveObj(opt.FeasibleIndex,:));
        %igd_array(opt.iter) = igd;
        %gd = igd_calculation(opt.archiveObj(opt.ParetoIndex,:), opt.paretofrontGD);
        %real_igd = igd_calculation(opt.paretofrontGD, opt.archiveObj(opt.FeasibleIndex,:));
        opt.funcEval = size(opt.archive,1);%number of function evaluations
        %disp(['Method: ' num2str(opt.methodology) ', Eval: ' num2str(opt.funcEval) ', real-IGD: ' num2str(real_igd) ...
        %    ',  GD: ' num2str(gd) ', Min(TD):' num2str(min(opt.TrustRadiusDeltaK)),  ', TD(var): ' num2str(opt.TrustDistVar) ', TD(obj): ' num2str(opt.TrustDistObj) ', Delta: ' num2str(opt.delta)]);
        
        opt.activeSetSize = min(opt.funcEval, opt.InitialActiveSetSize);% opt.funcEval;%min(opt.funcEval, 1000);

        opt.iter = opt.iter + 1;

        %==================PLOT============================================
        %if mod(opt.iter,10)==0
        %    plot_points(opt);
            %plot_trust_region(opt);
        %end
        
        %================SWITCH METHOD=====================================
%         if opt.methodology == 12 && opt.funcEval<opt.nsga2.totalFuncEval
%             opt.methodology = 12;
%         else
%             opt.methodology = 3;
%         end
        
    end %----------------------END OF OPTIMIZATION-------------------------
    
    %disp(igd_array);
    %dlmwrite('igd.txt', igd_array, 'delimiter','\n','precision','%.10f','-append');

    %--------FIND UNIQUE POINTS FOR ALL APPROACHES-------------------------
    [opt] = unique_population(opt); 
   
    
end



%------------------------------END OF -FILE--------------------------------