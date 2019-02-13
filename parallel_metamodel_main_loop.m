
function [opt] = parallel_metamodel_main_loop(opt)

%-------------------INITIAL HI-FI EVALUATION-------------------------------

    [pop, popObj, popCons] = initialization(opt); 
    
%-------------------CALCULATE & SAVE RESULTS-------------------------------    
    
    opt = store_results(opt, pop, popObj, popCons);
    
%-------------------MAIN LOOP----------------------------------------------

    opt.funcEval = size(pop,1);
    opt.iter = 1;
    
    while(opt.funcEval<opt.totalFuncEval)
        
        %--------------------PARALLEL CROSS-VALIDATION---------------------
        if opt.num_of_frameworks > 1
            opt.methodology = runCrossValidation(opt);
            disp(['Selected Framework: ' num2str(opt.methodology)]);
        end
        
        %----------------------RUN-----------------------------------------
        
        %clf(opt.fig)        
        switch(opt.methodology)
            case {11, 21, 31, 41, 5}%Model M1-1, M2-1
                opt = methodology11(opt);
            case {12, 22}%Model M1-2, M2-2
                opt = methodology12(opt);
            %case {31, 41}%Model M3-1 or M4-1
            %    opt = methodology31(opt);
            case {32, 42}%Model M3-2 or M4-2
                opt = methodology32(opt);
%             case 5%Model M5
%                 opt = methodology5(opt);
            case 6%Model M6
                if opt.selection_function_option==2
                    opt = methodology6_memo(opt);%Model Multimodal EMO
                elseif opt.selection_function_option==1
                    opt = methodology6_kktpm(opt);%Model KKTPM
                end
            otherwise
                input('Inside paralel_metamodel_main_loop Methodology Not Implemeted');
        end
        
        opt.funcEval = size(opt.archive,1);%number of function evaluations
        
        opt = compute_performance_metric(opt, opt.archiveObj, opt.archiveCV);
        disp(['Framework: ' num2str(opt.methodology) ', Eval: ' num2str(opt.funcEval) ',  IGD: ' num2str(opt.igd) ',  GD: ' num2str(opt.gd)]);
        %disp(['TrustRegion: ' num2str(opt.TrustDistVar) ', Delta: ' num2str(opt.delta) ', TrustRegion: min: ' num2str(min(opt.TrustRadiusDeltaK)) ', avg: ' num2str(mean(opt.TrustRadiusDeltaK)), ' max: ' num2str(max(opt.TrustRadiusDeltaK))]);
        opt.iter = opt.iter + 1;
        opt = update_trust_region(opt);
        
        if opt.plotOption==1
            plot_points(opt);
        end

    end %----------------------END OF OPTIMIZATION-------------------------
   
    
end



%------------------------------END OF -FILE--------------------------------
