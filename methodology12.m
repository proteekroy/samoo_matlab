function opt = methodology12(opt)%using nsga2
    
    opt = metamodel(opt, 1:size(opt.archive, 1));%build model
    
    %-------------------RUN NSGA-II/NSGA-III-------------------------------
    if opt.methodology12_option==1
        [opt, pop] = nsga2_main(opt);
    else
        [opt, pop] = nsga3_main(opt);
    end
    %-----------------EVALUATE POP-----------------------------------------
    if(size(pop,1)+opt.funcEval>opt.totalFuncEval)%discard excess solution in last iteration
        pop = pop(1:opt.totalFuncEval-opt.funcEval,:);
        opt.nsga2.pop =  opt.nsga2.pop(1:opt.totalFuncEval-opt.funcEval,:);
        opt.nsga2.popObj = opt.nsga2.popObj(1:opt.totalFuncEval-opt.funcEval,:);
        opt.nsga2.popCV = opt.nsga2.popCV(1:opt.totalFuncEval-opt.funcEval,:);
        opt.nsga2.popCons = opt.nsga2.popCons(1:opt.totalFuncEval-opt.funcEval,:);
    end
                                 
    [popObj, popCons] = evaluate_pop(opt, pop);%high fidelity computation    
    [popCV, ~]= evaluateCV(popCons);
    
    %-----------------UPDATE TRUST RADIUS----------------------------------
    if opt.trust_region_update_option==3
        if opt.adaptive_trust_region_option==1
            [opt, tempTrustRadiusDeltaK] = adjust_trust_region(opt, opt.nsga2.pop, opt.nsga2.popObj, opt.nsga2.popCV, popObj, popCV, opt.archive, opt.archiveObj, opt.archiveCV);
        else
            [opt, tempTrustRadiusDeltaK] = adjust_trust_region_asf(opt, opt.nsga2.pop, opt.nsga2.popObj, opt.nsga2.popCV, popObj, popCV, opt.archive, opt.archiveObj, opt.archiveCV);
        end

        opt.TrustRadiusDeltaK = horzcat(opt.TrustRadiusDeltaK, tempTrustRadiusDeltaK);
    end
    %-----------------STORE RESULTS & PLOT---------------------------------
    opt = store_results(opt, pop, popObj, popCons);

end
