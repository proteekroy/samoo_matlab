function opt = high_fidelity_nsga2(opt)

    if(size(opt.archive,1)+size(opt.nsga2.pop,1) > opt.totalFuncEval)%discard excess solution in last iteration
        tempI = randsample(size(opt.nsga2.pop,1), opt.totalFuncEval-size(opt.archive,1));
        opt.nsga2.pop =  opt.nsga2.pop(tempI,:);
        opt.nsga2.popObj = opt.nsga2.popObj(tempI,:);
        opt.nsga2.popCV = opt.nsga2.popCV(tempI,:);
        opt.nsga2.popCons = opt.nsga2.popCons(tempI,:);
        opt.nsga2.CD = opt.nsga2.CD(tempI,:);
        temppop = opt.nsga2.pop;
    else
        temppop = opt.nsga2.pop;
    end
    
    [temppopObj, temppopCons] = evaluate_pop(opt, temppop);
    temppopCV = evaluateCV(temppopCons);
    
    opt.nsga2.funcEval = opt.nsga2.funcEval + size(temppop,1);
            
    %======================UPDATE TRUST RADIUS=============================
%     if opt.trust_region_option_nsga2==4
%         [opt, tempTrustRadiusDeltaK] = adjust_trust_region(opt, opt.nsga2.pop, opt.nsga2.popObj, opt.nsga2.popCV, temppopObj, temppopCV, opt.archive, opt.archiveObj, opt.archiveCV);
%         opt.TrustRadiusDeltaK = horzcat(opt.TrustRadiusDeltaK, tempTrustRadiusDeltaK);
%     end
    
            
    %=============================SAVE AND STORE===========================
    opt.nsga2.hifipop = vertcat(opt.nsga2.hifipop, temppop);
    opt.nsga2.hifipopObj = vertcat(opt.nsga2.hifipopObj, temppopObj);
    opt.nsga2.hifipopCons = vertcat(opt.nsga2.hifipopCons, temppopCons);     
    opt.nsga2.hifipopCV = vertcat(opt.nsga2.hifipopCV, temppopCV);
    opt = buildModel(opt);
    [opt] = store_results(opt, temppop, temppopObj, temppopCons);
    opt.funcEval = size(opt.archive, 1);%number of function evaluations    
    
    
end