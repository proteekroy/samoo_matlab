function [opt] = update_trust_region(opt)



    switch(opt.trust_region_update_option)
      
        case 1
            %%{
            if opt.TrustDistVar*opt.TrustVarLR < opt.TrustDistVarMin
                opt.TrustDistVar = opt.TrustDistVarMin;
            else
                opt.TrustDistVar = opt.TrustDistVar*opt.TrustVarLR;
            end
            %}
            %opt.TrustDistVar = opt.TrustVarLR * opt.TrustDistVar;
            %opt.TrustVarLR = opt.TrustVarMin + 0.5*(opt.TrustVarMax - opt.TrustVarMin)*(1 + cos((opt.iter/opt.maxIter)*pi));
        case 2
            
            if opt.delta*0.9 < 0.0001
                opt.delta = 0.0001;
            else
                opt.delta = opt.delta*0.90;
            end
            %opt.TrustDelta = opt.TrustDeltaLR + opt.TrustDelta;
            %opt.TrustDeltaLR = opt.TrustMinDelta + 0.5*(opt.TrustMaxDelta - opt.TrustMinDelta)*(1 + cos((opt.iter/opt.maxIter)*pi));
            if opt.writeFlagOption==1
                dlmwrite(opt.TrustDeltaFileName, [opt.TrustDelta opt.TrustDeltaLR], 'delimiter','\t','precision','%.10f');        
            end
            
        case 3
            
            %don't update in a predefined way, update is adaptive
%             [opt, tempTrustRadiusDeltaK] = adjust_trust_region(opt, opt.nsga2.pop, opt.nsga2.popObj, opt.nsga2.popCV, temppopObj, temppopCV, opt.archive, opt.archiveObj, opt.archiveCV);
%             opt.TrustRadiusDeltaK = horzcat(opt.TrustRadiusDeltaK, tempTrustRadiusDeltaK);
            
        otherwise
            
            input('inside update_trust_region, option invalid');
            
        
    end
    
    if opt.trust_region_option_active==1
        if opt.funcEval>=opt.totalFuncEval
            if opt.writeFlag==1
                dlmwrite(opt.TrustDeltaFileName, opt.TrustRadiusDeltaK, 'delimiter','\n','precision','%.10f');        
            end
        end
    end



end