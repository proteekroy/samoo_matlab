function [opt] = crossover(opt)

    switch(opt.crossoverOption)
        case 1
            
            [opt.popChild, opt.nrealcross] = sbx(opt.popChild, opt.p_cross, opt.nrealcross, opt.eta_c, opt.bound(1,:), opt.bound(2,:));
    
        otherwise
            
            [opt.popChild, opt.nrealcross] = sbx(opt.popChild, opt.p_cross, opt.nrealcross, opt.eta_c, opt.bound(1,:), opt.bound(2,:));
    end

end