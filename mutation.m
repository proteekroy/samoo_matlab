function [opt] = mutation(opt)

    switch(opt.mutationOption)
        case 1
            
            [opt.popChild, opt.nrealmut] = pol_mut(opt.popChild, opt.p_mut, opt.nrealmut,  opt.eta_m,  opt.bound(1,:), opt.bound(2,:) );
             
        otherwise
            
            [opt.popChild, opt.nrealmut] = pol_mut(opt.popChild, opt.p_mut, opt.nrealmut,  opt.eta_m,  opt.bound(1,:), opt.bound(2,:) );
    end

end