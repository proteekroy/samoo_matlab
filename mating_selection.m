function [opt] = mating_selection(opt)

    switch(opt.matingselectionOption)
        case 1
            
            opt.popChild = constrained_tournament_selection(opt, opt.N);
    
        otherwise
            
            opt.popChild = constrained_tournament_selection(opt, opt.N);
    end

end