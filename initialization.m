function [pop, popObj, popCons] = initialization(opt)

    switch(opt.objfunction)
        
        case {'ddmop1','ddmop2','ddmop3','ddmop4','ddmop5','ddmop6','ddmop7'}
            pop = DDMOP1('init');
            popObj = DDMOP1('value', pop);
            popCons = zeros(size(pop,1),1);
        otherwise
            
            switch(opt.initOption)
                case 1
                    %pop = load(opt.varfilename);
                    pop = lhsamp_model(opt.initpopsize, opt);%latin hypercube
                    [popObj, popCons] = evaluate_pop(opt, pop); 
                case 2
                    [pop , popObj, popCons] = incremental_density_based_sampling(opt.initpopsize, opt);%roulette wheel

                otherwise
                    disp('Undefined Initialization Option');
            end
    end

end