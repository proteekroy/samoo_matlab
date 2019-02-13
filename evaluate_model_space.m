function [popObj, popCV, popCons] = evaluate_model_space(opt, pop)

    %--------------EVALUATE OBJECTIVES/AGGREGATION-------------------------
    
    switch(opt.methodology)
        
        case {11, 12, 21, 22} %objectives are separately metamodeled
            
            if all(strcmpi(opt.current_obj_metamodel, 'GP'))
                
                [popObj, ~] = predictor(pop(:,1:opt.V), opt.dmodel_obj);
                
            else
                
                popObj = zeros(size(pop,1), opt.M);
                
                for i=1:opt.M
                    
                    switch(upper(opt.current_obj_metamodel{i}))
                        
                        case 'GP' 
                            
                            [popObj(:,i), ~] = predictor(pop, opt.dmodel_obj{i});
                            
                        case 'RBF'
                            
                        case 'SVR'
                    end
                end
            end    
            if opt.M>1 && (opt.methodology==11 || opt.methodology==21) %compute the single objective to optimize

                switch(upper(opt.generative_framework_acquisition_func{1}))

                    case 'ASF'

                        popObj = max( (popObj-repmat(opt.curdir, size(popObj,1),1)),[],2);
                end
            end
            
            
        case {31, 41} %%aggregated function of objectives is metamodeled
            
            switch(upper(opt.current_obj_aggregation_metamodel{1}))
                
                case 'GP'
                    
                    [popObj,~] = predictor(pop(:,1:opt.V), opt.dmodel_asf);
                    
                case 'RBF'
                            
                case 'SVR'
                    
            end
            
        case {32, 42}
            
            switch(upper(opt.current_obj_aggregation_metamodel{1}))
                
                case 'GP' 
                    
                    [popObj,~] = predictor(pop(:,1:opt.V), opt.dmodel_asfall);
                    
                case 'RBF'
                            
                case 'SVR'
                    
            end
            
            
        case 5
            
            switch(upper(opt.current_m5_metamodel{1}))
                
                case 'GP' 
                    
                    [popObj, ~] = predictor(pop(:,1:opt.V), opt.dmodel_asfcv);
                    
                case 'RBF'
                    
                case 'SVR'
                    
            end
            
            
        case 6%---------------model selection function---------------------
            
            switch(upper(opt.current_m6_metamodel{1}))
                
                case 'GP' 
                    
                    [popObj, ~] = predictor(pop(:,1:opt.V), opt.dmodel_asfcvall);
                    popObj = min(popObj, [], 2);
                    
                case 'RBF'
                    
                case 'SVR'
                    
                case 'NN'
                    
                    popObj2 = opt.net(pop(:,1:opt.V)'); %prediction
                    popObj = popObj2';
                    
            end
            
    end
    
    %---------------EVALUATE CONSTRAINTS/AGGREGATION-----------------------
    
    if opt.C>0
        
        switch(opt.methodology)
        
            case {11, 12, 31, 32}
                
                if all(strcmpi(opt.current_cons_metamodel, 'GP'))
                
                    [popCons,~] = predictor(pop, opt.dmodel_cons);
                    
                else
                    popCons = zeros(size(pop,1), opt.C);
                    for i=1:opt.C
                        
                        switch(upper(opt.current_cons_metamodel{i}))
                            
                            case 'GP'
                                
                                [popCons(:,i),~] = predictor(pop, opt.dmodel_cons{i});
                            
                            case 'RBF'
                                
                                
                            
                            case 'SVR'
                        end
                    end
                end
                
                [popCV,~] = evaluateCV(popCons);
                
            case {21, 22, 41, 42}
                
                switch(upper(opt.current_cons_aggregation_option{1}))
                
                    case 'CV'
                        
                        [popCV,~] = predictor(pop(:,1:opt.V), opt.dmodel_cv);
                        [popCV,~] = evaluateCV(popCV);
                        
                    case 'ACV'
                        
                        [popACV,~] = predictor(pop(:,1:opt.V), opt.dmodel_acv);
                        [popCV,~] = evaluateCV(popACV);
                        
                    case 'ELCV'
                        
                        [popELCV,~] = predictor(pop(:,1:opt.V), opt.dmodel_elcv);
                        [popCV,~] = evaluateCV(popELCV);
                        
                    case 'AELCV'
                        [popAELCV,~] = predictor(pop(:,1:opt.V), opt.dmodel_aelcv);
                        [popCV,~] = evaluateCV(popAELCV);
                        
                    case 'TANH'                        
                        [popTHCV,~] = predictor(pop(:,1:opt.V), opt.dmodel_thcv);
                        [popCV,~] = evaluateCV(popTHCV);
                        
                    case 'ATANH'
                        [popATHCV,~] = predictor(pop(:,1:opt.V), opt.dmodel_athcv);
                        [popCV,~] = evaluateCV(popATHCV);
                        
                end
                popCons = zeros(size(pop,1), 1);
                
            case {5,6}
                
                popCons = zeros(size(pop,1), 1);
                popCV = zeros(size(pop,1), 1);
        end
    else
        popCons = zeros(size(pop,1), 1);
        popCV = zeros(size(pop,1), 1);
    end
    
    
    
end