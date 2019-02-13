function opt = metamodel(opt, TrainIndex)

    x = opt.archive(TrainIndex,:);
    f = opt.archiveObj(TrainIndex, :);
    normalized_f = normalize(opt, opt.archiveObj(TrainIndex,:), opt.min_val, opt.max_val);
    %normalized_f = opt.normalizedObj(TrainIndex, :);
    cons = opt.archiveCons(TrainIndex,:);
    normalized_cons = opt.normalizedCons(TrainIndex,:);
    cv = opt.archiveCV(TrainIndex,:); %archive constrain violation
    acv = opt.archiveACV(TrainIndex,:);
    normalized_cv = opt.normalizedCV(TrainIndex,:);
    normalized_acv = opt.normalizedACV(TrainIndex,:);
    y_asf = cell2mat(opt.archiveASFAll);
    y_asfcv = cell2mat(opt.archiveASFCVAll);
    y_asf = y_asf(TrainIndex,:);
    y_asfcv = y_asfcv(TrainIndex,:);
    
    [~,ia,~] = unique(x,'rows');
    x = x(ia,:);
    cons = cons(ia,:);
    normalized_cons = normalized_cons(ia,:);
    f = f(ia,:);
    
    normalized_f = normalized_f(ia, :);
    cv = cv(ia,:);
    acv = acv(ia,:);
    normalized_cv = normalized_cv(ia, :);
    normalized_acv = normalized_acv(ia, :);
    y_asf = y_asf(ia,:);
    y_asfcv = y_asfcv(ia,:);
    
    
    %----------------MODEL OBJECTIVES/AGGREGATION--------------------------
    
    switch(opt.methodology)
        
        case {11, 12, 21, 22} %objectives are separately metamodeled
            
            if all(strcmpi(opt.current_obj_metamodel, 'GP'))
                
                if opt.methodology==11 || opt.methodology==21
                    [opt.dmodel_obj, ~] = dacefit(x, f, @regpoly1, @corrgauss, opt.theta, opt.lob, opt.upb);%model all objectives with GP
                else
                    [opt.dmodel_obj, ~] = dacefit(x, f, @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);
                end
                %[opt.dmodel_obj, ~] = dacefit(x, f, @regpoly1, @corrgauss, opt.theta, opt.lob, opt.upb);%model obj
                %[opt.dmodel_obj, ~] = dacefit(x, normalized_f, @regpoly1, @corrgauss, opt.theta, opt.lob, opt.upb);%model obj
                
            else
                for i=1:opt.M
                    
                    switch(upper(opt.current_obj_metamodel{i}))
                        
                        case 'GP' 
                            
                            [opt.dmodel_obj{i}, ~] = dacefit(x, f(:,i), @regpoly1, @corrgauss, opt.theta, opt.lob, opt.upb);%model each objective
                            
                        case 'RBF'
                            
                        case 'SVR'
                    end
                end
            end
            
        case {31, 41} %%aggregated function of objectives is metamodeled
            
            switch(upper(opt.current_obj_aggregation_option{1}))
                
                case 'ASF'
                    
                    [opt.dmodel_asf, ~] = dacefit(x, y_asf(:, opt.curcluster), @regpoly1, @corrgauss, opt.theta, opt.lob, opt.upb);%model asf for current ref point
                    
                case 'WS'
                    
                case 'TCHEBYCHEFF'
                    
            end
            
        case {32, 42}
            
            switch(upper(opt.current_obj_aggregation_option{1}))
                
                case 'ASF'
                    
                    [opt.dmodel_asfall, ~] = dacefit(x, y_asf, @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model asf for all ref points
                    
                case 'WS'
                    
                case 'TCHEBYCHEFF'
                    
            end
            
        case 5
            
            switch(upper(opt.current_selection_function_option_m5{1}))
                
                case 'ASFCV' %asf+cv for specific ref points
                    [opt.dmodel_asfcv, ~] = dacefit(x, y_asfcv(:, opt.curcluster), @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model asfcv for current ref point
                    
            end
        
        case 6%---------------model selection function---------------------
            
            switch(upper(opt.current_selection_function_option_m6{1}))
                
                case 'MEMO' %minimum of all asf+cv for different ref points
                
                    SF_VAL = opt.archiveSF_VAL(TrainIndex,:);
                    SF_VAL = SF_VAL(ia,:);
                    opt.net = train_neural_network(x, SF_VAL);
                    
                case 'ASFCV'%asf+cv for each ref points
                    
                    [opt.dmodel_asfcvall, ~] = dacefit(x, y_asfcv, @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model asfcv for all ref points
                    
            end
    end
    
    %----------------MODEL CONSTRAINTS/AGGREGATION-------------------------
    
    if opt.C>0
        
        switch(opt.methodology)
        
            case {11, 12, 31, 32}
                %feasible_index = find(cv==0);
                %g2 = sum(cons, 2);
                
                if all(strcmpi(opt.current_cons_metamodel, 'GP'))
                
                    [opt.dmodel_cons, ~] = dacefit(x, cons, @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model each constraint separately
                    %[opt.dmodel_cons, ~] = dacefit(x, normalized_cons, @regpoly1, @corrgauss, opt.theta, opt.lob, opt.upb);%model each constraint separately
                    %{
                    elcv = normalized_cons;
                    alpha = 0.5;
                    for i=1:opt.C
                        elcv(feasible_index,i) = alpha.*(exp(g2(feasible_index))-1); 
                    end
                    [opt.dmodel_cons, ~] = dacefit(x, elcv, @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model each constraint separately
                    %}
                    
                
                else
                    for i=1:opt.C
                        switch(upper(opt.current_cons_metamodel{i}))
                            
                            case 'GP'
                                [opt.dmodel_cons{i}, ~] = dacefit(x, cons(:,i), @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model each constraint separately
                            case 'RBF'
                            
                            case 'SVR'
                        end
                    end
                end
                
            case {21, 22, 41, 42}
                
                %training set
                feasible_index = find(cv==0);
                g2 = sum(cons, 2);
                
                switch(upper(opt.current_cons_aggregation_option{1}))
                
                    case 'CV'
                        
                        %[opt.dmodel_cv, ~] = dacefit(x, cv, @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model cv function
                        [opt.dmodel_cv, ~] = dacefit(x, normalized_cv, @regpoly1, @corrgauss, opt.theta, opt.lob, opt.upb);%model normalized cv function
                        
                        
                    case 'ACV'
                        
                        %[opt.dmodel_acv, ~] = dacefit(x, acv, @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model acv function
                        [opt.dmodel_cv, ~] = dacefit(x, normalized_acv, @regpoly1, @corrgauss, opt.theta, opt.lob, opt.upb);%model normalized cv function
                        
                    case 'ELCV'
                        elcv = cons;
                        alpha = 0.5;
                        for i=1:opt.C
                            elcv(feasible_index,i) = alpha.*(exp(g2(feasible_index))-1); 
                        end
                        [opt.dmodel_elcv, ~] = dacefit(x, elcv, @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model aelcv function
                                            
                    case 'AELCV'
                        aelcv = cv;
                        alpha = 0.5;
                        aelcv(feasible_index) = alpha*(exp(g2(feasible_index))-1);
                        
                        [opt.dmodel_aelcv, ~] = dacefit(x, aelcv, @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model aelcv function

                    case 'TANH'
                        thcv = tanh(cons);
                        [opt.dmodel_thcv, ~] = dacefit(x, thcv, @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model thcv function
                                                
                    case 'ATANH'
                        athcv = tanh(cv_modified);
                        [opt.dmodel_athcv, ~] = dacefit(x, athcv, @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model athcv function
                        
                        
                end               
        
        end
    end
    
    
    
end