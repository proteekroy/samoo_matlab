function opt = compute_performance_metric(opt, popObj, popCV)

    temp_popObj = popObj(popCV<=0, :);%feasible Paretofront
    
    if size(temp_popObj,1)>0
        if size(temp_popObj,2)==1
            temp_popObj = min(temp_popObj);
        else
            [R,~] = bos(temp_popObj);
            temp_popObj = temp_popObj(R==1, :);
        end
            
        problem = lower(opt.objfunction);
        switch(problem)

            case {'sdtlz1', 'sdtlz2'}
                max_val = max(opt.PF);
                min_val = min(opt.PF);
                temp_popObj = (temp_popObj-min_val)./repmat(max_val-min_val, size(temp_popObj,1),1);
                ParetoFront =  (opt.PF-min_val)./repmat(max_val-min_val, size(opt.PF,1),1);
            otherwise
                ParetoFront = opt.PF;
        end
        opt.igd = IGD(temp_popObj, ParetoFront);
        opt.hv = HV(temp_popObj, ParetoFront);
        opt.gd = GD(temp_popObj, ParetoFront);
    else
        disp('No feasible solution');
        opt.igd = opt.Inf;
        opt.hv = 0;
        opt.gd = opt.Inf;
        disp(['Constraint Violation: ', num2str(min(popCV))])
    end

end