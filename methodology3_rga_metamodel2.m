function [pop, opt] = methodology3_rga_metamodel2(opt)

    pop = lhsamp_model(opt.N, opt);%LHC initialize
    [popObj, popCV, popCons] = evaluate_model_space(opt, pop);%objective value prediction
    
    %---------SIMPLE REAL PARAMETER GENETIC ALGORITHM----------------------
    opt.gen = 1;
    totalpopObj = zeros(2*opt.N,1);
    totalpop = zeros(2*opt.N,opt.V);
    totalpopCV = zeros(2*opt.N,1);
    
    while opt.gen <= opt.G

        %------------Mating Selection and Child Creation-------------------
        popChild = constrained_tournament_selection(opt, pop, popObj, popCV);
        [popChild, opt.nrealcross] = sbx( popChild, opt.pcross, opt.nrealcross, opt.eta_c, opt.bound(1,:), opt.bound(2,:));%simulated binary crossover
        [popChild, opt.nrealmut] = pol_mut( popChild, opt.pmut, opt.nrealmut, opt.eta_m, opt.bound(1,:), opt.bound(2,:));%polynomial mutation

        %--------------Evaluate--------------------------------------------
        [popChildObj, popChildCV, popChildCons] = evaluate_model_space(opt, popChild);%objective value prediction
        %popChildCV = evaluateCV(opt, popChildCons);

        %--------------Merge-----------------------------------------------
        totalpopObj(1:opt.N,:) = popChildObj;
        totalpopObj(opt.N+1:end,:) = popObj;
        totalpop(1:opt.N,:) = popChild;
        totalpop(opt.N+1:end,:) = pop;
        totalpopCV(1:opt.N,:) = popChildCV;
        totalpopCV(opt.N+1:end,:) = popCV;
        
        %---------------Survival Selection---------------------------------
        [pop, popObj, popCV]  = trustedObj_survival_selection(opt, totalpop, totalpopObj, totalpopCV);
        opt.gen = opt.gen+1;

    end
    %-----------UNIQUE SOLUTION--------------------------------------------
    [~,ia,~] = unique(pop,'rows');%find unique points
    pop = pop(ia,:);
    popObj = popObj(ia,:);
    popCV = popCV(ia,:);
    
    %------------RETURN ONE SOLUTION---------------------------------------
    index = find(popCV<=0);
    if isempty(index)%all are infeasible
        [~,index] = min(popCV);
    else %more than one feasible
        [~,index2] = sort(popObj(index,:));
        index = index(index2(1));
    end
    
    pop = pop(index,1:opt.V);
    
    %{
    opt.temp.pop = pop;
    
    if opt.methodology == 11 || opt.methodology == 21
        [opt.temp.popObj, ~] = predictor(pop, opt.dmodel_obj);
        opt.temp.popASF = calculate_Asf(opt.temp.popObj, opt.curdir);%using M1-1
        opt.temp.popObj = opt.temp.popObj';
    else
        [opt.temp.popObj, ~] = predictor(pop, opt.dmodel_asfcv);
        opt.temp.popASF = opt.temp.popObj;%using M1-1 
    end
    
    if opt.C>0
        opt.temp.popCons = zeros(size(pop,1), opt.C);
        if opt.methodology==11 || opt.methodology==31
            [opt.temp.popCons,~] = predictor(pop, opt.dmodel_cons);
            [opt.temp.popCV,~] = evaluateCV(popCons);
        elseif opt.methodology==5
            opt.temp.popCV = zeros(size(pop,1),1);
            opt.temp.popCons = zeros(size(pop,1), 1);
        else
            [opt.temp.popCV,~] = predictor(pop(:,1:opt.V), opt.dmodel_cv);  
            [opt.temp.popCV,~] = evaluateCV(popCV);
        end
    else
        opt.temp.popCV = zeros(size(pop,1),1);
        opt.temp.popCons = zeros(size(pop,1), 1);
    end    
    %}
end


%------------------------------END OF -FILE--------------------------------
