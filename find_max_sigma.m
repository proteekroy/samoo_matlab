function find_max_sigma(opt)

    
    pop = lhsamp_model(opt.N, opt);%LHC initialize
    
    
    popObj = zeros(size(pop,1), opt.M);
    for i=1:opt.M
        [~, popObj] = predictor(pop(:,1:opt.V), opt.dmodel_obj{i});
    end
                           
    
    if opt.C>0
        sz = size(pop, 1);
        popCons = zeros(sz, opt.C);
        popConsigma = zeros(sz, opt.C);
        popConsCDF = zeros(sz, opt.C);
        mu = zeros(sz, 1);
        for i=1:opt.C
            [popCons(:, i), popConsigma(:, i)] = predictor(pop, opt.dmodel_cons{i});
            popConsCDF(:, i) = normcdf(popCons(:, i), mu , popConsigma(:, i));
        end
        
        %popCV = evaluateCV(opt, popCons);
    else
        popCV = zeros(sz,1);
        popCons = zeros(sz, 1);
    end
    
    index = find(popCV<=0);
    temp = max(popObj);
    
    
    [popObj, popCV, popCons] = evaluate_pop_metamodel(opt, pop);%objective value prediction   
    
    %tic
    %---------SIMPLE REAL PARAMETER GENETIC ALGORITHM----------------------
    opt.gen = 1;
    totalpopObj = zeros(2*opt.N,1);
    totalpop = zeros(2*opt.N,opt.V);
    totalpopCV = zeros(2*opt.N,1);
    
    %tic
    while opt.gen <= opt.G

        %------------Mating Selection and Child Creation-------------------
        popChild = constrained_tournament_selection(opt, pop, popObj, popCV);
        [popChild, opt.nrealcross] = sbx( popChild, opt.pcross, opt.nrealcross, opt.eta_c, opt.bound(1,:), opt.bound(2,:));%simulated binary crossover
        [popChild, opt.nrealmut] = pol_mut( popChild, opt.pmut, opt.nrealmut, opt.eta_m, opt.bound(1,:), opt.bound(2,:));%polynomial mutation

        %--------------Evaluate--------------------------------------------
        [popChildObj, popChildCV, popChildCons] = evaluate_pop_metamodel(opt, popChild);%objective value prediction
        %popChildCV = evaluateCV(opt, popChildCons);

        %--------------Merge-----------------------------------------------
        %totalpopObj = vertcat(popChildObj, popObj);
        %totalpop = vertcat(popChild, pop);
        %totalpopCV = vertcat(popChildCV, popCV);
        totalpopObj(1:opt.N,:) = popChildObj;
        totalpopObj(opt.N+1:end,:) = popObj;
        totalpop(1:opt.N,:) = popChild;
        totalpop(opt.N+1:end,:) = pop;
        totalpopCV(1:opt.N,:) = popChildCV;
        totalpopCV(opt.N+1:end,:) = popCV;
        
        %---------------Survival Selection---------------------------------
        %if opt.phase==1
        %[pop, popObj, popCV]  = survival_tournament_selection(opt, totalpop, totalpopObj, totalpopCV);
        %    [pop, popObj, popCV]  = trustedObj_survival_selection2(opt, totalpop, totalpopObj, totalpopCV, opt.N);
        %else
        [pop, popObj, popCV]  = trustedObj_survival_selection(opt, totalpop, totalpopObj, totalpopCV);
        %end
        opt.gen = opt.gen+1;

    end

    %t = toc;
    %disp(t);
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



end