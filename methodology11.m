function opt = methodology11(opt)%generative, one by one

    repeated_times = 1;
    sign  = -1;
    l = 21;
    u = 1;%opt.numdir;
    
    %for dir=1:opt.numdir%:-1:1
    for dir = 1:1: size(opt.dirs,1) 
    %for dir=l:sign:u   
        
        sign = sign*(-1);
        temp = u;
        u = l;
        l = temp;
        %------------ASF WITH CURRENT DIRECTION----------------------------
                    
        opt.curdir = opt.dirs(dir,:);%current direction
        opt.curcluster = dir;%current cluster number
        
        for repeat = 1:repeated_times
            
            opt = evaluateASFAll(opt);            
            opt = metamodel(opt, 1:size(opt.archive, 1));%build model
            
            %--------------------RGA MODEL---------------------------------
            [pop, opt] = methodology3_rga_metamodel2(opt);
            
            %-----------------STORE RESULTS & PLOT-----------------------------
            [popObj, popCons] = evaluate_pop(opt, pop);%high fidelity computation
            opt = store_results(opt, pop, popObj, popCons);%HI-FI+ASF+KKTPM+CLUSTER

            %----------------APPLY TRUST REGION----------------------------
            
            %-------------FUNCTION EVALUATION CHECK------------------------
            opt.funcEval = size(opt.archive,1);%number of function evaluations
            if opt.funcEval>=opt.totalFuncEval
                break;
            end
        end
        
    end    
    
end

%------------------------------END OF -FILE--------------------------------
