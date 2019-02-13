function opt = methodology32(opt)%generative, one by one
    
    
%     [~,ia,~] = unique(opt.archive,'rows');
%     temp_archive = opt.archive(ia,:);
    opt.normalizedObj = normalize(opt, opt.archiveObj, opt.min_val, opt.max_val);

%     sign  = 1;
%     l = 1;
%     u = opt.numdir;
%     opt.dmodel_asf = cell(opt.numdir,1);
    opt = evaluateASFAll(opt);            
    opt = metamodel(opt, 1:size(opt.archive, 1));%build model
    
%     for dir=l:sign:u   
%         
%         sign = sign*(-1);
%         temp = u;
%         u = l;
%         l = temp;
%         %------------ASF WITH CURRENT DIRECTION----------------------------
%                     
%         opt.curdir = opt.dirs(dir,:);%current direction
%         opt.curcluster = dir;%current cluster number
%         opt.normalizedObj = normalize(opt, opt.archiveObj, opt.min_val, opt.max_val );
%         %opt.normalizedObj =  nsga3_normalization(opt, opt.archiveObj);
%         [opt.archiveASFAll{dir}, opt.archiveCluster] = evaluateASF(opt);%evaluate ASF for all solutions, previous and new
%         
%         %----------------FIND ACTIVE SET-----------------------------------
%         %opt = find_active_set2(opt);
% 
%         %--------------UNIQUES POPULATION----------------------------------
%         %opt = unique_population(opt);
%         
% 
%         
%     end
    
    %temp_archiveASF = cell2mat(opt.archiveASFAll);%opt.archiveASF(ia,:);  
    %temp_archiveASF = temp_archiveASF(ia,:);
    
    %-------------------MODEL ASF FOR ALL DIRECTION------------------------
        
    %[opt.dmodel_asf, ~] = dacefit(temp_archive, temp_archiveASF, @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model asf
    %plot_points(opt);
    
    %---------------------BUILD MODEL FOR CONSTRAINTS OR CV----------------
%     if opt.C>0
%         if opt.methodology==32            
%             [opt.dmodel_cons, ~] = dacefit(temp_archive, opt.archiveCons(ia,:), @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model con
%         elseif opt.methodology==42
%             [opt.dmodel_cv, ~] = dacefit(temp_archive, opt.archiveCVModified(ia, :), @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model cv
%             %[opt.dmodel_cv{dir}, ~] = dacefit(opt.activeArchive, opt.archiveCVModified, @regpoly1, @corrgauss, opt.theta, opt.lob, opt.upb);%model cv
%         end
%     end
    
    %----------------------------------------------------------------------
    %-----------------RUN RGA BASED ON THOSE DIRECTION---------------------
    %----------------------------------------------------------------------
    
    pop = methodology32_rga_metamodel(opt);%it will return one good solution which maximizes current reference direction
    
    if(size(pop,1)+opt.funcEval>opt.totalFuncEval)%discard excess solution in last iteration
        pop = pop(1:opt.totalFuncEval-opt.funcEval,:); 
    end
                             
    %-----------------STORE RESULTS & PLOT-----------------------------
    [popObj, popCons] = evaluate_pop(opt, pop);%high fidelity computation    
    opt = store_results(opt, pop, popObj, popCons); %HI-FI+ASF+KKTPM+CLUSTER
    
    if opt.plotOption==1
        plot_points(opt);  
    end
    
    
end


%------------------------------END OF -FILE--------------------------------