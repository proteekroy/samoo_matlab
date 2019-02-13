
function [opt] = store_results(opt, pop, pop_obj, pop_cons)


    %======================================================================
    %=========EVALUATE OBJ, CV, ACV, CONS,CV-Matrix========================
    %======================================================================
    
    opt.archive = vertcat(opt.archive, pop); %archive population variables  
    opt.archiveObj = vertcat(opt.archiveObj, pop_obj); %archive objective values
    
    if opt.C>0
        [pop_cv, pop_acv] = evaluateCV(pop_cons);%find constraint violation (CV) for new solutions
        opt.archiveCV = vertcat(opt.archiveCV, pop_cv); %archive constrain violation
        opt.archiveACV = vertcat(opt.archiveACV, pop_acv);
        opt.archiveCons = vertcat(opt.archiveCons, pop_cons);%archive constraint values 
        opt.normalizedCons = normalize_cons(opt.archiveCons, opt.archiveCV);
        [opt.normalizedCV, opt.normalizedACV] = evaluateCV(opt.normalizedCons);
        
        opt.archiveCVMatrix = opt.archiveCons;
        if opt.consOption==1 %negative values treated as zero
            opt.archiveCVMatrix(opt.archiveCVMatrix<0)=0;
        end
    else
        pop_cv = zeros(size(pop,1),1);
        opt.archiveCV = vertcat(opt.archiveCV, pop_cv); %archive constrain violation
        opt.archiveACV =  vertcat(opt.archiveACV, pop_cv); 
        opt.archiveCons = vertcat(opt.archiveCons, pop_cons);%archive constraint values  
        opt.normalizedCons = opt.archiveCons;
        opt.normalizedCV  = opt.archiveCV;
        opt.normalizedACV = opt.archiveACV;
    end
    
    %======================================================================
    %=======================WRITE TO FILE==================================
    %======================================================================
    
    if opt.writeFlag==1
        dlmwrite(opt.varfilename, pop, 'delimiter',' ','precision','%.10f','-append');
        dlmwrite(opt.objfilename, pop_obj, 'delimiter',' ','precision','%.10f','-append');
        dlmwrite(opt.cvfilename, pop_cv, 'delimiter', ' ','precision','%.10f','-append');
        dlmwrite(opt.consfilename,  pop_cons, 'delimiter', ' ','precision','%.10f','-append');
    end
    
    %======================================================================
    %====================FIND FEASIBLE PARETO FRONT========================
    %======================================================================
    %--------------------------------------------
    
    if opt.C>0 %Feasible Pareto front of Constraint Problems
        index = find(opt.archiveCV<=0);
    else
        index = (1:size(opt.archive,1))';
    end
    
    opt.FeasibleIndex = index;
    
    if size(index,1)>0 % there are some feasible solutions
        index2 = paretoFront(opt.archiveObj(index,:));
        opt.ParetoIndex = index(index2)';       
        opt.nadir_point = opt.archiveObj(opt.ParetoIndex,:);
        if size(opt.nadir_point,1)>1
            opt.nadir_point =  max(opt.nadir_point);
        end
    else % no feasible solution yet
        opt.ParetoIndex = [];
        opt.nadir_point = max(opt.archiveObj);
    end
        
    %-----------------Find Feasible Solution-------------------------------
    
    if size(opt.FeasibleIndex,1)>0
        opt.ideal_point = opt.archiveObj(opt.FeasibleIndex,:);
        if size(opt.FeasibleIndex,1)>1
            opt.ideal_point = min(opt.archiveObj(opt.FeasibleIndex,:));            
        end
    else
        opt.ideal_point = min(opt.archiveObj);
    end
    
    %---------NORMALIZE AND FIND ASF IN NORMALIZED SPACE-------------------
    opt.normalizedObj = normalize(opt, opt.archiveObj,  opt.min_val, opt.max_val);
    opt.normalizedPop = normalize(opt, opt.archive, opt.bound(1,:), opt.bound(2,:));
    [opt] = evaluateASFAll(opt);
    [opt.archiveASF, opt.archiveCluster] = evaluateASF(opt);%evaluate ASF for all solutions, previous and new
    [opt.archiveSF_VAL] = selectionFunction(opt, pop);
    opt.activeArchive = opt.archive;
    opt.activeArchiveObj = opt.archiveObj;
    opt.activeArchiveASF = opt.archiveASF;
    opt.activeArchiveCons = opt.archiveCons;
    opt.activeArchiveCV = opt.archiveCV;
    opt.activeSetSize = min(opt.funcEval, opt.InitialActiveSetSize);% no need
    %------------STORE UNIQUE RESULTS--------------------------------------
    %opt = unique_population(opt);
    %--------------Find Leader---------------------------------------------
    %opt = find_leader(opt);    
    opt.activeSetSize = min(opt.InitialActiveSetSize, size(opt.activeArchive, 1));
    
end

