%--------------------EVALUATE BY METAMODEL---------------------------------
function [popObj, popCV, popCons] = evaluate_pop_metamodel(opt, pop)

  
    sz = size(pop, 1);
    popObj = [];
    popCV = [];
    
    
    switch(opt.metamodelOption)
        
        
        case 1 %Kriging/ Gaussian Process Modeling
                    
            switch(opt.methodology)  
                
                %--------------METHODOLOGY 1,2-----------------------------
                case {11, 21}
                    
                    [popObj, ~] = predictor(pop(:,1:opt.V), opt.dmodel_obj);
                                        
                    %n = size(popObj, 1);
                    %opt.normalizedObj =  nsga3_normalization(opt, opt.archiveObj);
                    %normalized_obj =  nsga3_normalization(opt, vertcat(popObj, opt.archiveObj(opt.ParetoIndex,:)));
                    %normalized_obj =  nsga3_normalization(opt, popObj);
                    %normalized_obj = normalize(opt, popObj, min(popObj), max(popObj));% normalize
                    normalized_obj = popObj;
                    %------------------------------------------------------
                    
                    popObj = calculate_Asf(normalized_obj, opt.curdir);%using M1-1
                    %popObj = popObj(1:n,:);
                    %[popObj2,~] = predictor(pop(:,1:opt.V), opt.dmodel_asf);
                    %disp(norm(popObj2));
                    
                    if opt.C>0
                        popCons = zeros(size(pop,1), opt.C);
                        if opt.methodology==11
                            [popCons,~] = predictor(pop, opt.dmodel_cons);
                            [popCV,~] = evaluateCV(popCons);
                        else
                            [popCV,~] = predictor(pop(:,1:opt.V), opt.dmodel_cv);  
                            [popCV,~] = evaluateCV(popCV);
                        end
                    else
                        popCV = zeros(sz,1);
                        popCons = zeros(size(pop,1), 1);
                    end          
                    %{
                    popCons = zeros(size(pop,1), opt.C);
                    popObj = zeros(size(pop,1), opt.M);
                    
                    popCons(:,1) = quadratic_predictor(opt.regC{1}, pop(:,1:opt.V));
                    popCons(:,2) = quadratic_predictor(opt.regC{2}, pop(:,1:opt.V));
                    popCons(:,3) = quadratic_predictor(opt.regC{3}, pop(:,1:opt.V));
                    popCons(:,4) = quadratic_predictor(opt.regC{4}, pop(:,1:opt.V));
                    popCons(:,5) = quadratic_predictor(opt.regC{5}, pop(:,1:opt.V));
                    popCons(:,6) = quadratic_predictor(opt.regC{6}, pop(:,1:opt.V));
                    
                    popObj(:,1) = quadratic_predictor(opt.regO{1}, pop(:,1:opt.V));
                    popObj(:,2) = quadratic_predictor(opt.regO{2}, pop(:,1:opt.V));
                    
                    popCV = evaluateCV(opt, popCons);
                    %}
                %--------------METHODOLOGY 3------------------------------- 
                case {12, 22}
                    
                    [popObj, ~] = predictor(pop(:,1:opt.V), opt.dmodel_obj);
                                                            
                    if opt.C>0
                        if opt.methodology==12
                            [popCons,~] = predictor(pop, opt.dmodel_cons);
                            [popCV,~] = evaluateCV(popCons);
                        else
                            popCons = zeros(size(pop,1), opt.C);
                            [popCV,~] = predictor(pop(:,1:opt.V), opt.dmodel_cv);  
                            [popCV,~] = evaluateCV(popCV);
                        end
                    else
                        popCV = zeros(sz,1);
                        popCons = zeros(size(pop,1), 1);
                    end   
                    
                %--------------METHODOLOGY 31, 41--------------------------    
                case {31, 41}
                    
                    [popObj,~] = predictor(pop(:,1:opt.V), opt.dmodel_asf);
                    
                    if opt.C>0
                        
                        if opt.methodology==31
                            [popCons,~] = predictor(pop, opt.dmodel_cons);
                            [popCV,~] = evaluateCV(popCons);
                        else
                            popCons = zeros(size(pop,1), opt.C);
                            [popCV,~] = predictor(pop(:,1:opt.V), opt.dmodel_cv);  
                            [popCV,~] = evaluateCV(popCV);
                        end
                    else
                        popCV = zeros(sz,1);
                        popCons = zeros(size(pop,1), 1);
                    end
                   
%                     popCons(:,1) = quadratic_predictor(opt.regC{1}, pop(:,1:opt.V));
%                     popCons(:,2) = quadratic_predictor(opt.regC{2}, pop(:,1:opt.V));
%                     popCons(:,3) = quadratic_predictor(opt.regC{3}, pop(:,1:opt.V));
%                     popCons(:,4) = quadratic_predictor(opt.regC{4}, pop(:,1:opt.V));
%                     popCons(:,5) = quadratic_predictor(opt.regC{5}, pop(:,1:opt.V));
%                     popCons(:,6) = quadratic_predictor(opt.regC{6}, pop(:,1:opt.V));
%                    
%                     %popObj(:,1) = quadratic_predictor(opt.regO{1}, pop(:,1:opt.V));
%                     %popObj(:,2) = quadratic_predictor(opt.regO{2}, pop(:,1:opt.V));
%
%                     %popObj = calculate_Asf(opt, popObj, opt.curdir);%using M1-1
%                     
%                     popCV = evaluateCV(opt, popCons);
                    
                %--------------METHODOLOGY 32,42---------------------------
                case {32, 42}

                    [popObj,~] = predictor(pop(:,1:opt.V), opt.dmodel_asf);
                    
                    if opt.C>0
                        if opt.methodology==32
                            [popCons,~] = predictor(pop, opt.dmodel_cons);
                            [popCV,~] = evaluateCV(popCons);
                        else
                            popCons = zeros(size(pop,1), opt.C);
                            [popCV,~] = predictor(pop(:,1:opt.V), opt.dmodel_cv);  
                            [popCV,~] = evaluateCV(popCV);
                        end
                    else
                        popCV = zeros(sz,1);
                        popCons = zeros(size(pop,1), 1);
                    end
                %--------------METHODOLOGY 5-------------------------------
                case 5
                    [popObj, ~] = predictor(pop(:,1:opt.V), opt.dmodel_asf);
                    %popObj = exp_imp(popObj, opt.bestASFCV, MSE);
                    popCons = zeros(sz,1);
                    popCV = zeros(sz, 1);
                %--------------METHODOLOGY 6-------------------------------
                case 6
                    
                    %popObj = opt.net(pop(:,1:opt.V)'); %prediction
                    %popObj = popObj';  
                    %---pop_cluster = niching(opt,pop);%decision space association
                    %---[popObj,MSE] = predictor(pop(:,1:opt.V), opt.dmodel_asf);
                    %---popObj = exp_imp(popObj, MSE, opt.bestpopASF, pop_cluster);%make expected improvement as objective
                    %---popObj1 = popObj1 + MSE;
                    %popCV = zeros(sz, 1);
                    popObj2 = opt.net(pop(:,1:opt.V)'); %prediction
                    popObj = popObj2';                  
                    popCons = zeros(size(pop,1), 1);
                    popCV = zeros(size(pop,1), 1);
%                   popObj = zeros(size(pop,1), opt.M);

                    
                %--------------METHODOLOGY 7-------------------------------   
                case 7
                    popCV = zeros(sz, 1);
                    [popObj,~] = predictor(pop(:,1:opt.V), opt.dmodel_kktpm);
                    popCons = zeros(sz,1);                
                    
                otherwise
                    
                    disp('Inside evaliMetamodel evaluation is not implemented');
                    
            end
            
        case 2 %Quadradic 
            
        case 3 %HDMR
            
        case 4 %RBF
    
        otherwise
            
            
    end
 
end


function asf = calculate_Asf(obj, dir)%utopian)

    %w = (dir-utopian)./norm(dir-utopian);%unit direction
    %w = dir;
    %asf  = max((obj-repmat(opt.min_val,size(obj,1),1))./repmat(dir,size(obj,1),1),[],2);
    
    %dir = dir./norm(dir);
    %asf = max(obj.*dir);
    %asf  = max(obj./dir); % max(w.*(obj-utopian));  
    
    %No Normalization is done
    %asf = max(obj - repmat(dir-1/size(obj,2),size(obj,1),1),[],2);%+0.01*sum(obj-(dir-0.5));%+0.0001*sum(obj.*dir);
    asf = max( (obj-repmat(dir, size(obj,1),1)),[],2);
    
end

function cv = evaluateCVhere(opt, pop_cons)
    
    g = pop_cons;
    if opt.consOption==1 %negative values treated as zero
        for i=1:size(pop_cons,2)
           g(g(:,i)<0,i)=0;
        end
    end  

    cv = abs(g);
    %cv = cv./max(cv);
    cv = cv./opt.maxCV;
    cv = max(cv,[],2);
end


