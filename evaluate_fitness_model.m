function pop = evaluate_fitness_model(opt, pop)

    % rGA PREDICTOR FUNCTION EVALUATIONS
    for i=1:size(pop,1)
        
        switch(opt.methodology)
            case 1
                disp('Inside evaluate fitness model, option 1 not implemented');
            case 2
                disp('Inside evaluate fitness model, option 2 not implemented');
            case 3
                [yhat,~,~,~] = predictor(pop(i,1:opt.V), opt.dmodel_asf);
                pop(i, opt.V+1) = yhat;%predicted asf
                if opt.C>0
                    
                    for j=1:opt.C
                        [cv(j),~,mse_cv(j),~] = predictor(pop(i,1:opt.V), opt.dmodel_cv(j));
                    end
                    I = cv<=0;
                    cv(I)=0;
                    %ei = exp_imp(cv,mse_cv,opt.ybest_cv);
                    total_cv = sum(cv);
                    if(total_cv>0)
                        pop(i, opt.V+2) = total_cv;%<= constrained
                    else
                        pop(i, opt.V+2) = 0;%predicted cv
                    end
                else
                    pop(i, opt.V+2) = 0;
                end
                %ei = exp_imp(yhat,mse_asf,opt.ybest_asf);
                %constrained_ego(opt, opt.ybest, yhat, cv, mse_cv)
               
                
            case 4
                 %%{
                [yhat,~,mse_asf,~] = predictor(pop(i,1:opt.V), opt.dmodel_asf);
                pop(i, opt.V+1) = yhat;%predicted asf
                %ei = exp_imp(yhat,mse_asf,opt.ybest_asf);
                %constrained_ego(opt, opt.ybest, yhat, cv, mse_cv);
                
                if opt.C>0
                    [cv,~,mse_cv,~] = predictor(pop(i,1:opt.V), opt.dmodel_cv);
                    %ei = exp_imp(cv,mse_cv,opt.ybest_cv);
                    if(cv>0)
                        pop(i, opt.V+2) = cv;%<= constrained
                    else
                        pop(i, opt.V+2) = 0;%predicted cv
                    end
                else
                    pop(i, opt.V+2) = 0;
                end
                %}
                %{ 
                when using actual fitness values
                [f,g] = high_fidelity_calculation(opt,pop(i,:));
                
                
                normalizedobj = (f-opt.min_val)./((opt.max_val+000001)- opt.min_val);      
                asf = evaluateASF(opt,normalizedobj);
                pop(i, opt.V+1) = asf;
                if opt.C>0
                    cv = evaluateCV(opt, pop(i, :), g);
                    pop(i, opt.V+2) = cv;
                else
                    pop(i, opt.V+2) = 0;
                end
                %}
            case 5
                disp('Inside evaluate fitness model, option 5 not implemented');
            case 6
                %%{
                [yhat,~,mse_kktpm,~] = predictor(pop(i,1:opt.V), opt.dmodel_kktpm);
                %index = dsearchn(opt.activearchive,pop(i,1:opt.V));
                %c = opt.archiveCluster(index);
                %ybest = opt.ybest(c);               
                %ei = exp_imp(yhat,mse_kktpm,ybest);
                pop(i, opt.V+1) = yhat;%ei;%yhat;%ei;%;yhat;%predicted kktpm value
                pop(i, opt.V+2) = 0;%no cv information
                %}
                %yreal = KKT(opt,pop(i,1:opt.V));
                %pop(i, opt.V+1) = yreal;
                %pop(i, opt.V+2) = 0;
            case 7
                
                %[yhat,~,mse_asf,~] = predictor(pop(i,1:opt.V), opt.dmodel_asf);
                [popObj, popCons] = evaluate_pop(opt, pop(i,1:opt.V));
                popObj = normalize(opt, popObj);

                [yhat, ~] = minimumASF(opt, popObj);
                pop(i, opt.V+1) = yhat;%ei;%yhat;%ei;%;yhat;%predicted kktpm value
                
                if opt.C>0
                    [cv,~,mse_cv,~] = predictor(pop(i,1:opt.V), opt.dmodel_cv);
                    %ei = exp_imp(cv,mse_cv,opt.ybest_cv);
                    if(cv>0)
                        pop(i, opt.V+2) = cv;%<= constrained
                    else
                        pop(i, opt.V+2) = 0;%predicted cv
                    end
                else
                    pop(i, opt.V+2) = 0;
                end
                
            otherwise
                disp('Inside evaluate fitness model, option >6 not implemented');
        end

    end   
end
