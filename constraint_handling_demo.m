function constraint_handling_demo()


    close all;
    test_function = ['tnk    '; 'osy    '; 'srn    ';'bnh    ';];%test functions
    opt.r = 1;%run number
    
    for func_no = 1:1
        
        opt.objfunction = strtrim(test_function(func_no,:));%remove whitespaces
        opt.func_no = func_no;%function number
        opt.methodology = 11; %methodology
        opt = load_parameters(opt);%parameters, 2obj 2var problem       
        if opt.func_no==2 || opt.func_no==3 || opt.func_no==4 
            opt.pareto_set = zeros(200, 2);
            opt.pareto_set(:,1) = linspace(0,1,200);
        elseif opt.func_no==1
            opt.pareto_set = opt.paretofrontGD;
        elseif opt.func_no==5
            
        end
        
        %-----------------create dataset-----------------------------------
        opt.initpopsize = 200;%LHS sampling
        [pop, popObj, popCons] = initialization(opt);%LHS sampling  
        opt = store_results_for_plot(opt, pop, popObj, popCons);%save results,asf,cv
        TrainIndex = 1:opt.initpopsize;
        opt.TrainIndex = TrainIndex;
        N = 50;%120;%25;
        [test_x, test_obj, test_cons] = generate_meshgrid_data(opt, N);
        opt = store_results_for_plot(opt, test_x, test_obj, test_cons);%save results,asf,cv
        index = randperm(size(opt.archive,1));
        opt.TrainIndex = index(1:30);
        opt.TestIndex = opt.initpopsize+1:size(opt.archive,1);
        %--------------create metamodels-----------------------------------
        opt = build_models(opt, opt.TrainIndex);
        
        %--------------compute selection values----------------------------
        opt = calculate_prediction(opt, opt.TestIndex);
        
        %--------------compute performance metric--------------------------
        opt = calculate_performance_metric(opt);
        
        %-------------plot model space for constraints---------------------
        opt = plot_model_function(opt, N);
        
    end
    
    disp('Done');


end

function opt = plot_model_function(opt, N)

    xtest = opt.archive(opt.TestIndex,:);
    cons_test = opt.archiveCons(opt.TestIndex,:);
    
    %{
    x1 = reshape(xtest(:, 1), N, N);
    x2 = reshape(xtest(:, 2), N, N);
    figure;
    hold all;
    surfc(x1, x2, reshape(opt.aelcv_pred, N, N));
    %}
    
    %[test_x, test_obj, test_cons] = generate_linearized_data(opt, 100, 1);
    
    index = xtest(:,2)==xtest(10, 2);
    %test_x = xtest(index,:);
    %[x, index] = sort(test_x(:,2));
    str = {'Constraint Function', 'CV    ','Aggregated CV      ','Hyperbolic Tangent','Exponential-Linear', 'Aggregated Hyperbolic Tangent', 'Aggregated Exponential-Linear'};
    for i=1:size(opt.pred_func,2)
        if size(opt.pred_func{i},2)==2
            temp = opt.pred_func{i}(index,1);
            temp2  = opt.actual_func{i}(index,1);
        else
            temp = opt.pred_func{i}(index);
            temp2  = opt.actual_func{i}(index);
        end
        figure;
        hold all;
        plot(xtest(index,1), temp, '--db','MarkerSize',7, 'MarkerFaceColor','b','LineWidth',3);
        plot(xtest(index,1), temp2, '-og','MarkerSize',4, 'MarkerFaceColor','g','LineWidth',3);
        plot([0,3.5],[0,0],'r-','LineWidth',3);
        legend(strcat('predicted:',num2str(opt.sep{i})), 'original');
        title(upper(str{i}));
        xlabel('x1');
        ylabel('Function Value');
        %surfc(x1, x2, elcv_pred);
    
    end
    
    %{
    for i=1:size(opt.cv_func,2)
        figure;
        hold all;
        temp = reshape(opt.cv_func{i}, N, N);
        surfc(x1, x2, temp);
        %h = plot((x1, x3opt.cv_func{i}, pareto_front(:, 2), 'ro','MarkerSize',4, 'MarkerFaceColor','r','LineWidth',3);
        %surfc(x1, x2, elcv_pred);
    
    end
    %}
    %{
    c = size();
    
    
    N = sqrt(size(test_obj, 1));
    F1 = reshape(test_obj(:, 1), N, N);
    F2 = reshape(test_obj(:, 2), N, N);
    figure;
    hold all;
    surfc(F1,F2, Z_rank_mesh);
%     plot3(test_obj(Z_rank_test==1, 1), test_obj(Z_rank_test==1, 2),Z_rank_test(Z_rank_test==1), 'bo','MarkerSize',10, 'MarkerFaceColor','b');
%     plot3(test_obj(Z_rank_test==2, 1), test_obj(Z_rank_test==2, 2),Z_rank_test(Z_rank_test==2), 'go','MarkerSize',10, 'MarkerFaceColor','g');
%     plot3(test_obj(Z_rank_test==3, 1), test_obj(Z_rank_test==3, 2),Z_rank_test(Z_rank_test==3), 'co','MarkerSize',10, 'MarkerFaceColor','c');
%     
    %plot3(pareto_front(:, 1), pareto_front(:, 2), Z_rank(I1),'ro','MarkerSize',5, 'MarkerFaceColor','r','LineWidth',4);
    %h = plot(pareto_front(:, 1), pareto_front(:, 2), 'ro','MarkerSize',4, 'MarkerFaceColor','r','LineWidth',3);
    h = plot(pareto_front(:, 1), pareto_front(:, 2), 'ro','MarkerSize',4, 'MarkerFaceColor','r','LineWidth',3);
    %h = plot3(pareto_front(:, 1), pareto_front(:, 2), Z_rank_test(I1), 'ro','MarkerSize',4, 'MarkerFaceColor','r','LineWidth',3);
    %plot3(train_obj(:, 1), train_obj(:, 2), Z_rank_test(I2), 'mo','MarkerSize',5, 'MarkerFaceColor','m','LineWidth',5);
    title(strcat(title_string,'-Exact Space'));
    %title("NSGA-II Exact Space")
    legend(h, 'Pareto-front');
    xlabel('f_1');
    ylabel('f_2');
    %xlim([0 3]);
    %ylim([0 3]);
    zlabel('Selection(x)');
    set(gca,'fontsize',fontsize);
    box on;
    grid on;
    drawnow;
    %}

end


function [opt] = calculate_performance_metric(opt)


    opt.Z_rank =  cell(1, 5);
    opt.Z_rank_pred = cell(1, 5);
    opt.sep = cell(1, 5);
    opt.mse = cell(1, 5);
    
    dummy_zero = zeros(size(opt.cv_actual,1),1);
    
    %------------------CALCULATE SEP ERROR---------------------------------
    
    [opt.Z_rank{1}, opt.Z_rank_pred{1}, opt.sep{1}] = compute_sep_error(opt.cons_actual, opt.cons_pred, dummy_zero, dummy_zero);
    
    [opt.Z_rank{2}, opt.Z_rank_pred{2}, opt.sep{2}] = compute_sep_error(opt.cv_actual, opt.cv_pred, dummy_zero, dummy_zero);
    
    [opt.Z_rank{3}, opt.Z_rank_pred{3}, opt.sep{3}] = compute_sep_error(opt.acv_actual, opt.acv_pred, dummy_zero, dummy_zero);
    
    [opt.Z_rank{4}, opt.Z_rank_pred{4}, opt.sep{4}] = compute_sep_error(opt.thcv_actual, opt.thcv_pred, dummy_zero, dummy_zero);
    
    [opt.Z_rank{5}, opt.Z_rank_pred{5}, opt.sep{5}] = compute_sep_error(opt.elcv_actual, opt.elcv_pred, dummy_zero, dummy_zero);     
    
    [opt.Z_rank{6}, opt.Z_rank_pred{6}, opt.sep{6}] = compute_sep_error(opt.athcv_actual, opt.athcv_pred, dummy_zero, dummy_zero);
    
    [opt.Z_rank{7}, opt.Z_rank_pred{7}, opt.sep{7}] = compute_sep_error(opt.aelcv_actual, opt.aelcv_pred, dummy_zero, dummy_zero);
    
    
    str = {'cons  ', 'cv    ','acv      ','thcv  ','elcv  ', 'athcv ', 'aelcv '};
    disp('SEP ERROR');
    disp('=====================');
    for i=1:7
        disp([strtrim(str{i}) ': ' num2str(opt.sep{i})]);
    end
    
    %------------------CALCULATE MSE ERROR---------------------------------
    
    opt.mse{1} = compute_mse_error(opt.cons_actual, opt.cons_pred);
    
    opt.mse{2} = compute_mse_error(opt.cv_actual, opt.cv_pred);
    
    opt.mse{3} = compute_mse_error(opt.acv_actual, opt.acv_pred);
    
    opt.mse{4} = compute_mse_error(opt.thcv_actual, opt.thcv_pred);
    
    opt.mse{5} = compute_mse_error(opt.elcv_actual, opt.elcv_pred);     
    
    opt.mse{6} = compute_mse_error(opt.athcv_actual, opt.athcv_pred);
    
    opt.mse{7} = compute_mse_error(opt.aelcv_actual, opt.aelcv_pred);
    
    disp('MSE ERROR');
    disp('=====================');
    for i=1:7
        disp([strtrim(str{i}) ': ' num2str(opt.mse{i})]);
    end

end


function mse = compute_mse_error(actual_value, predicted_value)

    n = size(actual_value,1);
%     d = vertcat(actual_value, predicted_value);
%     a = (actual_value-repmat(min(d),n,1))./((repmat(max(d),n,1)+1e-16)-repmat(min(d),n,1));
%     b = (predicted_value-repmat(min(d),n,1))./((repmat(max(d),n,1)+1e-16)-repmat(min(d),n,1)); 
%     mse = mean(vecnorm(abs(a-b),2,2));
    %mse = mean(vecnorm(abs(actual_value-predicted_value)./(actual_value+1e-16), 2, 2));
    mse = mean(vecnorm(abs((actual_value-predicted_value)./(repmat(max(actual_value),n,1)-repmat(min(actual_value),n, 1))), 2, 2));
end

function [d] = constrained_domination(obj1, obj2, cv1, cv2)
    
    if (cv1<=0 && cv2<=0) || (cv1-cv2)<1e-16 %both are feasible or same CV

        dom = lex_dominate(obj1, obj2);
        
        if dom == 1 %p1 dominates p2
            d = 1;%1 dominates 2
        elseif dom == 3 % p2 dominates p1
            d = 2;%2 dominates 1
        else
            d = 3;%both nondominated
        end
        %}
    else
        if cv1 < cv2%p1 less constraint violation
            d=1;%1 dominates 2
        else
            d=2;%2 dominates 1 
        end
    end
end


%This function computes selection function values for frameworks
function [Z_rank, Z_rank_pred, sep] = compute_sep_error(actual_obj, predicted_obj, actual_value, predicted_value)

    n = size(actual_value,1);
    Z_rank = zeros(n, 1);
    Z_rank_pred = zeros(n, 1);

    %---------------------COMPUTE SELECTION FUNCTION-----------------------
    selection_err = 0;
    count = 0;
    
    for i=1:n-1
        for j=i+1:n
            d1 = constrained_domination(actual_obj(i,:), actual_obj(j,:), actual_value(i,:), actual_value(j,:));%1, 2, 3
            d2 = constrained_domination(predicted_obj(i,:), predicted_obj(j,:), predicted_value(i,:), predicted_value(j,:));
            
            if d1==1
                Z_rank(j) = Z_rank(j)+1;
            elseif d1==2
                Z_rank(i) = Z_rank(i)+1;
            end

            if d2==1
                Z_rank_pred(j) = Z_rank_pred(j)+1;
            elseif d2==2
                Z_rank_pred(i) = Z_rank_pred(i)+1;
            end
            if d1~=d2
                selection_err = selection_err + 1;
            end
            count = count+1;
        end
    end
    
    sep = selection_err/count;

end


%==========================================================================
%=======================SELECTION VALUE====================================
%==========================================================================


function [opt] = calculate_prediction(opt, TestIndex)
    
    xtest = opt.archive(TestIndex,:);
    %test_obj = opt.archiveObj(TestIndex,:);
    cons_test = opt.archiveCons(TestIndex,:);
    
    %Predicted Functions
    opt.cons_pred = predictor(xtest, opt.dmodel_cons);%predict constraints
    opt.cv_pred = predictor(xtest, opt.dmodel_cv);%predict constraint violation
    opt.acv_pred = predictor(xtest, opt.dmodel_acv);%predict aggregated cv
    opt.elcv_pred = predictor(xtest, opt.dmodel_elcv);%predict exp-linear cv
    opt.thcv_pred = predictor(xtest, opt.dmodel_thcv);%predict hyperbolic tangent cv
    opt.aelcv_pred = predictor(xtest, opt.dmodel_aelcv);%predict exp-linear cv
    opt.athcv_pred = predictor(xtest, opt.dmodel_athcv);%predict hyperbolic tangent cv
    
    opt.pred_func = cell(1,7);
    opt.pred_func{1} = opt.cons_pred;
    opt.pred_func{2} = opt.cv_pred;
    opt.pred_func{3} = opt.acv_pred;
    opt.pred_func{4} = opt.elcv_pred;
    opt.pred_func{5} = opt.thcv_pred;
    opt.pred_func{6} = opt.aelcv_pred;
    opt.pred_func{7} = opt.athcv_pred;
    
    %Actual Functions
    g = cons_test;
    opt.cons_actual = cons_test;
    g(g<0) = 0;
    opt.cv_actual = sum(g, 2);
    
    %aggregated cv
    opt.acv_actual = opt.cv_actual;
    feasible_index = find(opt.cv_actual==0);
    g2 = sum(cons_test, 2);
    opt.acv_actual(feasible_index) = g2(feasible_index);
    
    %aggregated hyperbolic tangent
    opt.athcv_actual = tanh(opt.acv_actual);
    
    %hyperbolic tangent
    opt.thcv_actual = tanh(cons_test);
    
    %aggregated exponential linear cv
    opt.aelcv_actual = opt.cv_actual;
    alpha = 0.5;
    opt.aelcv_actual(feasible_index) = alpha*(exp(g2(feasible_index))-1);
    
    
    %exponential linear cv
    opt.elcv_actual = cons_test;
    alpha = 0.5;
    opt.elcv_actual(feasible_index) = alpha.*(exp(g2(feasible_index))-1); 
    
    opt.actual_func = cell(1,7);
    opt.actual_func{1} = opt.cons_actual;
    opt.actual_func{2} = opt.cv_actual;
    opt.actual_func{3} = opt.acv_actual;
    opt.actual_func{4} = opt.elcv_actual;
    opt.actual_func{5} = opt.thcv_actual;
    opt.actual_func{6} = opt.aelcv_actual;
    opt.actual_func{7} = opt.athcv_actual;
    
    %We only care about CV of those functions,
    %Therefore, calculate CV for those functions
    %%{
    opt.cv_func = cell(1,8);
    opt.cv_func{1} = opt.cv_actual;
    g = opt.cons_pred;
    g(g<0) = 0;
    opt.cv_func{2} = sum(g, 2);
    
    g = opt.cv_pred;
    g(g<0) = 0;
    opt.cv_func{3} = sum(g, 2);
    
    
    g = opt.acv_pred;
    g(g<0) = 0;
    opt.cv_func{4} = sum(g, 2);
    
    g = opt.elcv_pred;
    g(g<0) = 0;
    opt.cv_func{5} = sum(g, 2);
    
    
    g = opt.thcv_pred;
    g(g<0) = 0;
    opt.cv_func{6} = sum(g, 2);
    
    g = opt.aelcv_pred;
    g(g<0) = 0;
    opt.cv_func{7} = sum(g, 2);
    
    
    g = opt.athcv_pred;
    g(g<0) = 0;
    opt.cv_func{8} = sum(g, 2);
    %}

end


function opt = build_models(opt, TrainIndex)

    %training set
    x = opt.archive(TrainIndex,:);
    %f = opt.archiveObj(TrainIndex, :);
    cons = opt.archiveCons(TrainIndex,:);
    
    %find unique solutions
    [~,ia,~] = unique(x,'rows');
    x = x(ia,:);
    cons = cons(ia,:);
    %f = f(ia,:);    
    
    %find CV 
    g = cons;
    g(g<0) = 0;
    cv = sum(g, 2);
    acv = cv;

    %find ACV
    feasible_index = find(cv==0);
    g2 = sum(cons, 2);
    acv(feasible_index) = g2(feasible_index);
    
    %tanh hyperbolic tangent constraint violation function
    thcv = tanh(cons);
    
    %aggregated hyperbolic tangent constraint violation function
    athcv = tanh(acv);
    
    %exponential linear cv
    elcv = cons;
    alpha = 0.5;
    elcv(feasible_index) = alpha.*(exp(g2(feasible_index))-1); 
    
    %aggregated exponential linear cv
    aelcv = cv;
    alpha = 0.5;
    aelcv(feasible_index) = alpha*(exp(g2(feasible_index))-1);     
   
    %----------------MODEL CONSTRAINTS-------------------------------------
    [opt.dmodel_cons, ~] = dacefit(x, cons, @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model cons
    
    
    %----------------MODEL CONSTRAINT VIOLATION FUNCTION-------------------
    [opt.dmodel_cv, ~] = dacefit(x, cv, @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model cv
    
    
    %----------------MODEL AGGREGATED CONSTRAINT FUNCTION------------------    
    [opt.dmodel_acv, ~] = dacefit(x, acv, @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model acv
    
    
    %----------------MODEL EXPONENTIAL LINEAR CONSTRAINT FUNCTION----------
    [opt.dmodel_elcv, ~] = dacefit(x, elcv, @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model elcv
    
    %-----MODEL AGGREGATED EXPONENTIAL LINEAR CONSTRAINT FUNCTION----------
    [opt.dmodel_aelcv, ~] = dacefit(x, aelcv, @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model elcv
    
    
    %--------------MODEL HYPERBOLIC CONSTRAINT VIOLATION FUNCTION----------
    [opt.dmodel_thcv, ~] = dacefit(x, thcv, @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model thcv
    
    
    %--------------MODEL HYPERBOLIC CONSTRAINT VIOLATION FUNCTION----------
    [opt.dmodel_athcv, ~] = dacefit(x, athcv, @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model thcv
    
    
    %---------------MODEL ASF FOR EACH DIRECTIONS--------------------------
    
    %[opt.dmodel_asf, ~] = dacefit(x, y_asf, @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model asf       
    %[opt.dmodel_asfcv, ~] = dacefit(x, y_asfcv, @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model asfcv
    
    
    %---------------MODEL SELECTION FUNCTION-------------------------------
    %opt.net = train_neural_network(x, SF_VAL);
    
end



%==========================================================================
%===================GENERATE TEST SET======================================
%==========================================================================

function [test_x, test_obj, test_cons] = generate_meshgrid_data(opt, N)
    %create test samples as 2D mesh grid
    x1 = linspace(opt.bound(1,1),opt.bound(2,1),N);
    x2 = linspace(opt.bound(1,2),opt.bound(2,2),N);
    [X, Y] = meshgrid(x1, x2);
    test_x = [X(:) Y(:)];
    test_x = test_x+1e-16;
    test_obj = zeros(size(test_x,1), opt.M);
    test_cons = zeros(size(test_x,1), max(1,opt.C));

    for i=1:size(test_x,1)
        if opt.C>0
            [test_obj(i,:), test_cons(i,:)] = high_fidelity_evaluation(opt, test_x(i,:));
        else
            [test_obj(i,:), ~] = high_fidelity_evaluation(opt, test_x(i,:)); 
        end
    end
end


%==========================================================================
%===================GENERATE LINEAR DATA===================================
%==========================================================================

function [test_x, test_obj, test_cons] = generate_linearized_data(opt, N, dim)
    %create test samples as 2D mesh grid
    test_x = linspace(opt.bound(1,dim),opt.bound(2,dim),N);
    test_x = test_x+1e-16;
    test_obj = zeros(size(test_x,1), opt.M);
    test_cons = zeros(size(test_x,1), max(1,opt.C));

    for i=1:size(test_x,1)
        if opt.C>0
            [test_obj(i,:), test_cons(i,:)] = high_fidelity_evaluation(opt, test_x(i,:));
        else
            [test_obj(i,:), ~] = high_fidelity_evaluation(opt, test_x(i,:)); 
        end
    end
end

%==========================================================================
%============================LOAD PARAMETERS===============================
%==========================================================================

function [opt] = load_parameters(opt)

    %-----METAMODEL PARAMETERS---------------------------------------------
    
    opt.theta = 10;% starting value of Kriging correlation parameter
    opt.lob = 1e-2;%lower bound of Kriging correlation parameter
    opt.upb = 20;%20;%upper bound of Kriging correlation parameter
    opt.consOption = 1; %1  = CV, 2 = mCV, option for modeling constraints
    opt.rho = 0;%0.001;%----rho value for augment asf MEMO approach
    %opt.metamodel = 1%1=Kriging Model
    opt.funcEval = 0;
    opt.iter = 1;
    
    
    %------OPTIMIZATION ALGORITHM PARAMETERS-------------------------------
     
    opt.eta_c = 15;%crossover index
    opt.eta_m = 20;%mutation index
    opt.G = 100;%100;%number of generations 
    opt.N = 100;%100;%100;%50;%100;%population size in optimization algorithm
    opt.pcross = 0.95; % Crossover probability
    opt.nrealcross = 0;%number of crossover performed
    opt.nrealmut = 0;%number of mutation performed
    opt.gen = 1;
    opt.pop = [];
    opt.popObj = [];
    opt.changeDelta = 200;%change delta in every ## function evaluation
    opt.epsilonTrust = 100;%controls smoothness of trusted region, larger means isolated
    opt.SampleSize = 10;%how many points are minimum needed around an expensively evaluated point to build the model
    opt.NumExplorationPoints = 3;
    opt.hybrid = 8;%it can be changed to value: 3
    opt.switch_method = 1000;%switch method after this much iteration
    opt.phase = 1;
    opt.Epsilon = 1e-14;
    
    
    %-----OBJECTIVE FUNCTION PARAMETERS------------------------------------
    
    switch(opt.objfunction)
        
        case {'c2dtlz2','c3dtlz2'}
            opt.M = 2;
            opt.V = 2; 
            
            opt.totalFuncEval = 1500;
            opt.utopian = (-0.05).*ones(1, opt.M);
            opt.min_val = zeros(1, opt.M);
            opt.initpopsize = 700;
            
            if strcmp(opt.objfunction,'c2dtlz2')
                opt.C = 1;
                opt.max_val = ones(1, opt.M);
            else
                opt.C = 3;
                opt.max_val = (2.01)*ones(1, opt.M);
            end
        
        case 'bnh'  
            opt.M = 2;
            opt.V = 2; 
            opt.C = 2;
            opt.totalFuncEval = 800;
            opt.utopian = [-0.05, -0.05];
            opt.min_val = [0 0];
            opt.max_val = [140 55];
            opt.initpopsize = 200;
            
        case 'osy'
            opt.M = 2;
            opt.V = 6;
            opt.C = 6;
            opt.utopian = [-300 0];
            opt.totalFuncEval =  800;%1250;
            opt.min_val = [-273.8 4];
            opt.max_val = [-42 76];
            opt.initpopsize = 200;
        case 'srn'
            opt.M = 2;
            opt.V = 2; 
            opt.C = 2;
            opt.totalFuncEval = 800;
            opt.utopian = [0 -300];
            opt.min_val = [0 -250];
            opt.max_val = [240 0];
            opt.initpopsize = 200;
        case 'tnk'
            opt.M = 2;
            opt.V = 2; 
            opt.C = 2;
            opt.totalFuncEval = 800;
            opt.utopian = [-0.001, -0.001];
            opt.min_val = [0 0];
            opt.max_val = [1.2 1.2];
            opt.initpopsize = 200;
        case 'water'
            opt.M = 5;
            opt.V = 3; 
            opt.C = 7;
            opt.totalFuncEval = 300;
            opt.utopian = (-0.05).*ones(1, opt.M);
            opt.min_val = [0.75 0 0 0 0];
            opt.max_val = [0.95 0.9 1.0 1.6 3.2];
            opt.initpopsize = 100;
        case 'carside'
            opt.M = 3;
            opt.V = 7; 
            opt.C = 10;
            opt.totalFuncEval = 2000;
            opt.utopian = [24.3180    3.5352   10.5610];
            opt.min_val = [24.3680    3.5852   10.6110];
            opt.max_val = [42.7620    4.0000   12.5210];
            opt.initpopsize = 700;
            
        case 'welded'
            opt.M = 2;
            opt.V = 4; 
            opt.C = 4;
            opt.totalFuncEval = 1000;
            opt.utopian = [2.3316   -0.0496];
            opt.min_val = [2.3816    0.0004];
            opt.max_val = [36.4403    0.0157];
            opt.initpopsize = 200;

        otherwise
            input('Function definition is not found');
    end
    
    %--------LOWER AND UPPER BOUND OF DECISION VARIABLE--------------------
    opt.bound = zeros(2,opt.V);
    opt.bound(2,:) = ones(1,opt.V);
    if(strcmpi(opt.objfunction,'zdt4'))
        opt.bound(1,2:end) = opt.bound(1,2:end)+(-5);
        opt.bound(2,2:end) = opt.bound(2,2:end)*5;
    elseif(strcmpi(opt.objfunction,'UF1') ||strcmpi(opt.objfunction,'UF2') || strcmpi(opt.objfunction,'UF5') || strcmpi(opt.objfunction,'UF6') || strcmpi(opt.objfunction,'UF7'))
        opt.bound(1,2:end) = opt.bound(1,2:end)+(-1);
        opt.bound(2,2:end) = opt.bound(2,2:end)*1;
    elseif(strcmpi(opt.objfunction,'UF4') )
        opt.bound(1,2:end) = opt.bound(1,2:end)+(-2);
        opt.bound(2,2:end) = opt.bound(2,2:end)*2;
    elseif(strcmpi(opt.objfunction,'UF8') || strcmpi(opt.objfunction,'UF9')|| strcmpi(opt.objfunction,'UF10'))
        opt.bound(1,3:end) = opt.bound(1,3:end)+(-2);
        opt.bound(2,3:end) = opt.bound(2,3:end)*2;
    elseif(strcmpi(opt.objfunction,'BNH'))
        opt.bound(2,1)=5;
        opt.bound(2,2)=3;
    elseif(strcmpi(opt.objfunction,'OSY'))
        opt.bound(2,1)=10;%x1
        opt.bound(2,2)=10;
        opt.bound(2,6)=10;
        opt.bound(1,3)=1;
        opt.bound(1,5)=1;
        opt.bound(2,3)=5;
        opt.bound(2,5)=5;
        opt.bound(2,4)=6;
    elseif(strcmpi(opt.objfunction,'SRN'))
        opt.bound(1,1:end) = opt.bound(1,1:end)+(-20);
        opt.bound(2,1:end) = opt.bound(2,1:end)*20;
    elseif(strcmpi(opt.objfunction,'TNK'))
        opt.bound(1,1:end) = opt.bound(1,1:end)+1e-12;
        opt.bound(2,1:end) = opt.bound(2,1:end)*pi;
    elseif(strcmpi(opt.objfunction,'WATER'))
        opt.bound(1,1:end) = opt.bound(1,1:end)+0.01;
        opt.bound(2,1) = 0.45;   
        opt.bound(2,2:end) = 0.10; 
    elseif(strcmpi(opt.objfunction,'carside'))
        %opt.bound(1,[1 3 4]) = 0.5; opt.bound(1,2) = 0.45; opt.bound(1,[6 7]) = 0.4; opt.bound(1,5) = 0.875;
        %opt.bound(2,[1 3 4]) = 1.5; opt.bound(2,2) = 1.35; opt.bound(2,[6 7]) = 1.2; opt.bound(2,5) = 2.625;  
        opt.bound(1,1:end) = [0.5 0.45 0.5 0.5 0.875 0.4 0.4];
        opt.bound(2,1:end) = [1.5 1.35 1.5 1.5 2.625 1.2 1.2];
    elseif(strcmpi(opt.objfunction,'welded'))
        opt.bound(1,1:end) = [0.125 0.1 0.1 0.125];
        opt.bound(2,1:end) = [5 10 10 5];
    elseif    (strcmpi(opt.objfunction,'wfg1') || strcmpi(opt.objfunction,'wfg2') || strcmpi(opt.objfunction,'wfg3') ...
            || strcmpi(opt.objfunction,'wfg4') || strcmpi(opt.objfunction,'wfg5') || strcmpi(opt.objfunction,'wfg6')...
            || strcmpi(opt.objfunction,'wfg7') || strcmpi(opt.objfunction,'wfg8') || strcmpi(opt.objfunction,'wfg9'))
        opt.bound(1,:) = zeros(1,opt.V);
        opt.bound(2,:) = ones(1,opt.V);
    elseif strcmpi(opt.objfunction,'toy')
        opt.bound(1) = 0;
        opt.bound(2) = 1;
    end
    
    %---------PARETO FRONT-------------------------------------------------
    if opt.M == 2
        opt.paretofront = load(strcat(upper(opt.objfunction),'.2D.pf'));%Pareto front for IGD calculation
        %opt.paretofront = calculate_target_pareto_points(strcat('IGD Calculation/GD/',upper(opt.objfunction),'.2D.new.pf'), opt.dirs);
        opt.paretofrontGD = load(strcat('GD/',upper(opt.objfunction),'.2D.pf'));%Pareto front for GD calculation
    elseif opt.M == 3
        opt.paretofront = load(strcat(upper(opt.objfunction),'.3D.pf'));
        opt.paretofrontGD = load(strcat('GD/',upper(opt.objfunction),'.3D.pf'));
    elseif opt.M == 5
        opt.paretofront = load(strcat(upper(opt.objfunction),'.5D.pf'));
        opt.paretofrontGD = load(strcat('GD/',upper(opt.objfunction),'.5D.pf'));
    end
    
    %---------DIRECTION----------------------------------------------------
    
    if opt.M==5
        opt.dirs = initweight(5, 210)';
        opt.refdiv = 6;
    elseif opt.M==3 %Three obj
        opt.dirs = initweight(3, 91)';
        opt.refdiv = 12;
    else
        %opt.dirs = [0.5 0.5];
        opt.dirs = initweight(2, 21)';%opt.dirs(11,:)=[];opt.dirs(1,:)=[];
        %opt.dirs = [0.75 0.25;0.5 0.5;0.25 0.75];%initweight(2, 11)';
        opt.refdiv = 20;
    end
    
    opt.numdir = size(opt.dirs,1);%number of reference direction
    opt.curdir = opt.dirs(1,:);%current direction
    opt.curcluster = 1;%current cluster number
    
    opt.pmut = 1.0/opt.V; % Mutation probability
    
    opt.ref_point = cell(1, opt.numdir);
    %-------------------ALLOCATE-------------------------------------------
    
    opt.archive = [];
    opt.archiveObj = [];
    opt.archiveCV = [];
    opt.archiveASF = [];
    opt.archiveCVModified = [];
    opt.archiveCluster = [];
    opt.archiveCons = [];
    opt.archiveKKTPM = [];
    opt.activearchive=[];
    opt.activearchiveCluster = [];
    opt.normalizedObj = [];
    opt.ybest = [];
    opt.Pareto = [];
    opt.ParetoVar = [];
    opt.Color = {'k','b','r','g','y',[.5 .6 .7],[.8 .2 .6]}; %Colors.
    opt.InitialActiveSetSize = 300; %opt.initpopsize;
    opt.activeSetSize = opt.InitialActiveSetSize;
    if opt.C>0
        opt.regC = cell(1,opt.C);
    end
    opt.regO = cell(1,opt.M);                
end


function [opt] = store_results_for_plot(opt, pop, pop_obj, pop_cons)


    %----------------EVALUATE POPULATION-----------------------------------
    %{
    if ~isempty(opt.archive)
        [~,~,index] = intersect(opt.archive,pop,'rows');
        if ~isempty(index)
            if ~size(index,1)==size(pop,1)
                index = ~index;
                pop = pop(index,:);
                pop_obj = pop_obj(index,:);
                pop_cons = pop_cons(index,:);
            end
        end
    end
    %}
    opt.archive = vertcat(opt.archive, pop); %archive population variables  
    opt.archiveObj = vertcat(opt.archiveObj, pop_obj); %archive objective values
    
    if opt.C>0
        [pop_cv, pop_cv_modified] = evaluateCV(pop_cons);%find constraint violation (CV) for new solutions
        opt.archiveCV = vertcat(opt.archiveCV, pop_cv); %archive constrain violation
        opt.archiveCVModified = vertcat(opt.archiveCVModified, pop_cv_modified);
        opt.archiveCons = vertcat(opt.archiveCons, pop_cons);%archive constraint values  
        opt.ybest_cv = min(opt.archiveCV);%minimum of all
        opt.archiveCVMatrix = evaluateCVhere(opt, opt.archiveCons);
    else
        pop_cv = zeros(size(pop,1),1);
        opt.archiveCV = vertcat(opt.archiveCV, pop_cv); %archive constrain violation
        opt.archiveCVModified =  vertcat(opt.archiveCVModified, pop_cv); 
        opt.archiveCons = vertcat(opt.archiveCons, pop_cons);%archive constraint values  
    end
    
    
    %---------------------WRITE TO FILE------------------------------------
%     opt.writeFlagOption = 1;
%     if opt.writeFlagOption==1
%         dlmwrite(opt.varfilename, pop, 'delimiter',' ','precision','%.10f','-append');
%         dlmwrite(opt.objfilename, pop_obj, 'delimiter',' ','precision','%.10f','-append');
%         dlmwrite(opt.cvfilename, pop_cv, 'delimiter', ' ','precision','%.10f','-append');
%     end
%     
    %-------------FIND FEASIBLE PARETO FRONT-------------------------------
    
    if opt.C>0 %Feasible Pareto front of Constraint Problems
        index = find(opt.archiveCV<=0);
    else
        index = (1:size(opt.archive,1))';
    end
    
    opt.FeasibleIndex = index;
    
    if size(index,1)>0 % there are some feasible solutions
        index2 = paretoFront(opt.archiveObj(index,:));
        opt.ParetoIndex = index(index2)';
        %opt.min_val = 0.9*min(opt.archiveObj(opt.ParetoIndex,:),[],1);% - 1;%1e-6;
        %opt.max_val = 1.1*max(opt.archiveObj(opt.ParetoIndex,:),[],1);
        
        opt.nadir_point = opt.archiveObj(opt.ParetoIndex,:);
        if size(opt.nadir_point,1)>1
            opt.nadir_point =  max(opt.nadir_point);
        end
            

    else % no feasible solution yet
        opt.ParetoIndex = [];
        %opt.min_val = 0.9*min(opt.archiveObj,[],1);% - 1;%1e-6;
        %opt.max_val = 1.1*max(opt.archiveObj,[],1);
        opt.nadir_point = max(opt.archiveObj);
    end

    
    
    
    for i=1:size(opt.dirs,1)
        w = opt.dirs(i,:);
        w_org = zeros(1,opt.M);
        for j=1:opt.M
            w_org(j) =  opt.min_val(j)+ w(j)*(opt.max_val(j)-opt.min_val(j));
        end
        %Z  = max((opt.archiveObj-repmat(opt.min_val,size(opt.archiveObj,1),1))./repmat(w_org,size(opt.archiveObj,1),1),[],2);
        %Z  = max(opt.normalizedObj./repmat(w_org,size(opt.normalizedObj,1),1),[],2);
        small_z = (opt.max_val-opt.min_val)/2;
        opt.ref_point{i} = w_org - small_z;
    end
    
    %opt = unique_population(opt);
    
    %-----------------Find Feasible Solution-------------------------------
    if size(opt.FeasibleIndex,1)>0
        if size(opt.FeasibleIndex,1)>1
            opt.ideal_point = min(opt.archiveObj(opt.FeasibleIndex,:));
        else
            opt.ideal_point = opt.archiveObj(opt.FeasibleIndex,:);
        end
    else
        opt.ideal_point = min(opt.archiveObj);
    end
    
    %---------NORMALIZE AND FIND ASF IN NORMALIZED SPACE-------------------
    min_val = min(opt.archiveObj);
    max_val = max(opt.archiveObj);
    opt.normalizedObj = opt.archiveObj;%normalize(opt, opt.archiveObj,  min_val, max_val);%opt.archiveObj;%normalize(opt, opt.archiveObj,  opt.min_val, opt.max_val);
    opt.normalizedPop = opt.archive;%normalize(opt, opt.archive, opt.bound(1,:), opt.bound(2,:));
    %opt.normalizedObj_nsga3 =  nsga3_normalization(opt, opt.archiveObj);
    
    
    [opt] = evaluateASFAll(opt);
    [opt.archiveASF, opt.archiveCluster] = evaluateASF(opt);%evaluate ASF for all solutions, previous and new
    
    opt.selection_function_option = 2;
    [opt.archiveSF_VAL] = selectionFunction(opt, pop);
    
    %{
    if opt.methodology == 6 || opt.methodology == 7 || opt.methodology == 8
        
        %-------------FIND BEST SOLUTION FOR EACH CLUSTER------------------
        index = cell(1,opt.numdir);
        opt.bestPop = cell(1,opt.numdir);%best population for each cluster
        opt.bestpopObj = cell(1,opt.numdir);%best population for each cluster
        opt.bestpopASF = zeros(1, opt.numdir);%best ASF for each cluster
        
        for i = 1:opt.numdir %number of clusters
            index{i} = find(opt.archiveCluster == i); %all solutions in cluster i
            
            if ~isempty(index{i})
                %obj = opt.archiveKKTPM(index{i},:); %objectives of cluster i
                obj = opt.archiveASF(index{i},:);

                [~,I] = sort(obj(:,1)); %should be non-dominated sort for multiple objectives
                opt.bestPop{i} = opt.archive(index{i}(I(1),:),:);
                opt.bestpopObj{i} = opt.archiveObj(index{i}(I(1),:),:);
                opt.bestpopASF(i) = opt.archiveASF(index{i}(I(1)));
            else
                opt.bestPop{i} = [];
                opt.bestpopObj{i} = [];
            end
        end

    else
        opt.bestpopObj = min(opt.archiveASF);
    end
    %}
             
    
    opt.activeArchive = opt.archive;
    opt.activeArchiveObj = opt.archiveObj;
    opt.activeArchiveASF = opt.archiveASF;
    opt.activeArchiveCons = opt.archiveCons;
    opt.activeArchiveCV = opt.archiveCV;
    %------------STORE UNIQUE RESULTS--------------------------------------
    %opt = unique_population(opt);
    %--------------Find Leader---------------------------------------------
    %opt = find_leader(opt);
    
    opt.activeSetSize = min(opt.InitialActiveSetSize, size(opt.activeArchive, 1));
    
end

function [opt] = evaluateASFAll(opt)
    
    opt.archiveASFAll = cell(1, opt.numdir);
    opt.archiveASFCVAll = cell(1, opt.numdir);
    for i=1:size(opt.dirs,1)
        w = opt.dirs(i,:);
        %asf = calculate_Asf(opt.archiveObj, w);
        asf = calculate_Asf(opt.normalizedObj, w);
        opt.archiveASFAll{i} = asf; 
    end
    
    opt.archiveASFCVAll = opt.archiveASFAll;
    feasbile_index = opt.archiveCV<=0;
    
    for i=1:size(opt.dirs,1)
        if opt.C>0
            feasibleASF =  opt.archiveASFAll{i}(feasbile_index,:);
            
            if ~isempty(feasibleASF)
                fmax = max(feasibleASF);%maximum feasible ASF
            else
                fmax = 0;
            end
            infeasible_index = find(opt.archiveCV>0);%infeasible indices
            if ~isempty(infeasible_index)
                opt.archiveASFCVAll{i}(infeasible_index) = fmax + opt.archiveCV(infeasible_index);
            end
        end
    end
    
end

function asf = calculate_Asf(obj, dir)   
    dir = dir+1e-16;
    %asf = max( (obj-repmat(dir, size(obj,1),1)),[],2);
    asf = max( (obj./dir),[],2);
end

function cv = evaluateCVhere(opt, pop_cons)
    
    g = pop_cons;
    if opt.consOption==1 %negative values treated as zero
        for i=1:size(pop_cons,2)
           g(g(:,i)<0,i)=0;
        end
    end  

    cv = g;
end


%==========================================================================
%===============================END========================================
%==========================================================================
