function plot_framework()

    scalable_function = {'DTLZ1','DTLZ2','DTLZ3', 'DTLZ4', 'SDTLZ1', 'SDTLZ2','CDTLZ2'...
                  'IDTLZ1','IDTLZ2','C2_DTLZ2', 'C3_DTLZ4', 'C1_DTLZ1', ...
                   'DTLZ5',  'DTLZ7', 'DTLZ5IM', 'DTLZ8',  'DTLZ6', 'DTLZ9',...
                   'ZDT1', 'ZDT2','ZDT3','ZDT4','ZDT6'};   
    test_function = ['osy    '; 'zdt1   ';'tnk    ';'zdt2   ';'zdt3   ';'zdt6   ';   'bnh    ';];%test functions to plot
    framework_name = ['M1-1';'M1-2';'M2-1';'M2-2';'M3-1';'M3-2';'M4-1';'M4-2';'M5  ';'M6  '];
    method_list = [11, 12, 21, 22, 31, 32, 41, 42];
    addpath(genpath('Problems'));
    addpath('Metrics');
    opt.r = 1;%run number
    param.M = 2;
    param.V = 2;
    param.C = 0;
    opt.initpopsize = 100;
    opt.totalFuncEval = 500;
    random_initial = 1;
     
    for methodology = 1:1:1
        
        opt.methodology = method_list(methodology);
        
        for func_no = 2:2
            
            %------------------Problem Parameters--------------------------
            
            opt.func_no = func_no;%function number
            opt.objfunction = strtrim(test_function(func_no,:));%remove whitespaces
            opt.Global = [];
            opt.Global = feval(str2func(upper(opt.objfunction)), 'init', opt.Global, []);

            if ismember(opt.objfunction, scalable_function)
                opt.M = param.M;
                opt.C = param.C;
                opt.V = param.V;
            else
                opt.V = 2;%opt.Global.D;
                opt.M = opt.Global.M;
                opt.C = opt.Global.C;
            end
            opt.bound = zeros(2, opt.V);           
            opt.bound(1,:) = opt.Global.lower(1:opt.V);
            opt.bound(2,:) = opt.Global.upper(1:opt.V);
            opt.writeFlag = 0;
            opt.plotOption = 0;
            opt.initOption = 1;
            opt = app_parameters(opt);%load other parameters
            
            
            %-----------------Create Dataset-------------------------------
            if random_initial==0
                opt.initpopsize = 100;%LHS sampling
                [pop, popObj, popCons] = initialization(opt);%LHS sampling  
                opt = store_results(opt, pop, popObj, popCons);%save results,asf,cv
                TrainIndex = 1:100;
                opt.TrainIndex = TrainIndex;
            else
                N = 20;
                [test_x, test_obj, test_cons] = generate_meshgrid_data(opt, N);
                opt = store_results(opt, test_x, test_obj, test_cons);%save results,asf,cv
                index = randperm(size(opt.archive,1));
                opt.TrainIndex = index(1:30);
                TrainIndex = opt.TrainIndex;
                TestIndex = 1:size(opt.archive,1);%opt.initpopsize+1:size(opt.archive,1);
                opt.TestIndex = TestIndex;
            end
            
            %--------------Create Metamodels-------------------------------
            
            opt = collect_models(opt, TrainIndex);
        
            %--------------compute selection values------------------------
            
            [opt, rank_test, rank_test_pred] = calculate_rank(opt, TestIndex);
            [opt, rank_train, rank_train_pred] = calculate_rank(opt, TrainIndex);
        
            %--------------------Plot Selection Function-------------------
            
            for i=6:2:10
                if mod(i, 2)==1
                    %plot for 11-th direction
                    for j=1:opt.numdir
                        plot_sf_framework(opt, rank_test{i}{j}, rank_test_pred{i}{j}, ...
                        rank_train{i}{j}, rank_train_pred{i}{j}, framework_name(i,:));
                    end
                else
                    plot_sf_framework(opt, rank_test{i}, rank_test_pred{i}, ...
                    rank_train{i}, rank_train_pred{i}, framework_name(i,:));
                end
            end
        end
        
    end

end
