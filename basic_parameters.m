function [opt] = basic_parameters(opt)

    
    %======================================================================
    %=================METAMODELING FRAMEWORK OPTIONS=======================
    %======================================================================
    %opt.methodology = 11;%<===================select method here
    opt.framework_list  = opt.methodology;
    %opt.framework_list = [11,12,21,22,31,32,41,42,5,6];
    %opt.framework_list = [11,12, 21,31,41,5];
    %opt.framework_list = [11];%[12,22,32,42,6];
    %delete(gcp('nocreate'));
    opt.num_of_frameworks = size(opt.framework_list,2);
    opt.total_num_of_frameworks = 10;
    
    if opt.num_of_frameworks==1
        opt.algo_name = strcat('framework',num2str(opt.methodology));
    else
        opt.algo_name = strcat('switching_framework_', regexprep(num2str(opt.framework_list), ' +', '_'));
    end
    
    
    %======================================================================
    %===================GENERAL PARAMETERS=================================
    %======================================================================
    
    opt.funcEval = 0;
    opt.gen = 1;
    opt.Epsilon = 1e-16;
    opt.Inf = 1e16;
    
    %======================================================================
    %=================LOW-FIDELITY OPTIMIZATION OPTIONS====================
    %======================================================================
    
    opt.initOption = 1;%1 = latin hypercube on decision space
    opt.metamodelOption = 1;% 1 = solutions will be evaluated by metamodel, 2 = high-fidelity
    opt.crossoverOption = 1;% 1 = simulated binary crossover
    opt.mutationOption = 1;% 1 = polynomial mutation
    opt.matingselectionOption = 1;%1 = binary constraint tournament selection
    opt.survivalselectionOption = 1;%1 = NSGA-II, 2 = NSGA-III
    opt.selection_function_option = 2;%1=ASFCV, 2=MEMO, 3=NSGA-II
    
    %======================================================================
    %=============LOW FIDELITY RGA ALGORITHM PARAMETERS====================
    %======================================================================
    
    opt.eta_c = 15;%crossover index
    opt.eta_m = 20;%mutation index
    opt.G = 100;%100;%number of generations 
    opt.N = 100;%100;%100;%50;%100;%population size in optimization algorithm
    opt.pcross = 0.95; % Crossover probability
    opt.nrealcross = 0;%number of crossover performed
    opt.nrealmut = 0;%number of mutation performed
    
    %======================================================================
    %=================LOW-FIDELITY NSGA-III PARAMETERS=====================
    %======================================================================
    
    
    %======================================================================
    %==============LOW-FIDELITY TRUST REGION PARAMETERS====================
    %======================================================================
    
    
    
    %======================================================================
    %==================TEST PROBLEMS PARAMETERS============================
    %======================================================================
    
    switch(opt.objfunction)
        
        case {'zdt1','zdt2','zdt3','zdt4','zdt6'}
            
            opt.M = 2;%number of objectives
            opt.V = 10;%;10;%number of variables
            opt.C = 0;%number of constraints
            opt.totalFuncEval = 500;%121;%high-fidelity function evaluation
            opt.utopian = [-0.05, -0.05];%un achivable point
            opt.min_val = [0 0];%minimum value for normalization
            opt.max_val = [1 1];%maximum objective 
            opt.initpopsize = opt.V*10;%initial sample size for high fidelity computation
            opt.bound = zeros(2,opt.V);
            opt.bound(2,:) = ones(1,opt.V);
            switch(opt.objfunction)
                case 'zdt3'
                    opt.min_val = [0 -1];
                    opt.max_val = [1 1];
                    opt.utopian = [-0.05, -1.1];
                case 'zdt4'
                    opt.V = 5;
                    opt.bound = zeros(2,opt.V);
                    opt.bound(2,:) = ones(1,opt.V);
                    opt.bound(1,2:end) = opt.bound(1,2:end)+(-5);
                    opt.bound(2,2:end) = opt.bound(2,2:end)*5;
                case 'zdt6'
                    opt.min_val = [0.25 0];
                    opt.max_val = [1 1];
            end
            
        case {'dtlz1', 'dtlz2', 'dtlz3','dtlz4', 'sdtlz1', 'sdtlz2','cdtlz2','c1_dtlz1','c1_dtlz3', 'c2_dtlz2'}
            
            opt.C = 0;
            opt.Dim     = [2, 3, 5, 8, 10, 15];
            opt.V_all{1} = [10, 10, 10, 10, 10, 10];%dtlz1
            opt.V_all{2} = [10, 10, 10, 10, 10, 10];%dtlz2
            opt.V_all{3} = [10, 10, 10, 10, 10, 10];%dtlz3
            opt.V_all{4} = [10, 10, 10, 10, 10, 10];%dtlz4
            opt.V_all{5} = [10, 10, 10, 10, 10, 10];%sdtlz1
            opt.V_all{6} = [10, 10, 10, 10, 10, 10];%sdtlz2
            opt.V_all{7} = [10, 10, 10, 10, 10, 10];%cdtlz2
            opt.V_all{8} = opt.Dim+4;%c1_dtlz1
            opt.V_all{9} = opt.Dim+9;%c1_dtlz3
            opt.V_all{10} = 7;%opt.Dim+9;%c1_dtlz3

            opt.M = opt.Dim(opt.dim);
            opt.PopSize = [100, 100, 100, 100, 100, 100, 100, 100, 100, 100];%initial popsize
            opt.initpopsize = opt.PopSize(opt.dim);
            opt.funcEval_all{1} = [200, 300, 500, 700, 800, 1000];%dtlz1
            opt.funcEval_all{2} = [200, 300, 500, 700, 800, 1000];%dtlz2
            opt.funcEval_all{3} = [200, 300, 500, 700, 800, 1000];%dtlz3
            opt.funcEval_all{4} = [200, 300, 500, 700, 800, 1000];%dtlz4
            opt.funcEval_all{5} = [200, 300, 500, 700, 800, 1000];%sdtlz1
            opt.funcEval_all{6} = [200, 300, 500, 700, 800, 1000];%sdtlz2
            opt.funcEval_all{7} = [200, 300, 500, 700, 800, 1000];%cdtlz2
            opt.funcEval_all{8} = [200, 500, 600, 800, 1000, 1500];%c1_dtlz1
            opt.funcEval_all{9} = [200, 1000, 1500, 2500, 3500, 5000];%c1_dtlz3
            opt.funcEval_all{10} = [1500, 1000, 1500, 2500, 3500, 5000];%c2_dtlz2
            
            switch(opt.objfunction)
                case 'dtlz1'
                    opt.V = opt.V_all{1}(opt.dim);
                    opt.totalFuncEval = opt.funcEval_all{1};
                case 'dtlz2'
                    opt.V = opt.V_all{2}(opt.dim);
                    opt.totalFuncEval = opt.funcEval_all{2};
                case 'dtlz3'
                    opt.V = opt.V_all{3}(opt.dim);
                    opt.totalFuncEval = opt.funcEval_all{3};
                case 'dtlz4'
                    opt.V = opt.V_all{4}(opt.dim);
                    opt.totalFuncEval = opt.funcEval_all{4};
                case 'sdtlz1'
                    opt.V = opt.V_all{5}(opt.dim);
                    opt.totalFuncEval = opt.funcEval_all{5};
                case 'sdtlz2'
                    opt.V = opt.V_all{6}(opt.dim);
                    opt.totalFuncEval = opt.funcEval_all{6};
                case 'cdtlz2'
                    opt.C =1;
                    opt.V = opt.V_all{7}(opt.dim);
                    opt.totalFuncEval = opt.funcEval_all{7};
                case 'c1_dtlz1'
                    opt.V = opt.V_all{8}(opt.dim);
                    opt.totalFuncEval = opt.funcEval_all{8};
                case 'c1_dtlz3'
                    opt.V = opt.V_all{9}(opt.dim);
                    opt.totalFuncEval = opt.funcEval_all{9};
                case 'c2_dtlz2'
                    opt.V = 7;
                    opt.initpopsize = 700;
                    opt.totalFuncEval = 2500;
                    opt.M = 5;
                otherwise
                    opt.V = 10 + opt.M-1;
                    opt.totalFuncEval = opt.funcEval_all{1};
            end
            
            opt.bound = zeros(2, opt.V);
            opt.bound(2,:) = ones(1,opt.V);    
        
        case {'dtlz5', 'dtlz7'}
            
            opt.M = 2;
            opt.V = 10;
            opt.C = 0;
            opt.totalFuncEval = 1000;
            opt.initpopsize = 200;
                
            if strcmp(opt.objfunction,'dtlz4')
                opt.initpopsize = 700;%initial sample size for high fidelity computation
                opt.totalFuncEval = 2000;
            elseif strcmp(opt.objfunction,'dtlz7')
                opt.M = 3;
            end
            opt.C = 0;
            opt.bound = zeros(2, opt.V);
            opt.bound(2,:) = ones(1,opt.V);  
            
        case {'wfg1', 'wfg2', 'wfg3','wfg4', 'wfg5', 'wfg6','wfg7','wfg8','wfg9'}
            
            opt.M = 3;
            opt.V = 8; 
            opt.C = 0;
            opt.totalFuncEval = 500;
            
            opt.initpopsize = 100;
            opt.bound(1,:) = zeros(1,opt.V);
            opt.bound(2,:) = ones(1,opt.V);
            
        case {'c3dtlz2'}
            
            opt.M = 3;
            opt.V = 7;
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
            
            opt.bound = zeros(2, opt.V);
            opt.bound(2,:) = ones(1,opt.V);   
        
        case 'bnh'  
            
            opt.M = 2;
            opt.V = 2; 
            opt.C = 2;
            opt.totalFuncEval = 500;
            opt.utopian = [-0.05, -0.05];
            opt.min_val = [0 0];
            opt.max_val = [140 55];
            opt.initpopsize = 100;
            opt.bound = zeros(2,opt.V);
            opt.bound(2,:) = ones(1,opt.V);
            opt.bound(2,1)=5;
            opt.bound(2,2)=3;
            
        case 'osy'
            
            opt.M = 2;
            opt.V = 6;
            opt.C = 6;
            opt.utopian = [-300 0];
            opt.totalFuncEval =  800;%1250;
            opt.min_val = [-273.8 4];
            opt.max_val = [-42 76];
            opt.initpopsize = 200;
            opt.bound = zeros(2,opt.V);
            opt.bound(2,:) = ones(1,opt.V);
            opt.bound(2,1)=10;%x1
            opt.bound(2,2)=10;
            opt.bound(2,6)=10;
            opt.bound(1,3)=1;
            opt.bound(1,5)=1;
            opt.bound(2,3)=5;
            opt.bound(2,5)=5;
            opt.bound(2,4)=6;
            
        case 'srn'
            
            opt.M = 2;
            opt.V = 2; 
            opt.C = 2;
            opt.totalFuncEval = 500;
            opt.utopian = [0 -300];
            opt.min_val = [0 -250];
            opt.max_val = [240 0];
            opt.initpopsize = 100;
            opt.bound = zeros(2,opt.V);
            opt.bound(2,:) = ones(1,opt.V);
            opt.bound(1,1:end) = opt.bound(1,1:end)+(-20);
            opt.bound(2,1:end) = opt.bound(2,1:end)*20;
            
        case 'tnk'
            
            opt.M = 2;
            opt.V = 2; 
            opt.C = 2;
            opt.totalFuncEval = 500;
            opt.utopian = [-0.001, -0.001];
            opt.min_val = [0 0];
            opt.max_val = [1.2 1.2];
            opt.initpopsize = 200;
            opt.bound = zeros(2,opt.V);
            opt.bound(2,:) = ones(1,opt.V);
            opt.bound(1,1:end) = opt.bound(1,1:end)+1e-12;
            opt.bound(2,1:end) = opt.bound(2,1:end)*pi;
            
        case 'water'
            
            opt.M = 5;
            opt.V = 3; 
            opt.C = 7;
            opt.totalFuncEval = 500;
            opt.utopian = (-0.05).*ones(1, opt.M);
            opt.min_val = [0.75 0 0 0 0];
            opt.max_val = [0.95 0.9 1.0 1.6 3.2];
            opt.initpopsize = 100;
            opt.bound = zeros(2,opt.V);
            opt.bound(2,:) = ones(1,opt.V);
            opt.bound(1,1:end) = opt.bound(1,1:end)+0.01;
            opt.bound(2,1) = 0.45;   
            opt.bound(2,2:end) = 0.10; 
            
        case 'carside'
            
            opt.M = 3;
            opt.V = 7; 
            opt.C = 10;
            opt.totalFuncEval = 500;
            opt.utopian = [24.3180    3.5352   10.5610];
            opt.min_val = [24.3680    3.5852   10.6110];
            opt.max_val = [42.7620    4.0000   12.5210];
            opt.initpopsize = 100;
            opt.bound = zeros(2,opt.V);
            opt.bound(2,:) = ones(1,opt.V);
            opt.bound(1,1:end) = [0.5 0.45 0.5 0.5 0.875 0.4 0.4];
            opt.bound(2,1:end) = [1.5 1.35 1.5 1.5 2.625 1.2 1.2];
            
        case 'welded'
            
            opt.M = 2;
            opt.V = 4; 
            opt.C = 4;
            opt.totalFuncEval = 500;
            opt.utopian = [2.3316   -0.0496];
            opt.min_val = [2.3816    0.0004];
            opt.max_val = [36.4403    0.0157];
            opt.initpopsize = 100;
            opt.bound = zeros(2,opt.V);
            opt.bound(2,:) = ones(1,opt.V);
            opt.bound(1,1:end) = [0.125 0.1 0.1 0.125];
            opt.bound(2,1:end) = [5 10 10 5];
            
        case {'g1','g2','g3','g4','g5','g6','g7','g8','g9','g10'}
            
            opt.G = 200;
            opt.M = 1;
            opt.N = 100;
            
            switch(opt.objfunction)
                case 'g1'
                    opt.V = 13;
                    opt.C = 9;
                    opt.bound = zeros(2, opt.V);
                    opt.bound(2,1:9) = ones(1,9);
                    opt.bound(2,10:12) = 100*ones(1,3);
                    opt.bound(2,13) = 1;
                    
                case 'g2'
                    opt.V = 5;%20;
                    opt.C = 2;
                    opt.bound = zeros(2, opt.V);
                    opt.bound(2,1:opt.V) = 10*ones(1,opt.V);
                    
                case 'g3'
                    opt.V = 10;
                    opt.C = 1;
                    opt.bound = zeros(2, opt.V);
                    opt.bound(2,1:opt.V) = ones(1,opt.V);

                case 'g4'
                    opt.V = 5;
                    opt.C = 6;
                    opt.bound = zeros(2, opt.V);
                    opt.bound(1,1)  = 78;
                    opt.bound(2,1)  = 102;
                    opt.bound(1,2)  = 33;
                    opt.bound(2,2)  = 45;
                    opt.bound(1,3:5)  = 27;
                    opt.bound(2,3:5)  = 45;

                case 'g5'
                    opt.V = 4;
                    opt.C = 8;
                    opt.bound = zeros(2, opt.V);
                    opt.bound(2,1)  = 1200;
                    opt.bound(2,2)  = 1200;
                    opt.bound(1,3:4)  = -0.55;
                    opt.bound(2,3:4)  =  0.55;
                    
                case 'g6'
                    opt.V = 2;
                    opt.C = 2;
                    opt.bound = zeros(2, opt.V);
                    opt.bound(1,1)  = 13;
                    opt.bound(2,1)  = 100;
                    opt.bound(1,2)  = 0;
                    opt.bound(2,2)  = 100;
                
                case 'g7'
                    opt.V = 10;
                    opt.C = 8;
                    opt.bound = zeros(2, opt.V);
                    opt.bound(1,1:opt.V) = -10;
                    opt.bound(2,1:opt.V) = 10;
                    
                case 'g8'
                    opt.V = 2;
                    opt.C = 2;
                    opt.bound = zeros(2, opt.V);
                    opt.bound(2,1:opt.V) = 10;
                    
                case 'g9'
                    opt.V = 7;
                    opt.C = 4;
                    opt.bound = zeros(2, opt.V);
                    opt.bound(1,1:opt.V) = -10;
                    opt.bound(2,1:opt.V) = 10;
                    
                case 'g10'
                    opt.V = 8; 
                    opt.C = 6;
                    opt.bound = zeros(2, opt.V);
                    opt.bound(1,1)  = 100;
                    opt.bound(2,1)  = 10000;
                    opt.bound(1,2:3)  = 1000;
                    opt.bound(2,2:3)  = 10000;
                    opt.bound(1,4:8)  = 10;
                    opt.bound(2,4:8)  = 1000;

            end
            opt.totalFuncEval = 500;%opt.V*10+21;%high-fidelity function evaluation
            opt.initpopsize = opt.V*10;%initial sample size for high fidelity computation
            
        case {'ddmop1','ddmop2','ddmop3','ddmop4','ddmop5','ddmop6','ddmop7'}
            functionHandle  = {@DDMOP1, @DDMOP2, @DDMOP3, @DDMOP4, @DDMOP5, @DDMOP6, @DDMOP7};  % test problems
            PopSize    = [256,105,105,256,105,100,100];	% population size
            ObjDim   = [9,3,3,10,3,2,2];				% number of objectives
            NumVar    = [11,5,6,13,15,10,17];			% number of decision variables
            FuncEval = [400,300,400,600,800,300,600];	% total number of real objective function evaluations
            opt.C = 1;   
            switch(opt.objfunction)
                case 'DDMOP1'
                    opt.V = NumVar(1);
                    opt.totalFuncEval = FuncEval(1);
                    opt.M = ObjDim(1);
                    opt.initpopsize = PopSize(1);
                    opt.functionHandle = functionHandle{1};
                case 'DDMOP2'
                    opt.V = NumVar(2);
                    opt.totalFuncEval = FuncEval(2);
                    opt.M = ObjDim(2);
                    opt.initpopsize = PopSize(2);
                    opt.functionHandle = functionHandle{2};
                case 'DDMOP3'
                    opt.V = NumVar(3);
                    opt.totalFuncEval = FuncEval(3);
                    opt.M = ObjDim(3);
                    opt.initpopsize = PopSize(3);
                    opt.functionHandle = functionHandle{3};
                case 'DDMOP4'
                    opt.V = NumVar(4);
                    opt.totalFuncEval = FuncEval(4);
                    opt.M = ObjDim(4);
                    opt.initpopsize = PopSize(4);
                    opt.functionHandle = functionHandle{4};
                case 'DDMOP5'
                    opt.V = NumVar(5);
                    opt.totalFuncEval = FuncEval(5);
                    opt.M = ObjDim(5);
                    opt.initpopsize = PopSize(5);
                    opt.functionHandle = functionHandle{5};
                case 'DDMOP6'
                    opt.V = NumVar(6);
                    opt.totalFuncEval = FuncEval(6);
                    opt.M = ObjDim(6);
                    opt.initpopsize = PopSize(6);
                    opt.functionHandle = functionHandle{6};
                case 'DDMOP7'
                    opt.V = NumVar(7);
                    opt.totalFuncEval = FuncEval(7);
                    opt.M = ObjDim(7);
                    opt.initpopsize = PopSize(7);
                    opt.functionHandle = functionHandle{7};
            end
            [lowerval, upperval] = DDMOP1('boundary');  		%Output the lower and upper bounary of the decision variables
            opt.bound(1,:) = lowerval;
            opt.bound(2,:) = upperval; 
             	%Evaluate the decision vectors and output the objective vectors
        case {'UF1'}
            
        otherwise
            input('Function definition is not found');
    end
    %======================================================================
    %===================REFERENCE DIRECTIONS/POINTS========================
    %======================================================================
    
    if opt.M>1
        %[opt.dirs, opt.numdir] = UniformPoint(opt.N, opt.M);
        if opt.M==15
            opt.dirs = initweight(15, 135)';
        elseif opt.M==10
            opt.dirs = initweight(10, 275)';
        elseif opt.M==8
            opt.dirs = initweight(5, 156)';
    %         p = 4;
    %         H1 = initweight(5, nchoosek(opt.M+p-1,p))';
    %         p = 3;
    %         H_temp = initweight(5, nchoosek(opt.M+p-1,p))';
    %         H2 = layered_weight(.85, H1);
    %         H3 = layered_weight(.6, H_temp);
    %         H4 = layered_weight(.3, H_temp);
    %         opt.dirs = vertcat(H1, H2, H3, H4);

        elseif opt.M==3
            opt.dirs = initweight(3, 91)';
            %opt.refdiv = 12;
            %{
            p = 2*opt.M-1;
            H = cell(1,opt.M);
            H{1} = initweight(opt.M, nchoosek(opt.M+p-1,p))';
            opt.dirs = H{1};
            alphas = 0:(1/opt.M):1;
            for i=2:opt.M
                H{i} = layered_weight(alphas(i), H{1});
                opt.dirs = vertcat(opt.dirs, H{i});
            end

            opt.dirs = unique(opt.dirs,'rows');
            %}
            %{
            figure;
            hold all;
            plot3(H{1}(:,1),H{1}(:,2),H{1}(:,3),'bo');
            plot3(H{2}(:,1),H{2}(:,2),H{2}(:,3),'ro');
            plot3(H{3}(:,1),H{3}(:,2),H{3}(:,3),'go');
            %}

        elseif  opt.M==5 %Five obj
            opt.dirs = initweight(5, 210)';        
        elseif opt.M==2
            opt.dirs = initweight(2, 21)';
        end
        opt.numdir = size(opt.dirs,1);%number of reference direction
    else
        opt.dirs = 1;%linspace(0,1, 21)';
        opt.numdir = 1;%21;
    end
    
    opt.dirs(opt.dirs<1e-16) = 1e-16;     
    opt.curdir = opt.dirs(1,:);%current direction
    opt.curcluster = 1;%current cluster number
    

    %======================================================================
    %==================VARIABLE MEMORY ALLOCATION==========================
    %======================================================================
    
    opt.archive = [];
    opt.archiveObj = [];
    opt.archiveCV = [];
    opt.archiveASF = [];
    opt.archiveACV = [];
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
    
    if opt.C>0
        opt.regC = cell(1,opt.C);
    end
    opt.regO = cell(1,opt.M);
    
    %======================================================================
    %====================OTHER PARAMETERS==================================
    %======================================================================
    
    opt.InitialActiveSetSize = 200;
    opt.activeSetSize = opt.InitialActiveSetSize;
    opt.pmut = 1.0/opt.V; % Mutation probability

    
    %======================================================================
    %===============LOW-FIDELITY NSGA-II PARAMETERS========================
    %======================================================================
    
    if opt.M<=2
        opt.methodology12_option=1;%1=NSGA2 for opt.M=2, 2=NSGA3 for opt.M>3
    else
        opt.methodology12_option=2;%1=NSGA2 for opt.M=2, 2=NSGA3 for opt.M>3
    end
    opt.nsga2.totalFuncEval = opt.totalFuncEval;%600;
    opt.nsga2.eta_c = 20;%crossover index
    opt.nsga2.eta_m = 20;%mutation index
    opt.nsga2.N = 100;%50;%population size in optimization algorithm
    
    if opt.N < opt.numdir %&& opt.M>1
        opt.N = opt.numdir;
%     elseif opt.N < opt.numdir && opt.M==1
%         opt.N = 10;
%         opt.nsga2.N = 10;
    end
    
    if opt.M>2
        opt.nsga2.N = opt.numdir;
        %opt.N = opt.numdir;
    end
    opt.nsga2.initpopsize = opt.nsga2.N;
    opt.nsga2.CD = zeros(opt.nsga2.N,1);%initial crowding distance
    opt.nsga2.G = 300;
    opt.nsga2.pcross = 0.9; % Crossover probability
    opt.nsga2.nrealcross = 0;%number of crossover performed
    opt.nsga2.nrealmut = 0;%number of mutation performed
    opt.nsga2.gen = 1;
    opt.nsga2.pop = [];
    opt.nsga2.popObj = [];
    opt.nsga2.Epsilon = 1e-16;
    opt.nsga2.Inf = 1e14;
    opt.nsga2.crossoverOption = 1;% 1 = simulated binary crossover
    opt.nsga2.mutationOption = 1;% 1 = polynomial mutation
    opt.nsga2.matingselectionOption = 1;%1 = binary constraint tournament selection, 2 = nsga3 selection, 3 = tournament with Knee based
    opt.nsga2.survivalselectionOption = 1;%1 = NSGA-II, 2 = NSGA-III, 3 = PageRank
    opt.nsga2.pmut = opt.pmut;
    opt.nsga2.nadir_point = [];
    opt.nsga2.ideal_point = repmat(realmax,1,opt.M);
    %======================================================================
    %===============SPECIFIC NSGA-III Parameters===========================
    %======================================================================
    
    opt.nsga3.dirs = opt.dirs;
    opt.nsga3.numdir = size(opt.nsga3.dirs,1);
    opt.nsga3.associationsReady = false;
    
    %======================================================================
    %=======================TRUST REGION PARAMETERS========================
    %======================================================================
    if opt.func_family_no==3
        opt.trust_region_option_active = 2;
    else
        opt.trust_region_option_active = 1;
    end
    opt.trust_region_update_option = 1;%1== continuous decreasing, 2 = non-linear/diffusion trust region, 3 - adaptive
    %opt.trust_region_option_active = 1;%1== on, 2 = off
    opt.adaptive_trust_region_option = 1;%1 = HV, 2=ASF
    opt.delta = 1;%sqrt(opt.V);%1.0;%0.5;%1.2;%1.0;%per variable
    %opt.TrustDistObj = sqrt(opt.M)/opt.refdiv;%repmat(sqrt(opt.M)/opt.refdiv, [1 opt.numdir]);
    %opt.TrustDistVar = sqrt(opt.V)/opt.refdiv;%repmat(sqrt(opt.V)/opt.refdiv, [1 opt.V]);
    %opt.TrustDistObj = 0.1*sqrt(sum((opt.max_val - opt.min_val).^2))/opt.refdiv;
    %sqrt(opt.V);%0.5;%sqrt(sum((opt.bound(2,:) - opt.bound(1,:)).^2))/opt.refdiv;
    %opt.TrustDistLearningRate = 0.05;% 10% grow or shrink after each hi-fi evaluations
    
    
    %continuous decrease
    opt.TrustDistVar =  sqrt(opt.V);%norm(opt.bound(2,:)-opt.bound(1,:));
    opt.TrustDistVarMin = 0.01*opt.TrustDistVar;
    opt.maxIter = floor(opt.totalFuncEval/opt.numdir);
    opt.iter = 1;
    opt.TrustVarMaxLR = 0.9;
    opt.TrustVarMinLR = 0.1;
    opt.TrustVarLR = opt.TrustVarMaxLR;
    if strcmpi(opt.objfunction,'zdt4')
        opt.TrustDistVarMin = 0.5*opt.TrustDistVar;%0.4*opt.TrustDistVar;
        opt.TrustVarLR = 0.95;
    end
     
    %nonlinear trust region
    opt.TrustDelta =  sqrt(opt.V);%norm(opt.bound(2,:)-opt.bound(1,:));
    opt.TrustDeltaMin = 0.1* opt.TrustDelta;
    opt.TrustMaxDeltaLR = 0.95;
    opt.TrustMinDeltaLR = 0.1;
    opt.TrustDeltaLR = opt.TrustMaxDeltaLR;
    opt.changeDelta = 200;%change delta in every ## function evaluation
    opt.epsilonTrust = 100;%controls smoothness of trusted region, larger means isolated
    opt.SampleSize = 10;%how many points are minimum needed around an expensively evaluated point to build the model
    opt.NumExplorationPoints = 3;
    opt.TrustDeltaFileName = strcat('methodology',num2str(opt.methodology),'_',opt.objfunction,'_delta.txt');
    
    %adaptive trust-region
    opt.TrustInit = sqrt(opt.V);
    opt.TrustRadiusDeltaK = repmat(opt.TrustInit, 1, opt.initpopsize);
    %opt.TrustLinearRate = opt.TrustRadiusDeltaK(1)/opt.TrustDistVar;
    opt.TrustDeltaStar = sqrt(opt.V);
    opt.TrustC1 = 0.75;%0.75;
    opt.TrustC2 = 1.10;
    opt.TrustR1 = 0.9;%0.1;
    opt.TrustR2 = 1.05;%0.75;
    
    
    %======================================================================
    %=================METAMODEL PARAMETERS=================================
    %======================================================================
    
    opt.obj_metamodel_list = {'GP','RBF','SVR','NN'};
    opt.cons_metamodel_list = {'GP'};
    opt.obj_aggregation_metamodel_list = {'GP'};
    opt.cons_aggregation_metamodel_list = {'GP'};
    opt.cons_aggregation_option_list = {'CV','ACV','ELCV','AELCV','TANH','ATANH'};
    
    opt.aggregated_obj = {{1,2}, {3,4,5}};
    opt.aggregated_cons = {{1}};%{{1,2}, {3,4,5}};
    opt.aggregated_obj_cons = {{{1,2},1}, {{3,4,5},2}};    
    
    opt.current_obj_metamodel = repmat({'GP'}, 1, opt.M);
    opt.current_cons_metamodel = repmat({'GP'}, 1, opt.M);
    opt.current_obj_aggregation_metamodel = {'GP'};
    opt.current_cons_aggregation_metamodel = {'GP'};
    opt.current_m5_metamodel = {'GP'};
    opt.current_m6_metamodel = {'GP'};%{'NN'}
    
    opt.current_obj_aggregation_option = {'ASF'};
    opt.current_cons_aggregation_option = {'CV'};%{'CV'};
    opt.current_selection_function_option_m5 = {'ASFCV'};
    opt.current_selection_function_option_m6 = {'ASFCV'};%{'MEMO'};
    
    
    %======================================================================
    %=================LOW FIDELITY OPTIONS=================================
    %======================================================================
    
    opt.generative_framework_acquisition_func = {'ASF'};
    
    %Gaussian Process Model Parameters
    opt.theta = 10;% starting value of Kriging correlation parameter
    opt.lob = 1e-2;%lower bound of Kriging correlation parameter
    opt.upb = 20;%20;%upper bound of Kriging correlation parameter
    opt.consOption = 1; %1  = CV, 2 = mCV, option for modeling constraints
    opt.rho = 0;%0.001;%----rho value for augment asf MEMO approach

    
    %======================================================================
    %======================OUTPUT OPTIONS==================================
    %======================================================================
    
    opt.plotOption = 2;%1 = on, 2 = off
    if opt.plotOption==1
        opt.fig = figure;
    end
    opt.writeFlag = 1;%1 = write information to file, 2 = don't write

    opt.varfilename = strcat('BenchmarkResults/',upper(opt.test_func_family{opt.func_family_no}),'/',opt.algo_name,'_',lower(opt.objfunction),'_varfile_',   num2str(opt.V),'_',num2str(opt.M), '_', num2str(opt.r),'.txt');%save variables
    opt.objfilename = strcat('BenchmarkResults/',upper(opt.test_func_family{opt.func_family_no}),'/',opt.algo_name,'_',lower(opt.objfunction),'_objfile_',   num2str(opt.V),'_',num2str(opt.M), '_', num2str(opt.r),'.txt');%save objectives
    opt.cvfilename  = strcat('BenchmarkResults/',upper(opt.test_func_family{opt.func_family_no}),'/',opt.algo_name,'_',lower(opt.objfunction),'_cvfile_',    num2str(opt.V),'_',num2str(opt.M), '_', num2str(opt.r),'.txt');%save constraint violation
    opt.consfilename  = strcat('BenchmarkResults/',upper(opt.test_func_family{opt.func_family_no}),'/',opt.algo_name,'_',lower(opt.objfunction),'_consfile_',num2str(opt.V),'_',num2str(opt.M), '_', num2str(opt.r),'.txt');%save constraint values
    
    dlmwrite(opt.varfilename, [], 'delimiter',' ','precision','%.10f');
    dlmwrite(opt.objfilename, [], 'delimiter',' ','precision','%.10f');
    dlmwrite(opt.cvfilename, [], 'delimiter',' ','precision','%.10f');
    dlmwrite(opt.consfilename, [], 'delimiter',' ','precision','%.10f');
    
    if opt.num_of_frameworks > 1
        opt.sep_filename = strcat('BenchmarkResults/',upper(opt.test_func_family{opt.func_family_no}),'/',opt.algo_name,'_error_data_sep_',lower(opt.objfunction),'_',num2str(opt.V),'_',num2str(opt.M), '_', num2str(opt.r),'.txt');
        opt.mse_filename = strcat('BenchmarkResults/',upper(opt.test_func_family{opt.func_family_no}),'/',opt.algo_name,'_error_data_mse_',lower(opt.objfunction),'_',num2str(opt.V),'_',num2str(opt.M), '_', num2str(opt.r),'.txt');
        opt.selected_method_filename = strcat('BenchmarkResults/',upper(opt.test_func_family{opt.func_family_no}),'/',opt.algo_name,'_selected_method_',lower(opt.objfunction),'_',num2str(opt.V),'_',num2str(opt.M), '_', num2str(opt.r),'.txt'); 
        dlmwrite(opt.sep_filename, [], 'delimiter',' ','precision','%.10f');
        dlmwrite(opt.mse_filename, [], 'delimiter',' ','precision','%.10f');
        dlmwrite(opt.selected_method_filename, [], 'delimiter',' ');
    end
    
    %---------------INTERFACE WITH PLATEMO---------------------------------
    opt.Global.M = opt.M;
    opt.Global.D = opt.V;
    opt.Global.lower = opt.bound(1,:);
    opt.Global.upper = opt.bound(2,:);
    opt.PF = feval(str2func(upper(opt.objfunction)), 'PF', opt.Global, opt.numdir);%pareto front for IGD
    opt.PFGD = feval(str2func(upper(opt.objfunction)), 'PF', opt.Global, 10*opt.numdir);
    %temp_objfunction = regexprep(opt.objfunction, '_', '' );
    %opt.PF = load(strcat(upper(temp_objfunction),'.', num2str(opt.M),'D.pf'));%Pareto front for IGD calculation
    %opt.PFGD = load(strcat('GD/',upper(temp_objfunction),'.', num2str(opt.M),'D.pf'));%Pareto front for GD calculation   
    
    opt.min_val = min(opt.PF);%problem information
    opt.max_val = max(opt.PF);%problem information
    opt.utopian = opt.bound(1,:)-opt.Epsilon;
    
    %----------------------------------------------------------------------
    
        
end

%------------------------------END OF -FILE--------------------------------

