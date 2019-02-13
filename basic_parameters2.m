% Copyright [2018] [Proteek Chandan Roy]
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

%
% Proteek Chandan Roy, 2018
% Contact: royprote@egr.msu.edu, proteek_buet@yahoo.com

function [opt] = basic_parameters2(opt)

    
    %======================================================================
    %==================OPTIMIZATION OPTIONS================================
    %======================================================================
    opt.writeFlag = 2;
    opt.initOption = 1;%1 = latin hypercube on decision space, 2 = roulette wheel based on objective space, add more options
    opt.metamodelOption = 1;% 1 = solutions will be evaluated by metamodel, 2 = high-fidelity
    opt.crossoverOption = 1;% 1 = simulated binary crossover
    opt.mutationOption = 1;% 1 = polynomial mutation
    opt.matingselectionOption = 1;%1 = binary constraint tournament selection
    opt.survivalselectionOption = 1;%1 = NSGA-II, 2 = NSGA-III
    opt.selection_function_option = 2;%1=ASFCV, 2=MEMO, 3=NSGA-II
    opt.switching_option = 1;
    
    opt.plotOption = 2;%plot option
    if opt.plotOption==1
        opt.fig = figure;
    end
    
    opt.switching_option = 4;%1=swithcing part1, 2= part2, 3=part3, 4=no switch
    delete(gcp('nocreate'));
    
    if opt.switching_option==1
        opt.num_of_frameworks = 6;
        %parpool(min(feature('numcores'),6));
    elseif opt.switching_option==2
        opt.num_of_frameworks = 5;
        %parpool(min(feature('numcores'),5));
    elseif opt.switching_option==3
        opt.num_of_frameworks=10;
        %parpool(min(feature('numcores'),10));
    else
        opt.num_of_frameworks=1;
    end
    opt.total_num_of_frameworks = 10;
    opt.writeFlagOption = 1;
    
    

    %======================================================================
    %=================METAMODEL PARAMETERS=================================
    %======================================================================
    
    opt.theta = 10;% starting value of Kriging correlation parameter
    opt.lob = 1e-2;%lower bound of Kriging correlation parameter
    opt.upb = 20;%20;%upper bound of Kriging correlation parameter
    opt.consOption = 1; %1  = CV, 2 = mCV, option for modeling constraints
    opt.rho = 0;%0.001;%----rho value for augment asf MEMO approach
    %opt.metamodel = 1%1=Kriging Model
    opt.funcEval = 0;
    opt.iter = 1;
    
    
    %======================================================================
    %==============BASIC RGA ALGORITHM PARAMETERS==========================
    %======================================================================
     
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
    opt.Epsilon = 1e-14;
    
    
    %======================================================================
    %=========================FILE NAMES FOR SAVE==========================
    %======================================================================
    
    if opt.switching_option == 4
        opt.varfilename = strcat('methodology',num2str(opt.methodology),'_',opt.objfunction,'_var_',num2str(opt.r),'.txt');
        opt.objfilename = strcat('methodology',num2str(opt.methodology),'_',opt.objfunction,'_obj_',num2str(opt.r),'.txt');
        opt.cvfilename  = strcat('methodology',num2str(opt.methodology),'_',opt.objfunction,'_cv_',num2str(opt.r),'.txt');
    else
        opt.varfilename = strcat('switching_framework_',opt.objfunction,'_var_',num2str(opt.r),'.txt');
        opt.objfilename = strcat('switching_framework_',opt.objfunction,'_obj_',num2str(opt.r),'.txt');
        opt.cvfilename  = strcat('switching_framework_',opt.objfunction,'_cv_',num2str(opt.r),'.txt');
        opt.sep_filename = strcat('error_data_sep_',lower(opt.objfunction),'_',num2str(opt.r),'.txt');
        opt.mse_filename = strcat('error_data_mse_',lower(opt.objfunction),'_',num2str(opt.r),'.txt');
        opt.selected_method_filename = strcat('selected_method_',lower(opt.objfunction),'_',num2str(opt.r),'.txt'); 
        
        dlmwrite(opt.sep_filename, [], 'delimiter',' ','precision','%.10f');
        dlmwrite(opt.mse_filename, [], 'delimiter',' ','precision','%.10f');
        dlmwrite(opt.selected_method_filename, [], 'delimiter',' ');
    end
       
    
    %======================================================================
    %==================TEST PROBLEMS PARAMETERS============================
    %======================================================================
    opt.nadir_point = [];
    opt.ideal_point = [];

    switch(opt.objfunction)
        case {'zdt1','zdt2','zdt3','zdt4','zdt6'}
            opt.M = 2;%number of objectives
            opt.V = 10;%;10;%number of variables
            opt.C = 0;%number of constraints
            opt.totalFuncEval = 121;%500;%high-fidelity function evaluation
            opt.utopian = [-0.05, -0.05];%ideal point
            opt.min_val = [0 0];%minimum value for normalization
            opt.max_val = [1 1];%maximum objective 
            opt.initpopsize = 100;%initial sample size for high fidelity computation
            if strcmp(opt.objfunction,'zdt6')
                opt.utopian = [-0.05, -0.05];
                opt.min_val = [0.25 0];
                opt.max_val = [1 1];
                opt.V = 10;
                %opt.initpopsize = 200;
            elseif strcmp(opt.objfunction,'zdt3')
                opt.utopian = [-0.05, -1.1];
                opt.min_val = [0 -1];
                opt.max_val = [1 1];
            elseif strcmp(opt.objfunction,'zdt4')
                opt.V = 5;
                %opt.totalFuncEval = 1000;
                opt.utopian = [-0.05, -0.05];
            end
            
        case {'dtlz1', 'dtlz2', 'dtlz3','dtlz4', 'dtlz5', 'dtlz7'}
            opt.M = 5;
            opt.V = 7;%10; 
            opt.C = 0;
            opt.totalFuncEval = 2000;
            opt.initpopsize = 700;
                
%             if strcmp(opt.objfunction,'dtlz4')
%                 opt.initpopsize = 700;%initial sample size for high fidelity computation
%                 opt.totalFuncEval = 2000;
%             elseif strcmp(opt.objfunction,'dtlz7')
%                 opt.M = 3;
%             end
            opt.utopian = (-0.05).*ones(1, opt.M);
            opt.min_val = zeros(1,opt.M);
            opt.max_val = ones(1, opt.M);
            
        case {'wfg1', 'wfg2', 'wfg3','wfg4', 'wfg5', 'wfg6','wfg7','wfg8','wfg9'}
            opt.M = 3;
            opt.V = 8; 
            opt.C = 0;
            opt.totalFuncEval = 500;
            opt.utopian = (-0.05).*ones(1, opt.M);
            opt.min_val = zeros(1,opt.M);
            opt.max_val = ones(1, opt.M);
            opt.initpopsize = 100;
            if strcmp(opt.objfunction,'wfg4')
                opt.initpopsize = 100;%initial sample size for high fidelity computation
                opt.totalFuncEval = 500;
            end
            
        case {'c2dtlz2','c3dtlz2'}
            opt.M = 5;
            opt.V = 7;
            opt.totalFuncEval = 2000;
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
            opt.totalFuncEval = 500;
            opt.utopian = [-0.05, -0.05];
            opt.min_val = [0 0];
            opt.max_val = [140 55];
            opt.initpopsize = 100;
            
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
            opt.totalFuncEval = 500;
            opt.utopian = [0 -300];
            opt.min_val = [0 -250];
            opt.max_val = [240 0];
            opt.initpopsize = 100;
        case 'tnk'
            opt.M = 2;
            opt.V = 2; 
            opt.C = 2;
            opt.totalFuncEval = 500;
            opt.utopian = [-0.001, -0.001];
            opt.min_val = [0 0];
            opt.max_val = [1.2 1.2];
            opt.initpopsize = 200;
        case 'water'
            opt.M = 5;
            opt.V = 3; 
            opt.C = 7;
            opt.totalFuncEval = 500;
            opt.utopian = (-0.05).*ones(1, opt.M);
            opt.min_val = [0.75 0 0 0 0];
            opt.max_val = [0.95 0.9 1.0 1.6 3.2];
            opt.initpopsize = 100;
        case 'carside'
            opt.M = 3;
            opt.V = 7; 
            opt.C = 10;
            opt.totalFuncEval = 500;
            opt.utopian = [24.3180    3.5352   10.5610];
            opt.min_val = [24.3680    3.5852   10.6110];
            opt.max_val = [42.7620    4.0000   12.5210];
            opt.initpopsize = 100;
            
        case 'welded'
            opt.M = 2;
            opt.V = 4; 
            opt.C = 4;
            opt.totalFuncEval = 500;
            opt.utopian = [2.3316   -0.0496];
            opt.min_val = [2.3816    0.0004];
            opt.max_val = [36.4403    0.0157];
            opt.initpopsize = 100;
        case 'toy'
            opt.M = 2;
            opt.V = 1;
            opt.C = 0;
            opt.totalFuncEval = 200;
            opt.utopian = [0    0];
            opt.min_val = [0    0];
            opt.max_val = [1.2    1.2];
            opt.initpopsize = 10;
            
        otherwise
            input('Function definition is not found');
    end
    
    %======================================================================
    %==============DECISION VARIABLE BOUNDS================================
    %======================================================================
    
    
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
    
    %======================================================================
    %==================LOAD PARETO FRONT===================================
    %======================================================================
    
    
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
    
    
    
    
    %======================================================================
    %===================REFERENCE DIRECTIONS/POINTS========================
    %======================================================================
    
    if opt.M==15
        opt.dirs = initweight(15, 135)';
        opt.refdiv = 6;
    elseif opt.M==10
        opt.dirs = initweight(10, 275)';
        opt.refdiv = 6;
    elseif opt.M==8
        opt.dirs = initweight(5, 156)';
        opt.refdiv = 6;
%     elseif opt.M==5
%         %opt.dirs = initweight(5, 126)';%initweight(5, 210)';
%         %opt.refdiv = 6;
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
        %{
        p = 15;
        H1 = initweight(3, nchoosek(opt.M+p-1,p))';
        p = 15;
        H_temp = initweight(3, nchoosek(opt.M+p-1,p))';
        H2 = layered_weight(.85, H1);
        H3 = layered_weight(.6, H_temp);
        H4 = layered_weight(.3, H_temp);
        opt.dirs = vertcat(H1, H2, H3, H4);
        
        figure;
        hold all;
        plot3(H1(:,1),H1(:,2),H1(:,3),'bo');
        plot3(H2(:,1),H2(:,2),H2(:,3),'ro');
        plot3(H3(:,1),H3(:,2),H3(:,3),'go');
        plot3(H4(:,1),H4(:,2),H4(:,3),'co');
        %}
        %H1 = initweight(opt.M, nchoosek(opt.M+4-1,4))';
        %H1_temp = initweight(opt.M, nchoosek(opt.M+2-1,2))';
        %H2 = layered_weight(0.5, H1);
        %H3 = layered_weight(0.3, H1_temp);%ones(1, opt.M)*(1/opt.M);
        %opt.dirs = vertcat(H1, H2);
        %opt.dirs = unique(opt.dirs,'rows');
    elseif  opt.M==5 %Five obj
        opt.dirs = initweight(5, 210)';
%         p = 15;
%         H1 = initweight(5, nchoosek(opt.M+1-1,1))';%5
%         H1_temp = initweight(5, nchoosek(opt.M+2-1,2))';%15
%         H2 = layered_weight(.75, H1_temp);%15
%         H3 = layered_weight(.50, H1_temp);%15
%         H4 = layered_weight(.35, H1_temp);%15
%         opt.dirs = vertcat(H1, H2, H3, H4);
        
    elseif opt.M==2
        opt.dirs = initweight(2, 21)';
        opt.refdiv = 20;
    end

    opt.numdir = size(opt.dirs,1);%number of reference direction
    opt.curdir = opt.dirs(1,:);%current direction
    opt.curcluster = 1;%current cluster number
    opt.pmut = 1.0/opt.V; % Mutation probability
    opt.ref_point = cell(1, opt.numdir);

    %======================================================================
    %==================MEMORY ALLOCATION===================================
    %======================================================================
    
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
    opt.InitialActiveSetSize = 200; %opt.initpopsize;
    opt.activeSetSize = opt.InitialActiveSetSize;
    if opt.C>0
        opt.regC = cell(1,opt.C);
    end
    opt.regO = cell(1,opt.M);

    %======================================================================
    %==================NSGA-II PARAMETERS==================================
    %======================================================================
    
    if opt.M==2
        opt.methodology12_option=1;%1=NSGA2 for opt.M=2, 2=NSGA3 for opt.M>3
    else
        opt.methodology12_option=2;%1=NSGA2 for opt.M=2, 2=NSGA3 for opt.M>3
    end
    opt.nsga2.totalFuncEval = opt.totalFuncEval;%600;
    opt.nsga2.nonlinearTrust = opt.nsga2.totalFuncEval;%800;
    opt.nsga2.buildModelIntervalGeneration = 100;
    opt.nsga2.eta_c = 15;%crossover index
    opt.nsga2.eta_m = 20;%mutation index
    opt.nsga2.N = 100;%50;%population size in optimization algorithm
    
    if opt.N < opt.numdir 
        opt.N = opt.numdir;
    end
    %{
    iter = 2;
    while opt.N < opt.numdir 
        opt.N = 100*iter;
        iter = iter+1;
    end
    %}
    
    if opt.M>2
        opt.nsga2.N = opt.numdir;
        %opt.N = opt.numdir;
    end
    opt.nsga2.initpopsize = opt.nsga2.N;
    opt.nsga2.CD = zeros(opt.nsga2.N,1);%initial crowding distance
    opt.nsga2.G = 300;
    %opt.nsga2.buildModelInterval = opt.nsga2.N*opt.nsga2.buildModelIntervalGeneration;
    opt.nsga2.pcross = 0.9; % Crossover probability
    opt.nsga2.nrealcross = 0;%number of crossover performed
    opt.nsga2.nrealmut = 0;%number of mutation performed
    opt.nsga2.gen = 1;
    opt.nsga2.pop = [];
    opt.nsga2.popObj = [];
    opt.nsga2.Epsilon = 1e-14;
    opt.nsga2.Inf = 1e14;
    opt.nsga2.crossoverOption = 1;% 1 = simulated binary crossover
    opt.nsga2.mutationOption = 1;% 1 = polynomial mutation
    opt.nsga2.matingselectionOption = 1;%1 = binary constraint tournament selection, 2 = nsga3 selection, 3 = tournament with Knee based
    opt.nsga2.survivalselectionOption = 1;%1 = NSGA-II, 2 = NSGA-III, 3 = PageRank
    opt.nsga2.pmut = opt.pmut;
    opt.nsga2.nadir_point = [];
    opt.nsga2.ideal_point = repmat(realmax,1,opt.M);
    %======================================================================
    %======================NSGA-III Parameters=============================
    %======================================================================
    
%     dirCount = inf;
%     tempN = opt.nsga2.N;
%     while dirCount > opt.nsga2.N
%         opt.nsga3.dirs = initweight(opt.M, tempN)';
%         opt.dirs = opt.nsga3.dirs;
%         dirCount = size(opt.dirs, 1);
%         tempN = tempN - 1;
%     end
    opt.nsga3.dirs = opt.dirs;
    opt.nsga3.numdir = size(opt.nsga3.dirs,1);
    opt.nsga3.associationsReady = false;
    
    %======================================================================
    %=======================TRUST REGION OPTIONS===========================
    %======================================================================
    
    opt.trust_region_update_option = 1;%1== continuous decreasing, 2 = non-linear/diffusion trust region, 3 - adaptive
    opt.trust_region_option_active = 1;%1== on, 2 = off
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
    
    
    %---------------INTERFACE WITH PLATEMO---------------------------------
    opt.Global.M = opt.M;
    opt.Global.D = opt.V;
    opt.Global.lower = opt.bound(1,:);
    opt.Global.upper = opt.bound(2,:);
    opt.PF = feval(str2func(upper(opt.objfunction)), 'PF', opt.Global, opt.numdir);%pareto front for IGD
    opt.PFGD = feval(str2func(upper(opt.objfunction)), 'PF', opt.Global, 10*opt.numdir);
    opt.min_val = min(opt.PF);%problem information
    opt.max_val = max(opt.PF);%problem information
    opt.utopian = opt.bound(1,:)-opt.Epsilon;
    
    
end

%------------------------------END OF -FILE--------------------------------

