function main()
    
    clear all;
    close all;
    test_function = {{'DTLZ1','DTLZ2','DTLZ3', 'DTLZ4', 'SDTLZ1', 'SDTLZ2','CDTLZ2'...
                      'IDTLZ1','IDTLZ2','C2_DTLZ2', 'C3_DTLZ4', 'C1_DTLZ1', ...
                       'DTLZ5',  'DTLZ7', 'DTLZ5IM', 'DTLZ8',  'DTLZ6', 'DTLZ9',},...
                     {'ZDT1', 'ZDT2','ZDT3','ZDT4','ZDT6'}...
                     {'G1','G2','G3','G4','G5','G6','G7','G8','G9','G10'},...
                     {'TNK','OSY','BNH','SRN'},...
                     {'DDMOP1','DDMOP2','DDMOP3','DDMOP4','DDMOP5','DDMOP6','DDMOP7'}};
    
    test_func_family = {'DTLZ', 'ZDT', 'G', 'PRACTICAL', 'DDMOP', 'CF', 'UF', 'WFG'};
    
    
%-------------------MAIN LOOP----------------------------------------------
    method_list = [6];
    
    for methodology = 1:1:1%for all methodologies
        
        for func_family_no=2:2   

            %addpath(strcat('Problems/',test_func_family{func_family_no}));
            addpath(genpath('Problems'));
            addpath('Metrics');
        
            for func_no = 1:1:1%size(test_function{func_family_no},2)%for test function no in that family

                for d=1:1:1%extra running parameters for DTLZ problems

                    for r = 1:1:1%number of runs
                        
                        opt.methodology = method_list(methodology);
                        opt.func_no = func_no;%function number
                        opt.func_family_no = func_family_no;
                        opt.test_function = test_function;
                        opt.test_func_family = test_func_family;
                        opt.dim = d;%number of objectives
                        opt.r = r;%run number
                        opt.objfunction = lower(strtrim(opt.test_function{opt.func_family_no}{opt.func_no}));%remove whitespaces
                        opt = basic_parameters(opt);%all parameters of frameworks
                        disp(['Problem: ' opt.objfunction, ', Framework: ' num2str(opt.methodology), ', Run: ' num2str(opt.r)  ', Obj: ' num2str(opt.M) ...
                            ', Var: ', num2str(opt.V)]);
                        %---------------- OPTIMIZE ------------------------
                        opt = parallel_metamodel_main_loop(opt);
                        %-------------COMPUTE PERFORMANCE METRIC-----------
                        opt = compute_performance_metric(opt, opt.archiveObj, opt.archiveCV);                    
                        clear opt;

                    end
                end
            end
        end
    end
    disp('============END OF ALL RUNS==============');
end


%------------------------------END OF -FILE--------------------------------
