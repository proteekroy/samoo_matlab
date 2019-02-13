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

% This Software runs NSGA-II and its variants for different testfunctions
% You are free to change, modify, share and distribute this software
% Please Acknowledge the author (if possible) for any use of this software
% @author: Proteek Chandan Roy, Department of CSE
% Michigan State University, USA
% email: proteek_buet@yahoo.com, royprote@egr.msu.edu 

function igd_for_single_objective_optimization()
    
    clear all;
    close all;
    test_function = {{'DTLZ1','DTLZ2','DTLZ3', 'DTLZ4', 'SDTLZ1', 'SDTLZ2','CDTLZ2'...
                      'IDTLZ1','IDTLZ2','C2_DTLZ2', 'C3_DTLZ4', 'C1_DTLZ1', 'C1_DTLZ3'...
                       'DTLZ5',  'DTLZ7', 'DTLZ5IM', 'DTLZ8',  'DTLZ6', 'DTLZ9',},...
                     {'ZDT1', 'ZDT2','ZDT3','ZDT4','ZDT6'},...
                     {'G1','G2','G3','G4','G5','G6','G7','G8','G9','G10'},...
                     {'TNK','OSY','BNH','SRN'}};
    
    test_func_family = {'DTLZ', 'ZDT', 'G', 'PRACTICAL', 'CF', 'UF', 'WFG'};
    
    
%-------------------MAIN LOOP----------------------------------------------
    for func_family_no=3:3   
        
        addpath(strcat('Problems/',test_func_family{func_family_no}));
        addpath('Metrics');
        
        for func_no = 1:1:size(test_function{func_family_no},2)%for test functions
            
            for d=2:1:2%extra running parameters for DTLZ problems
                
                for r = 1:1:1%number of runs
                    opt.methodology = 11;
                    opt.func_no = func_no;%function number
                    opt.func_family_no = func_family_no;
                    opt.test_function = test_function;
                    opt.test_func_family = test_func_family;
                    opt.dim = d;%number of objectives
                    opt.r = r;%run number
                    opt.objfunction = lower(strtrim(opt.test_function{opt.func_family_no}{opt.func_no}));%remove whitespaces
                    %disp(opt.objfunction);
                    opt = basic_parameters(opt);%basic parameters MOEA
                    opt.archiveObj = load(opt.objfilename);
                    opt.archiveCV = load(opt.cvfilename);
%                     temp_popObj = min(opt.archiveObj(opt.archiveCV<=0, :));
%                     opt.igd = IGD(temp_popObj, opt.PF);
%                     opt.hv = HV(temp_popObj, opt.PF);
%                     opt.gd = GD(temp_popObj, opt.PF);
                    opt = compute_performance_metric(opt, opt.archiveObj, opt.archiveCV);
                    %clear opt;
                    disp(['Problem: ',opt.objfunction,', IGD:' num2str(opt.igd), ', GD: ', num2str(opt.gd), ', HV: ', num2str(opt.hv)]);
                end
            end
            

        end
    end
    disp('============END OF ALL RUNS==============');
end

%------------------------------END OF -FILE--------------------------------


