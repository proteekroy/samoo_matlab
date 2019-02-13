% Copyright [2017] [Proteek Chandan Roy]
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
% Proteek Chandan Roy, 2017
% Contact: royprote@egr.msu.edu, proteek_buet@yahoo.com


function opt = methodology6_memo(opt)%target multiple optimal solution

    %-------------------FIND ASF-------------------------------------------      
    
%     [opt.archiveASF, opt.archiveCluster] = evaluateASF(opt);%evaluate ASF for all solutions, previous and new
    
    %----------------FIND ACTIVE SET---------------------------------------
    %opt = find_active_set2(opt);
    
    %--------------UNIQUES POPULATION--------------------------------------
    %opt = unique_population(opt);        
%     [~,ia,~] = unique(opt.archive,'rows');
%     temp_archive = opt.archive(ia,:);
%     temp_archiveSF_VAL = opt.archiveSF_VAL(ia, :);
%                     
    %----------MODEL COMBINED (ASF+CV)-------------------------------------
    
    
    %[opt.dmodel_asf, ~] = dacefit(opt.archive, opt.archiveASF, @regpoly1, @corrgauss, opt.theta, opt.lob, opt.upb);%model asf               
    %[opt.dmodel_asf, ~] = dacefit(opt.activeArchive, opt.activeArchiveASF, @regpoly1, @corrgauss, opt.theta, opt.lob, opt.upb);%model asf               
    %[opt.dmodel_asf, ~] = dacefit(opt.archive, opt.archiveASF, @regpoly2, @corrcubic, opt.theta, opt.lob, opt.upb);%model asf       
    
    %if opt.funcEval <500
%        opt.net = train_neural_network(temp_archive, temp_archiveSF_VAL); 

%         opt.regC{1} = quadratic_model(opt.activeArchive, opt.activeArchiveCons(:,1),1);%model constraints
%         opt.regC{2} = quadratic_model(opt.activeArchive, opt.activeArchiveCons(:,2),1);%model constraints
%         opt.regC{3} = quadratic_model(opt.activeArchive, opt.activeArchiveCons(:,3),1);%model constraints
%         opt.regC{4} = quadratic_model(opt.activeArchive, opt.activeArchiveCons(:,4),1);%model constraints
%         opt.regC{5} = quadratic_model(opt.activeArchive, opt.activeArchiveCons(:,5),2);%model constraints
%         opt.regC{6} = quadratic_model(opt.activeArchive, opt.activeArchiveCons(:,6),2);%model constraints
% 
%         opt.regO{1} = quadratic_model(opt.activeArchive, opt.activeArchiveObj(:,1),2);%model objectives
%         opt.regO{2} = quadratic_model(opt.activeArchive, opt.activeArchiveObj(:,2),2);%model objectives
    %end
    opt = metamodel(opt, 1:size(opt.archive, 1));%build model
    pop = memo_niched_rga_metamodel(opt);%one or two good solution
                
                
    if(size(pop,1)+opt.funcEval>opt.totalFuncEval)%discard excess solution in last iteration
        pop = pop(1:opt.totalFuncEval-opt.funcEval,:);
    end
                              
    %-----------------STORE RESULTS & PLOT-----------------------------
    [popObj, popCons] = evaluate_pop(opt, pop);%high fidelity computation    
    opt = store_results(opt, pop, popObj, popCons); %HI-FI+ASF+KKTPM+CLUSTER
    
end


%------------------------------END OF -FILE--------------------------------

