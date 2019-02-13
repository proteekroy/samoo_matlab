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

function opt = methodology5(opt)%generative, one by one
    

    %--------------UNIQUES POPULATION--------------------------------------
    [~,ia,~] = unique(opt.archive,'rows');
    temp_archive = opt.archive(ia,:);
            
    repeated_times = 1;
    sign  = 1;
    l = 1;
    u = opt.numdir;
    
    %for dir=1:opt.numdir%:-1:1
    
    %for dir = 1:1: size(opt.dirs,1) 
    for dir=l:sign:u   
        
        sign = sign*(-1);
        temp = u;
        u = l;
        l = temp;
        %------------ASF WITH CURRENT DIRECTION----------------------------
                    
        opt.curdir = opt.dirs(dir,:);%current direction
        opt.curcluster = dir;%current cluster number
        
        %opt.curdir = opt.dirs(21,:);%current direction
        %opt.curcluster = 21;%current cluster number
        
        for repeat = 1:repeated_times
            
            opt.normalizedObj = normalize(opt, opt.archiveObj, opt.min_val, opt.max_val);
            %opt.normalizedObj =  nsga3_normalization(opt, opt.archiveObj);
            [opt.archiveASF, opt.archiveCluster] = evaluateASF(opt);%evaluate ASF for all solutions, previous and new
            %opt.archiveASF = opt.archiveASFAll{dir};
            %opt.archiveCluster = dir*ones(1,size(opt.archive,1));

            %----------------FIND ACTIVE SET-----------------------------------
            %opt = find_active_set2(opt);
            %opt = unique_population(opt);
            %temp_normobj = opt.normalizedObj(ia,:);
            temp_ASF = opt.archiveASF(ia,:);
%             temp_cons = opt.archiveCons(ia, :);
%             temp_cvmodified = opt.archiveCVModified(ia, :);
%             temp_activeArchive = opt.activeArchive(ia,:);
%             temp_activeArchiveObj = opt.activeArchiveObj(ia,:);
%             temp_activeArchiveASF = opt.activeArchiveASF(ia,:);
%             temp_activeArchiveCons = opt.activeArchiveCons(ia,:);
%             temp_activeArchiveCV = opt.activeArchiveCV(ia,:);

            %----------MODEL ASF FOR SPECIFIC DIRECTION------------------------
            
            %[opt.dmodel_asf, ~] = dacefit(opt.activeArchive, opt.activeArchiveASF, @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model asf
            [opt.dmodel_asf, ~] = dacefit(temp_archive, temp_ASF, @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model asf
                 
            %--------------------RGA MODEL-------------------------------------

            pop = methodology3_rga_metamodel2(opt);%it will return one good solution which maximizes current reference direction


            %-----------------STORE RESULTS & PLOT-----------------------------

            [popObj, popCons] = evaluate_pop(opt, pop);%high fidelity computation
            opt = store_results(opt, pop, popObj, popCons); %HI-FI+ASF+KKTPM+CLUSTER
            
            %-------------FUNCTION EVALUATION CHECK----------------------------
            opt.funcEval = size(opt.archive,1);%number of function evaluations
            if opt.funcEval>=opt.totalFuncEval
                break;
            end
        end
        
        if opt.plotOption==1
            plot_points(opt);
        end
    end
    
end


%------------------------------END OF -FILE--------------------------------