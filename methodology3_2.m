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

function opt = methodology3_2(opt)%generative, one by one
    
    repeated_times = 1;
    sign  = -1;
    l = 21;
    u = 1;%opt.numdir;
    
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
            
            opt.normalizedObj = normalize(opt, opt.archiveObj, opt.min_val, opt.max_val );
            [opt.archiveASF, opt.archiveCluster] = evaluateASF(opt);%evaluate ASF for all solutions, previous and new


            %----------------FIND ACTIVE SET-----------------------------------
            opt = find_active_set(opt);

            %--------------UNIQUES POPULATION----------------------------------
            opt = unique_population(opt);

            %--------------PLOT DIRECTION--------------------------------------
            
            
            
            %plot([opt.min_val(1)-w(1)*opt.min_val(1) w_org(1)], [opt.min_val(2)-w(2)*opt.min_val(1) w_org(2)],'-r','LineWidth',2);
            %plot_points(opt);
            %%{
            %if mod(opt.funcEval,100)==0 %&& opt.funcEval < 600
                %plot_points(opt);
                %hold on;
                %w = opt.curdir;
                %{
                for i=1:opt.numdir
                    w = opt.dirs(i,:);
                    w_org = zeros(1,2);
                    w_org(1) =  opt.min_val(1)+ w(1)*(opt.max_val(1)-opt.min_val(1));
                    w_org(2) =  opt.min_val(2)+ w(2)*(opt.max_val(2)-opt.min_val(2));
                    %w_org(1,:) = opt.min_val + (opt.max_val-opt.min_val).*[(w(1)-0.5) w(1)+1];
                    %w_org(2,:) = opt.min_val + (opt.max_val-opt.min_val).*[(w(2)-0.5) w(2)+1];
                    %plot([(w(1)-0.5) w(1)+1],[(w(2)-0.5) w(2)+1],'-r','LineWidth',2);
                    plot([opt.min_val(1) w_org(1)], [opt.min_val(2) w_org(2)],'-r','LineWidth',2);
                    %plot([opt.min_val(1)-w(1)*opt.min_val(1) w_org(1)], [opt.min_val(2)-w(2)*opt.min_val(1) w_org(2)],'-r','LineWidth',2);
                    %w = w./norm(w);
                    %plot([0 w(1)*2],[0 w(2)*2],'-r', 'LineWidth',2);
                end
                %}
                %opt.min_val = [-274 4];
                %opt.max_val = [-42 76];

                %plot([(w(1)-0.5) w(1)+1],[(w(2)-0.5) w(2)+1],'-r','LineWidth',2);            
                %plot(opt.activeArchiveObj(:,1), opt.activeArchiveObj(:,2), 'go', 'MarkerSize', 5,'MarkerFaceColor',[0.95 0.4 0.8]);
                %hold off;
            %end
            %}
                       
            %----------MODEL ASF FOR SPECIFIC DIRECTION------------------------
            
            [opt.dmodel_asf, ~] = dacefit(opt.activeArchive, opt.activeArchiveASF, @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model asf
                 
            if opt.C>0
                if opt.methodology==3
                    
                    [opt.dmodel_cons, ~] = dacefit(opt.activeArchive, opt.activeArchiveCons, @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model cv
                    
                    %{
                    for i=1:opt.C
                        %[opt.dmodel_cons{i}, ~] = dacefit(opt.archive, opt.archiveCons(:,i), @regpoly1, @corrgauss, opt.theta, opt.lob, opt.upb);%model cv
                        [opt.dmodel_cons{i}, ~] = dacefit(opt.activeArchive, opt.activeArchiveCons(:,i), @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model cv
                    end
                    %}
%                     [opt.dmodel_cons{1}, ~] = dacefit(opt.activeArchive, opt.activeArchiveCons(:,1), @regpoly1, @corrgauss, opt.theta, opt.lob, opt.upb);
%                     [opt.dmodel_cons{2}, ~] = dacefit(opt.activeArchive, opt.activeArchiveCons(:,2), @regpoly1, @corrgauss, opt.theta, opt.lob, opt.upb);
%                     [opt.dmodel_cons{3}, ~] = dacefit(opt.activeArchive, opt.activeArchiveCons(:,3), @regpoly1, @corrgauss, opt.theta, opt.lob, opt.upb);
%                     [opt.dmodel_cons{4}, ~] = dacefit(opt.activeArchive, opt.activeArchiveCons(:,4), @regpoly1, @corrgauss, opt.theta, opt.lob, opt.upb);
%                     [opt.dmodel_cons{5}, ~] = dacefit(opt.activeArchive, opt.activeArchiveCons(:,5), @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);
%                     [opt.dmodel_cons{6}, ~] = dacefit(opt.activeArchive, opt.activeArchiveCons(:,6), @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);
                    %}
                    %[opt.dmodel_cons{i}, ~] = dacefit(opt.archive, opt.archiveCons(:,i), @regpoly1, @corrgauss, opt.theta, opt.lob, opt.upb);%model cv
                    %[opt.dmodel_cons, ~] = dacefit(opt.activeArchive, opt.activeArchiveCons, @regpoly1, @corrgauss, opt.theta, opt.lob, opt.upb);%model cv
% 
%                     opt.regC{1} = quadratic_model(opt.activeArchive, opt.activeArchiveCons(:,1),1);%model constraints
%                     opt.regC{2} = quadratic_model(opt.activeArchive, opt.activeArchiveCons(:,2),1);%model constraints
%                     opt.regC{3} = quadratic_model(opt.activeArchive, opt.activeArchiveCons(:,3),1);%model constraints
%                     opt.regC{4} = quadratic_model(opt.activeArchive, opt.activeArchiveCons(:,4),1);%model constraints
%                     opt.regC{5} = quadratic_model(opt.activeArchive, opt.activeArchiveCons(:,5),2);%model constraints
%                     opt.regC{6} = quadratic_model(opt.activeArchive, opt.activeArchiveCons(:,6),2);%model constraints
% 
%                     opt.regO{1} = quadratic_model(opt.activeArchive, opt.activeArchiveObj(:,1),2);%model constraints
%                     opt.regO{2} = quadratic_model(opt.activeArchive, opt.activeArchiveObj(:,2),2);%model constraints

                elseif opt.methodology==4
                    %[opt.dmodel_cv, ~] = dacefit(opt.archive, opt.archiveCV, @regpoly1, @corrgauss, opt.theta, opt.lob, opt.upb);%model cv
                    [opt.dmodel_cv, ~] = dacefit(opt.activeArchive, opt.archiveCV, @regpoly1, @corrgauss, opt.theta, opt.lob, opt.upb);%model cv
                end
            end

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
        plot_points(opt);
    end
    
    %==========UPDATE HYPER-SPHERICAL TRUST REGION=========================
    if opt.trust_region_option_m3==1
        if opt.TrustDistObj*0.95 < 0.03
            opt.TrustDistObj = 0.03;
        else
            opt.TrustDistObj = opt.TrustDistObj*0.95;
        end

        if (0.95*opt.TrustDistVar)<0.03
            opt.TrustDistVar = 0.03;
        else
            opt.TrustDistVar =  0.95*opt.TrustDistVar;
        end
    elseif opt.trust_region_option_m3==2
        if opt.delta*0.95<0.005
            opt.delta = 0.005;
        else
            opt.delta = 0.95*opt.delta;
        end
    end
end


%------------------------------END OF -FILE--------------------------------