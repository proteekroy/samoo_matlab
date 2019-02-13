function opt = methodology31(opt)%generative, one by one
    
    %--------------UNIQUES POPULATION------------------------------
    [~,ia,~] = unique(opt.archive,'rows');
    temp_archive = opt.archive(ia,:);
    %temp_normobj = opt.normalizedObj(ia,:);
    
    temp_cons = opt.archiveCons(ia, :);
    temp_cvmodified = opt.archiveCVModified(ia, :);
    
    if opt.C>0
        if opt.methodology==31
            %[opt.dmodel_cons{i}, ~] = dacefit(temp_activeArchive, temp_activeArchiveCons(:,i), @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model cv
            [opt.dmodel_cons, ~] = dacefit(temp_archive, temp_cons, @regpoly1, @corrgauss, opt.theta, opt.lob, opt.upb);%model cv
        elseif opt.methodology==41
            [opt.dmodel_cv, ~] = dacefit(temp_archive, temp_cvmodified, @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model cv
        end
    end
    
    repeated_times = 1;
    sign  = 1;
    l = 1;
    u = opt.numdir;
    opt.temp = [];
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
        
        for repeat = 1:repeated_times
            
            opt.normalizedObj = normalize(opt, opt.archiveObj, opt.min_val, opt.max_val );
            %opt.normalizedObj =  nsga3_normalization(opt, opt.archiveObj);
            [opt.archiveASF, opt.archiveCluster] = evaluateASF(opt);%evaluate ASF for all solutions, previous and new
            %opt.archiveASF = opt.archiveASFAll{dir};
            
            temp_ASF = opt.archiveASF(ia,:);
            
            %----------------FIND ACTIVE SET-------------------------------
            %opt = find_active_set2(opt);

            
%             temp_activeArchive = opt.activeArchive(ia,:);
%             temp_activeArchiveObj = opt.activeArchiveObj(ia,:);
%             temp_activeArchiveASF = opt.activeArchiveASF(ia,:);
%             temp_activeArchiveCons = opt.activeArchiveCons(ia,:);
%             temp_activeArchiveCV = opt.activeArchiveCV(ia,:);
            %--------------PLOT DIRECTION----------------------------------
            
            
            
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
            
            %[opt.dmodel_asf, ~] = dacefit(temp_activeArchive, temp_activeArchiveASF, @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model asf
            [opt.dmodel_asf, ~] = dacefit(temp_archive, temp_ASF, @regpoly1, @corrgauss, opt.theta, opt.lob, opt.upb);%model asf
            
            

            %--------------------RGA MODEL---------------------------------

            [pop, opt] = methodology3_rga_metamodel2(opt);%it will return one good solution which maximizes current reference direction

            %-----------------STORE RESULTS & PLOT-------------------------

            [popObj, popCons] = evaluate_pop(opt, pop);%high fidelity computation
            popASF = calculate_Asf(popObj, opt.curdir);
            [popCV, ~]= evaluateCV(popCons);
            if opt.trust_region_update_option==3
                %if opt.adaptive_trust_region_option==1
                %    [opt, tempTrustRadiusDeltaK] = adjust_trust_region(opt, opt.temp.pop, opt.temp.popASF, opt.temp.popCV, popASF, popCV, opt.archive, opt.archiveObj, opt.archiveCV);
                %else
                [opt, tempTrustRadiusDeltaK] = adjust_trust_region_asf(opt, opt.temp.pop, opt.temp.popASF, opt.temp.popCV, popASF, popCV, opt.archive, opt.archiveObj, opt.archiveCV);
                %end

                opt.TrustRadiusDeltaK = horzcat(opt.TrustRadiusDeltaK, tempTrustRadiusDeltaK);
            end
            
            opt = store_results(opt, pop, popObj, popCons); %HI-FI+ASF+KKTPM+CLUSTER
            
            %-------------FUNCTION EVALUATION CHECK------------------------
            opt.funcEval = size(opt.archive,1);%number of function evaluations
            if opt.funcEval>=opt.totalFuncEval
                break;
            end
        end
        
        if opt.plotOption==1
            plot_points(opt);
        end
    end
    
    %-----------------UPDATE TRUST RADIUS----------------------------------
    
end


%------------------------------END OF -FILE--------------------------------