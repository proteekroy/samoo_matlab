function opt = find_active_set2(opt)


    %------------FIND ACTIVE ARCHIVE---------------------------------------
    switch(opt.methodology)
            
        case {11, 21, 31, 41, 5, 6}%when we find active set according to current direction
            
            d = zeros(size(opt.archiveObj,1),1);
            w = opt.curdir;
            %%{
            w_org = zeros(1,opt.M);
            for j=1:opt.M
                w_org(j) =  opt.min_val(j)+ w(j)*(opt.max_val(j)-opt.min_val(j));
            end
            for i=1:size(opt.archiveObj,1)
                
                p = opt.archiveObj(i,:);
                a = opt.ref_point{opt.curcluster};%opt.min_val;
                b = w_org;
                ab = b-a;
                ap = p-a;
                p_cos_theta = dot(ab,ap)/norm(ab);
                d(i) = sqrt(norm(ap)^2- p_cos_theta^2);
                
                %{
                obj = opt.archiveObj(i,:);
                if(opt.M==2)
                    pt = [obj 0];
                    v1 = [opt.min_val 0];
                    v2 = [w_org 0];
                else
                    pt = obj;
                    v1 = opt.min_val;
                    v2 = w_org;
                end
                %}
                %d(i) = point_to_line(pt, v1, v2);
                %d(i) =
                %norm(obj-((w_org*obj')*w_org)/(norm(w_org)*norm(w_org)));%This
                %is accurate when  min value is zero
            end
            
            %}
            %{
            w_org = w;
            for i=1:size(opt.archiveObj,1)
                obj = opt.normalizedObj(i,:);
                d(i) = norm(obj-((w_org*obj')*w_org)/(norm(w_org)*norm(w_org)));
            end
            %}
            
            index = find(d<=opt.TrustDistObj);
            %{
            opt.activeArchiveObj = opt.archiveObj(index,:);
            plot_points(opt);
            hold all;
            plot([opt.min_val(1) w_org(1)], [opt.min_val(2) w_org(2)],'-r','LineWidth',2);
            plot(opt.archiveObj(index,1), opt.archiveObj(index,2), 'co', 'MarkerEdgeColor','c','MarkerFaceColor','c','MarkerSize', 5);
            %}
            if(size(index,1)>opt.activeSetSize)%many solutions are within trust region
                %I = randi(size(index,1), [1 opt.activeSetSize]);%choose any random set (may be the set which improves condition number) 
                %%{
                leader = opt.archive(opt.LeaderIndex{opt.curcluster}, :);
                leader = normalize(opt, leader, opt.bound(1,:), opt.bound(2,:));
                norm_pop = normalize(opt, opt.archive(index,:), opt.bound(1,:), opt.bound(2,:)); 
                div = pdist2(norm_pop, leader);
                [~, index2] = sort(div,'ascend');
                I = index(index2(1:opt.activeSetSize));
                %}
                
            else % less than required number of points are within trust region
                [~,index] = sort(d, 'ascend');
                I = index(1:opt.activeSetSize);
            end
            %plot(opt.archiveObj(I,1), opt.archiveObj(I,2), 'co', 'MarkerEdgeColor','c','MarkerFaceColor','c','MarkerSize', 5);
            
            
            %{
            if opt.C>0
                I1 = find(opt.archiveCV<=0);
                %if size(I1,1)==1%reshape
                %   I1 =  I1';
                %end
                I3 = find(opt.archiveCV>0);
                %if size(I3,1)==1
                %   I3 =  I3';
                %end
            else
                I1 = (1:size(opt.archiveASF,1))';
                I3 = [];
            end
            
            if ~isempty(I1)
                [~,I2] = sort(opt.archiveASF(I1));
                I1 = I1(I2);
            end
            
            if ~isempty(I3)
                [~,I4] = sort(opt.archiveASF(I3));
                I3 = I3(I4);
            end
            
            I = vertcat(I1,I3);   
            p = min(opt.activeSetSize, size(I,1));
            I = I(1:p);
            %}
            opt.activeArchive = opt.archive(I,:);
            opt.activeArchiveObj = opt.archiveObj(I,:);
            opt.activeArchiveASF = opt.archiveASF(I,:);
            %if opt.C>0
            opt.activeArchiveCons = opt.archiveCons(I,:);
            opt.activeArchiveCV = opt.archiveCV(I,:);
            %opt.activeArchiveKKTPM = opt.archiveKKTPM(I,:);
            %end
            %}
%         case 6%Model KKTPM        
%             %[opt.activeArchive, opt.activeArchiveObj, opt.activeArchiveCV, opt.activeArchiveKKTPM, opt.activeArchiveCons] = select_best_individual(opt, opt.archive, opt.archiveObj, opt.archiveCV, opt.archiveKKTPM,  opt.archiveCons, opt.activeSetSize);
%             %index = select_best_individual(opt, opt.archive, opt.archiveObj, opt.archiveCV, opt.archiveASF, opt.archiveCons, opt.activeSetSize);  
%             index = select_best_individual(opt, opt.InitialActiveSetSize);
%             opt.activeArchive = opt.archive(index,:);
%             opt.activeArchiveObj = opt.archiveObj(index,:);
%             opt.activeArchiveCV = opt.archiveCV(index,:);
%             opt.activeArchiveCons = opt.archiveCons(index,:);
%             opt.archiveKKTPM = opt.archiveKKTPM(index,:);
%         case 7 %MODEL MEMO        
%             %[opt.activeArchive, opt.activeArchiveObj, opt.activeArchiveCV, opt.activeArchiveASF, opt.activeArchiveCons] = 
%             index = select_best_individual(opt, opt.InitialActiveSetSize);%opt.archive, opt.archiveObj, opt.archiveCV, opt.archiveASF, opt.archiveCons, opt.activeSetSize);  
%             opt.activeArchive = opt.archive(index,:);
%             opt.activeArchiveObj = opt.archiveObj(index,:);
%             opt.activeArchiveCV = opt.archiveCV(index,:);
%             opt.activeArchiveCons = opt.archiveCons(index,:);
%             opt.activeArchiveASF = opt.archiveASF(index,:);

        case 8 %MODEL Trust Region
            
            if opt.funcEval<opt.switch_method
                %------------Find Div value for all solution-------------------
                nd_set = opt.archive(opt.ParetoIndex,:);
                nd_set = normalize(opt, nd_set, opt.bound(1,:), opt.bound(2,:));
                norm_pop = normalize(opt, opt.archive, opt.bound(1,:), opt.bound(2,:)); 
                opt.div = trust_distance(opt, norm_pop, nd_set);

                %-----------Find Leaders for each direction--------------------
                opt.ExplorationIndex = find_diverse_points(opt);
                opt.OutreachSolutions = outreach_inside(opt.archive(opt.ExplorationIndex,:), 1);%find one solution for each non-dominated solution

                %-----------Choose Active Set for metamodeling-----------------
                I1 = choose_active_set_TrustRegion(opt);
                I2 = find(opt.div<opt.delta);%all the elements which has
                I = horzcat(I1, I2');
                I = unique(I);

                opt.activeArchive = opt.archive(I,:);
                opt.activeArchiveObj = opt.archiveObj(I,:);
                opt.activeArchiveCV = opt.archiveCV(I,:);
                opt.activeArchiveASF = opt.archiveASF(I,:);
                opt.activeArchiveCons = opt.archiveCons(I,:);
            else
                if opt.C>0
                    I1 = find(opt.archiveCV<=0);
                    I3 = find(opt.archiveCV>0);
                else
                    I1 = (1:size(opt.archiveASF,1))';
                    I3 = [];
                end

                if ~isempty(I1)
                    [~,I2] = sort(opt.archiveASF(I1));
                    I1 = I1(I2);
                end

                if ~isempty(I3)
                    [~,I4] = sort(opt.archiveASF(I3));
                    I3 = I3(I4);
                end
                
                if ~isempty(I1)%there exists feasible solutions
                    feasible_pop_cluster = opt.archiveCluster(I1);
                    index = find(feasible_pop_cluster == opt.curcluster, 1); 
                    if ~isempty(index)%at least one solution is in current cluster
                        index_best = I1(1);
                        I1 = knnsearch(opt.archive, opt.archive(index_best,:), 'K', opt.activeSetSize)';%I1 is now changed
                        %{
                        plot_points(opt);
                        hold on;
                        w = opt.curdir;
                        plot([(w(1)-0.5) w(1)+1],[(w(2)-0.5) w(2)+1],'-r','LineWidth',2);
                        plot(opt.archive(index_best, 1), opt.archive(index_best, 2),'ko','MarkerEdgeColor','y','MarkerFaceColor','y','MarkerSize', 10);
                        hold off;
                        %}
                    end
                end

                I = vertcat(I1,I3);   
                p = min(opt.activeSetSize, size(I,1));
                I = I(1:p);

                opt.activeArchive = opt.archive(I,:);
                opt.activeArchiveObj = opt.archiveObj(I,:);
                opt.activeArchiveASF = opt.archiveASF(I,:);
                opt.activeArchiveCons = opt.archiveCons(I,:);
                opt.activeArchiveCV = opt.archiveCV(I,:);
            end
                       
            %----plot nearest points---------------------------------------
            %if opt.funcEval>opt.changeDelta
            %    plot_points(opt);
            %end
        otherwise
            disp('Methodology Out Of Bound Exception');
    end
end

function d = point_to_line(pt, v1, v2)

      a = v1 - v2;
      b = pt - v2;
      d = norm(cross(a,b)) / norm(a);
      
end


function d = perpendicular_distance(p,a, b)
    ab = b-a;
    ap = p-a;
    p_cos_theta = dot(ab,ap)/norm(ab);
    d = sqrt(norm(ap)^2- p_cos_theta^2);
end


