    

function plot_points(opt)

    %h = figure('Visible','Off');
    h = opt.fig;
    box on;
    grid on;
    %clf(h)
    hold all;
    currentPareto = opt.archiveObj(opt.ParetoIndex,:);
    %currentLeader = opt.archiveObj(opt.LeaderIndex{opt.curcluster},:);
    if opt.M==3
          plot3(opt.archiveObj(:,1), opt.archiveObj(:,2), opt.archiveObj(:,3), 'bo', 'MarkerSize', 5);
          plot3(currentPareto(:,1), currentPareto(:,2),currentPareto(:,3), 'ro','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize', 5);
%         for i=1:opt.numdir
%             plot3(opt.archiveObj(opt.LeaderIndex{i},1), opt.archiveObj(opt.LeaderIndex{i},2),opt.archiveObj(opt.LeaderIndex{i},3),'mo','MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize', 5);
%         end
        plot3(opt.PF(:,1), opt.PF(:,2), opt.PF(:,3), 'ko', 'MarkerFaceColor','k', 'MarkerSize', 5);
        %plot3(opt.activeArchiveObj(:,1), opt.activeArchiveObj(:,2), opt.activeArchiveObj(:,3), 'go', 'MarkerSize', 5);
    elseif opt.M==2
        plot(opt.archiveObj(:,1), opt.archiveObj(:,2), 'bo', 'MarkerSize', 2);
        plot(currentPareto(:,1), currentPareto(:,2),'ro','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize', 5);
%         for i=1:opt.numdir
%             plot(opt.archiveObj(opt.LeaderIndex{i},1), opt.archiveObj(opt.LeaderIndex{i},2),'mo','MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize', 3);
%         end
        plot(opt.PF(:,1), opt.PF(:,2), 'ko', 'MarkerFaceColor','k', 'MarkerSize', 3);
        %plot(opt.activeArchiveObj(:,1), opt.activeArchiveObj(:,2), 'go', 'MarkerSize', 5);
    elseif opt.M>3        
        
        obj_all = vertcat(opt.archiveObj, opt.paretofront);
        I = true(size(obj_all,1), 1);
        I(1:size(opt.archive, 1))=false;
        species = cell(size(obj_all,1),1);
        species(~I) = {'Obtained Solutions'};
        species(I) = {'Non-dominated Solutions'};
        label = cell(1,opt.M);
%         label{1} = 'f_1';
%         label{2} = 'f_2';
        for i=1:opt.M
            label{i} = strcat('f_',num2str(i));
        end
        %species = fliplr(species);
        %P = fliplr(P);
        parallelcoords(obj_all,'group', species, 'labels',label,'LineWidth',2);
        ylabel('Objective Value');
        
    end
    
    %xlim([opt.min_val(1) opt.max_val(1)])
    %ylim([opt.min_val(2) opt.max_val(2)])
    %{
    plot(opt.normalizedObj(:,1), opt.normalizedObj(:,2), 'bo', 'MarkerSize', 5);
    
    for i = 1:size(opt.dirs,1)
        w = opt.dirs(i,:);
        w = w./norm(w);
        plot([0 w(1)*2],[0 w(2)*2],'-k');
        plot([(w(1)-0.5) w(1)+1],[(w(2)-0.5) w(2)+1],'-k');
    end
    %}
    %{
    leader = opt.archive(opt.LeaderIndex,:);
    norm_pop = normalize(opt, opt.archive, opt.bound(1,:), opt.bound(2,:)); 
    div = trust_distance(opt, norm_pop, leader);
    div = div/opt.V;%normalize by the number of variables
    index = div<opt.delta;
    
    plot(opt.archiveObj(index,1), opt.archiveObj(index,2), 'co', 'MarkerSize', 5);
    %}
    drawnow;
    %{
    name = strcat('video/methodology1_1_',num2str(opt.funcEval),'.jpg');
    saveas(h,name);
    name = strcat('video/methodology1_1_',num2str(opt.funcEval),'.fig');
    saveas(h,name);
    close(h);
    %}
end