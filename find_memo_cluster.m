
function [S, cluster] = find_memo_cluster(opt, popObj)

    %popObj = normalize(opt, popObj,  opt.min_val, opt.max_val);
    popObj =  nsga3_normalization(opt, popObj);
    %--------------FIND CLUSTER--------------------------------------------
    S = zeros(size(popObj,1),1);
    S(1:end) = intmax;
    
    cluster = -1*ones(size(popObj,1),1);
%     figure;
%     hold all;
%     for i=1:opt.numdir
%         w = opt.dirs(i,:);
%         w = w + 1e-10;
%         w = w./norm(w);
%         plot([opt.min_val(1) w(1)], [opt.min_val(2) w(2)],'-r','LineWidth',2);
%     end
%     hold off;
    
    for i=1:size(opt.dirs,1)
        w = opt.dirs(i,:);
        w = w + 1e-10;
        w = w./norm(w);
        
%         w_org = zeros(1,opt.M);
%         for j=1:opt.M
%             w_org(j) =  opt.min_val(j)+ w(j)*(opt.max_val(j)-opt.min_val(j));
%         end


        
        %Z = max(      (popObj-repmat(opt.ref_point{i}, size(popObj,1),1)),    [],2);
        Z = max((popObj./repmat(w, size(popObj,1), 1)),  [], 2);
        %asf  = max((obj-utopian)./w);
        %dir = dir./norm(dir);
            
        for j = 1:size(popObj, 1)
            if Z(j)<S(j)
                cluster(j) = i;
                S(j) = Z(j);
            end
        end
    end
    
    %plot_cluster(opt, popObj, cluster);
    
end