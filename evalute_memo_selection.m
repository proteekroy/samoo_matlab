function [asf, cluster] = evalute_memo_selection(opt, obj, cv)

    
    %--------------FIND CLUSTER--------------------------------------------
    S = zeros(size(obj,1),1);
    S(1:end) = intmax;

    cluster = ones(size(obj,1),1);
    minASFcluster = zeros(opt.numdir,1);
    minASFcluster(1:end) = intmax;

    for i=1:size(opt.dirs,1)
        w = opt.dirs(i,:);
        %w = w./norm(w);
        w_org = w;%zeros(1,2);
        %w_org(1) =  opt.min_val(1)+ w(1)*(opt.max_val(1)-opt.min_val(1));
        %w_org(2) =  opt.min_val(2)+ w(2)*(opt.max_val(2)-opt.min_val(2));
        Z  = max((obj-repmat(opt.min_val,size(obj,1),1))./repmat(w_org,size(obj,1),1),[],2);
        for j = 1:size(obj,1)
            if Z(j)<S(j)
                cluster(j) = i;
                S(j) = Z(j);
            end
            if Z(j)<minASFcluster(i)
                minASFcluster(i) = Z(j); 
            end
        end
    end

    asf =  S;

    %-------------COMBINE CV WITH ASF--------------------------------------
    if opt.C>0 

        index = find(cv<=0);
        feasibleASF = asf(index,:);
        fmax = max(feasibleASF);
        if isempty(fmax)
            fmax = max(asf);
        end
        %asf = asf./fmax;
        index = find(cv>0);
        if ~isempty(index)
            temp_asf = fmax + cv(index);
            asf(index) = temp_asf;
        end
    end


end