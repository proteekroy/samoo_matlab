function opt = find_leader(opt)


    %insert Non-dominated points as a leader
    cls = opt.archiveCluster(opt.ParetoIndex,:);
    opt.LeaderIndex = cell(1, opt.numdir);
    
    %======================================================================
    %first put Non-dominate front into cluster, if multiple points, choose
    %the one with minimum ASF w.r.t that direction
    %======================================================================
    
    for i=1:opt.numdir
        temp = opt.ParetoIndex(cls == i);%all solutions of class "i"
        if~isempty(temp)
            asfval = opt.archiveASF(temp);
            [~,index] = min(asfval);
            opt.LeaderIndex{i} = temp(index);
        else
           opt.LeaderIndex{i} = [];
        end
    end

    %======================================================================
    %find leaders for empty directions, ASF should not matter
    %======================================================================
    
    
    D = zeros(size(opt.archiveObj,1),opt.numdir);%matrix of orthogonal distances
     
    for j=1:opt.numdir
        w = opt.dirs(j,:);
        w_org = zeros(1,opt.M);
        for i=1:opt.M
            w_org(i) =  opt.min_val(i)+ w(i)*(opt.max_val(i)-opt.min_val(i));
        end
        
        for i=1:size(opt.archiveObj,1)
            p = opt.archiveObj(i,:);
            a = opt.min_val;
            b = w_org;
            ab = b-a;
            ap = p-a;
            p_cos_theta = dot(ab,ap)/norm(ab);
            D(i,j) = sqrt(norm(ap)^2- p_cos_theta^2);
       
        end
    end
       
    opt.D = D;
    
    for i=1:opt.numdir    
        
        if isempty(opt.LeaderIndex{i})%if some cluster is empty
            trust_index = find(opt.D(:,i)<=opt.TrustDistObj);
            if ~isempty(trust_index)%any solution within trust region
                asfval = opt.archiveASF(trust_index);
                [~,index] = min(asfval);
                opt.LeaderIndex{i} = trust_index(index);
            else
                [~, I] = min(opt.D(:,i));
                opt.LeaderIndex{i} = I;%closest to the line
            end
        end
    end

end