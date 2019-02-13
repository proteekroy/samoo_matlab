function [cluster] = find_euclidean_cluster(opt, popObj)

    
    d = repmat(intmax, size(popObj, 1),1);%zeros(size(opt.archiveObj,1),1);
    cluster = -1*ones(size(popObj,1), 1);
    for k=1:opt.numdir
        w = opt.dirs(k,:);

        w_org = zeros(1,opt.M);
        for j=1:opt.M
            w_org(j) =  opt.min_val(j)+ w(j)*(opt.max_val(j)-opt.min_val(j));
        end
                
        for i=1:size(popObj,1)
            p = popObj(i,:);
            a = opt.ref_point{k};
            b = w_org;
            ab = b-a;
            ap = p-a;
            p_cos_theta = dot(ab,ap)/norm(ab);
            temp_d = sqrt(norm(ap)^2- p_cos_theta^2);
            if temp_d < d(i)
                cluster(i) = k;
                d(i) = temp_d;
            end
        end
    end
    
end