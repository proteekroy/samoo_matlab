function [cluster] = find_euclidean_niche(opt, popObj)

    
    d = zeros(size(popObj, 1), opt.numdir);  
    
    for j=1:opt.numdir
        w = opt.dirs(j,:);
        w_org = zeros(1, opt.M);
        for k=1:opt.M % translation in original objective space
            w_org(k) =  opt.min_val(k)+ w(k)*(opt.max_val(k)-opt.min_val(k));
        end
        
        for i=1:size(popObj,1)
            p = popObj(i,:);
            a = opt.ref_point{j};%use previous reference point, don't update
            b = w_org;
            ab = b - a;
            ap = p - a;
            p_cos_theta = dot(ab,ap)/norm(ab);
            d(i, j) = sqrt(norm(ap)^2- p_cos_theta^2);
        end
    end
    
    [~,cluster] = min(d,[],2);%cluster number for new points

end