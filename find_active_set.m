
%This function finds solutions sorted by their distance from reference
%lines
function activeSetIndex = find_active_set(opt, normalized_obj)

    n = size(normalized_obj,1);
    activeSetIndex = cell(1,size(opt.dirs, 1));
    
    for i=1:size(opt.dirs,1)
        w = opt.dirs(i,:);
        w = w./norm(w);
        obj_cos_theta = dot(normalized_obj,repmat(w,n,1),2);
        d = sqrt(vecnorm(normalized_obj,2,2).^2- obj_cos_theta.^2);
        
        [~,activeSetIndex{i}] = sort(d, 'ascend');
    end
    
end

