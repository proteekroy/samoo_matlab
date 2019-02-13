function asf = calculate_Asf(obj, dir, utopian)

    %w = (dir-utopian)./norm(dir-utopian);%unit direction
    %w = dir;
    %asf  = max((obj-utopian)./w);
    
    %dir = dir./norm(dir);
    %asf = max(obj.*dir);
    %asf  = max(obj./dir); % max(w.*(obj-utopian));  
    
    %No Normalization is done
    %asf = max(obj - repmat(dir-1/size(obj,2),size(obj,1),1),[],2);%+0.01*sum(obj-(dir-0.5));%+0.0001*sum(obj.*dir);
    %asf = max(obj - repmat(dir + (1/size(obj,2)),size(obj,1),1),[],2);%+0.01*sum(obj-(dir-0.5));%+0.0001*sum(obj.*dir);
    
    asf = max( (obj-repmat(dir, size(obj,1),1)),[],2);
end