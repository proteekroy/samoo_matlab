function [asf, cluster] = minimumASF(opt, popObj)

    
    
    

    S = intmax;
    cluster = 0;
    
    for i=1:size(opt.dirs,1)
        w = opt.dirs(i,:);
        w = w./norm(w);
                
               
        %Z = max((popObj-opt.utopian)./w);
        Z = max(popObj-(w-0.5));
        if(Z<S)
            cluster = i;
            S = Z;
        end
    end
    
    asf =  S;
    
    
end

