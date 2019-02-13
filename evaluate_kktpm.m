function [kktpm] = evaluate_kktpm(opt, pop)

    sz = size(pop, 1);
    kktpm = zeros(sz,1);
    
    for i=1:sz
        %disp(pop(i,:))
        kktpm(i) = KKT(opt, pop(i,:));
    end

end