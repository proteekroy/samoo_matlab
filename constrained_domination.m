function [d] = constrained_domination(obj1, obj2, cv1, cv2)
    
    if (cv1<=0 && cv2<=0) || (cv1-cv2)<1e-16 %both are feasible or same CV
        dom = lex_dominate(obj1, obj2);
        
        if dom == 1 %p1 dominates p2
            d = 1;%1 dominates 2
        elseif dom == 3 % p2 dominates p1
            d = 2;%2 dominates 1
        else
            d = 3;%both nondominated
        end
    else
        if cv1 < cv2%p1 less constraint violation
            d=1;%1 dominates 2
        else
            d=2;%2 dominates 1 
        end
    end
end