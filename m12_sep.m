



function m12_sep(opt)

    [~,ia,~] = unique(opt.archive,'rows');
    solution = opt.archive(ia,:);
    solutionObj =  opt.archiveObj(ia,:);
    solutionCV = opt.archiveCV(ia,:);
    solutionCons = opt.archiveCons(ia,:);
    
    n = size(solution, 1);
    
    for i=1:n-1
        for j=2:n
            %{
            if solutionCV(i)<0 && solutionCV(j)<0 %both feasible
                
            elseif solutionCV(i)<0 && solutionCV(j)<0%one feasible, one infeasible
               
            elseif solutionCV(i)>0 && solutionCV(j)>0%both infeasible
            
            end
            %}
            d1 = constrained_domination(solutionObj(i,:), solutionObj(j,:), solutionCV(i,:), solutionCV(j,:));
            d2 = constrained_domination(solutionObj(i,:), solutionObj(j,:), solutionCV(i,:), solutionCV(j,:));
        end
    end
end


function [d] = constrained_domination(obj1, obj2, cv1, cv2)
    
    if (cv1<=0 && cv2(p2)<=0) || (cv1-cv2)<1e-16 %both are feasible or same CV
        dom = lex_dominate(obj1, obj2);
        
        if dom == 1 %p1 dominates p2
            d = 1;%1 dominates 2
        elseif dom == 3 % p2 dominates p1
            d = 2;%2 dominates 1
        else
            d = 3;%both nondominated
        end
    else
        if cv1<cv2%p1 less constraint violation
            d=1;%1 dominates 2
        else
            d=2;%2 dominates 1 
        end
    end
   

end