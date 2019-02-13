function [R,F] = bos(objective)
    
    if isempty(objective)
       R = [];
       F = [];
       return;
    elseif size(objective,1)==1
        R = 1;
        F(1).f = 1;
        return;
    end
    
    %objective = load('cloud_10000_20_1.txt');
    %--------INITIALIZATION---------------
    [n,m]= size(objective);
    R = ones(n,1);
    
    if(m<4)
        m1 = m;
    else
        m1 = floor(min(ceil(log2(n)),m));
    end
    %---------SORTING PART-----------------
    Q = zeros(n, m1);
    lex_order = zeros(n,1);
    
    %--------FIND LEX ORDER----------------
    [~,Q(:,1)] = sortrows(objective);
    
    for i = 1:n 
        lex_order(Q(i,1)) = i;
    end
    
    for i = 2:m1
       %H = horzcat(objective(:,i), lex_order);%in case of tie use lex order
       H = horzcat(objective(:,i), objective(:,1));
       [~,Q(:,i)] = sortrows(H);
    end
    
   

    
    %-------RANKING PART-------------------
    done = zeros(n, 1);
    total = 0;
    totalfront = 1;
    L = cell(m1, n);
    for i = 1:n
        for j = 1:m1 
            s  = Q(i,j);
            if done(s) == 1
                L{j, R(s)} = horzcat(s, L{j, R(s)});
                continue;
            end
            total = total + 1;
            done(s) = 1;
            
            
            for k = 1: totalfront %for all front
                d = 0;
                sz = size(L{j,k},2);
                for l = 1:sz % for all elements
                    d = lex_dominate(objective(L{j,k}(l),:), objective(s,:));
                    if d == 1
                       break; 
                    end
                end
                if d == 0 %not dominates
                    R(s) = k;
                    L{j, k} = horzcat(s, L{j, k});    
                    break;
                elseif d==1 && k==totalfront
                    totalfront = totalfront + 1;
                    R(s) = totalfront;
                    L{j,totalfront} = horzcat(s, L{j,totalfront});
                    break;
                
                end
            end
            
            
            if total==n
               break; 
            end
        end
    end
    
    F(n).f = [];
        
    for i=1:n
        F(R(i)).f = [F(R(i)).f i];
    end
    
    
    
    
    %{
    count = zeros(n,1);
    for i=1:n
        count(R(i)) = count(R(i))+1;
    end
    
    
    for i = 1:totalfront
        str = strcat('#Front ',num2str(i),' = ',num2str(count(i)));
        disp(str); 
    end
    %}
end



function [d] = lex_dominate(obj1, obj2)
    
    
    equal = 1;
    d = 1;
    sz = size(obj1,2);
    for i = 1:sz	
        if obj1(i) > obj2(i)
            d = 0;
            break;
        elseif(equal==1 && obj1(i) < obj2(i))
			equal = 0;
        end
        
    end
    if d ==1 && equal==1 %check if both solutions are equal
        d = 0;
    end
end