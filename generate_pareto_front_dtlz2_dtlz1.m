
function generate_pareto_front_dtlz2_dtlz1()

    close all;
    %======================Layer Wise======================================
    %{
    M=5;
    w = layer_wise(M);
    scatter3(w(:,1), w(:,2), w(:,3),'ob');
    norm_w = sqrt(sum(abs(w).^2,2));
    w_unit = w./norm_w;
    
    if M==3
        scatter3(w_unit(:,1), w_unit(:,2), w_unit(:,3),'ob');
    else
        parallelcoords(w_unit);
    end
    dlmwrite(strcat('DTLZ2.',num2str(M),'D.pf'), w_unit,'delimiter','\t','precision',6);
    %}
    %===========================IGD========================================
    %%{
    M = 5;
    N = 210;%715;%126;
    w = initweight(M, N)';
    norm_w = sqrt(sum(abs(w).^2,2));
    norm_w = repmat(norm_w, 1, M);
    w_unit = w./norm_w;
    disp(sum(w_unit.^2, 2));
    figure;
    if M==3
        scatter3(w_unit(:,1), w_unit(:,2), w_unit(:,3),'ob');
    else
        parallelcoords(w_unit);
    end
    g = zeros(N, 1);
    for i=1:size(w_unit,1)
        g(i) = c2dtlz2_cons(w_unit(i,:));
    end
    index = g<=0;
    w_unit_c2dtlz2 = w_unit(index, :);
    
    figure;
    if M==2
        scatter(w_unit_c2dtlz2(:,1), w_unit_c2dtlz2(:,2),'ob');
    elseif M==3
        scatter3(w_unit_c2dtlz2(:,1), w_unit_c2dtlz2(:,2), w_unit_c2dtlz2(:,3),'ob');
    else
        parallelcoords(w_unit_c2dtlz2);
    end
    w_dtlz1 = w/2;
    disp(sum(w_dtlz1, 2));
    figure;
    if M==2
        scatter(w_dtlz1(:,1), w_dtlz1(:,2),'ob');
    elseif M==3
        scatter3(w_dtlz1(:,1), w_dtlz1(:,2), w_dtlz1(:,3),'ob');
    else
        parallelcoords(w_dtlz1);
    end
    
    
    dlmwrite(strcat('DTLZ2.',num2str(M),'D.pf'), w_unit,'delimiter','\t','precision',6);
    dlmwrite(strcat('DTLZ4.',num2str(M),'D.pf'), w_unit,'delimiter','\t','precision',6);
    dlmwrite(strcat('DTLZ3.',num2str(M),'D.pf'), w_unit,'delimiter','\t','precision',6);
    dlmwrite(strcat('C2DTLZ2.',num2str(M),'D.pf'), w_unit_c2dtlz2,'delimiter','\t','precision',6);
    dlmwrite(strcat('DTLZ1.',num2str(M),'D.pf'), w_dtlz1,'delimiter','\t','precision',6);   
    
    %===========================GD=========================================
    N = 715;%500;%715;%126;
    w = initweight(M, N)';
    norm_w = sqrt(sum(abs(w).^2,2));
    norm_w = repmat(norm_w, 1, M);
    w_unit = w./norm_w;
    disp(sum(w_unit.^2, 2));
    figure;
    
    if M==2
        scatter(w_unit(:,1), w_unit(:,2),'ob');
    elseif M==3
        scatter3(w_unit(:,1), w_unit(:,2), w_unit(:,3),'ob');
    else
        parallelcoords(w_unit);
    end
    g = zeros(N, 1);
    for i=1:size(w_unit,1)
        g(i) = c2dtlz2_cons(w_unit(i,:));
    end
    index = g<=0;
    w_unit_c2dtlz2 = w_unit(index, :);
    
    figure;
    
    if M==2
        scatter(w_unit_c2dtlz2(:,1), w_unit_c2dtlz2(:,2),'ob');
    elseif M==3
        scatter3(w_unit_c2dtlz2(:,1), w_unit_c2dtlz2(:,2), w_unit_c2dtlz2(:,3),'ob');
    else
        parallelcoords(w_unit_c2dtlz2);
    end
    w_dtlz1 = w/2;
    disp(sum(w_dtlz1, 2));
    figure;
    if M==2
        scatter(w_dtlz1(:,1), w_dtlz1(:,2),'ob');
    elseif M==3
        scatter3(w_dtlz1(:,1), w_dtlz1(:,2), w_dtlz1(:,3),'ob');
    else
        parallelcoords(w_dtlz1);
    end
    dlmwrite(strcat('GD/','DTLZ2.',num2str(M),'D.pf'), w_unit,'delimiter','\t','precision',6);
    dlmwrite(strcat('GD/','DTLZ4.',num2str(M),'D.pf'), w_unit,'delimiter','\t','precision',6);
    dlmwrite(strcat('GD/','DTLZ3.',num2str(M),'D.pf'), w_unit,'delimiter','\t','precision',6);
    dlmwrite(strcat('GD/','C2DTLZ2.',num2str(M),'D.pf'), w_unit_c2dtlz2,'delimiter','\t','precision',6);
    dlmwrite(strcat('GD/','DTLZ1.',num2str(M),'D.pf'), w_dtlz1,'delimiter','\t','precision',6);
    
    %}
    
    %{
    N = 40;%1000;%45;%18;
    x1 = linspace(0, 1, N)';
    %x2 = sqrt(1-x1);
    %{
    [x1,x2] = meshgrid(x1);
    x = zeros(N*N,3);
    x(:,1) = x1(:);
    x(:,2) = x2(:);
    %}
    x2 = zeros(N,1);
    x(:,1) = x1(:);
    x(:,2) = x2(:);
    f = dtlz7(x, 2);
    f = unique(f,'rows');
    [R] = bos(f);
    f = f(R==1,:);
    %{
    figure;hold all;box on; grid on;
    plot3(f(:,1),f(:,2),f(:,3),'ok');
    n1 = size(f,1);
    n2 = size(f,1)-91;
    %(f-repmat(min(f), n1, 1))./((repmat(max(f),n1,1)+1e-16)-repmat(min(f),n1,1));  
    for i=1:n2
        f_normalized = f;
        [Idx, d] = knnsearch(f_normalized,f_normalized,'K',2, 'Distance', 'euclidean');
        Idx = Idx(:,2);
        %d = mean(d(:, 1:5),2);
        d = d(:, 2);
        [d, I2] = sort(d, 'ascend');
        f(Idx(I2(1)),:) = [];
        n1 = size(f,1);
    end
    %}
    figure;hold all;box on; grid on;
    plot(f(:,1),f(:,2),'ok');
    %plot3(f(:,1),f(:,2),f(:,3),'ok');
    dlmwrite(strcat('DTLZ7',num2str(M),'D.pf'), f,'delimiter','\t','precision',6);
    %dlmwrite(strcat('GD/DTLZ7.',num2str(M),'D.pf'), f,'delimiter','\t','precision',6);
    %}
    
end


function dir = layer_wise(M)
    %p = 15;
    H1 = initweight(M, nchoosek(M+1-1,1))';%5
    H1_temp = initweight(M, nchoosek(M+2-1,2))';%15
    H2 = layered_weight(.75, H1_temp);%15
    H3 = layered_weight(.50, H1_temp);%15
    H4 = layered_weight(.35, H1_temp);%15
    dir = vertcat(H1, H2, H3, H4);
end


function [f] = dtlz7(x, M)

    n=size(x,2);
    
    f = zeros(size(x,1), M);
    
    nfunc = M;
    k = n - nfunc + 1;

    g = sum(x(:,nfunc:n),2);
    g = 1 + 9 * g ./ k;


    for i = 1:nfunc-1
        f(:,i) = x(:, i);
    end
    
    h2 = zeros(size(x,1), 1);
    
    for j = 1:nfunc-1
       h2 = h2 +  x(:,j) ./ (1 + g) .* (1 + sin(3 * 3.141592654 *  x(:,j))); 
    end

    h2 = nfunc - h2;
    f(:,nfunc) = (1 + g) .* h2;

end




function g = c2dtlz2_cons(f)

    nfunc = size(f,2);
    if(nfunc>3)
        r = 0.5;
    else
        r = 0.4;
    end
    
    v1 = realmax;
    v2 = 0.0;

    for i = 1: nfunc 
        sum1 = power(f(i)-1.0, 2.0);
        for j = 1: nfunc 
            if i ~= j 
                sum1 = sum1 + power(f(j), 2.0);
            end
        end

        v1 = min(v1, sum1 - power(r, 2.0));
        v2 = v2 + power(f(i) -  (1.0 / sqrt(nfunc)), 2.0);
    end

    c = min(v1, v2 - power(r, 2.0));
            
    g = [];
    if(c<=0)
        g(1) = c;
    else
        g(1) = c;
    end

end