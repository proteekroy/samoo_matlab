function show_selection_function_m6()

    
    
    %-------necessary parameters-------------------------------------------
    opt.objfunction = 'tnk';
    opt.C = 2;
    opt.modality = 3;
    opt.rho = 0;
    opt.dirs = [0.01 0.99; 0.1 0.9; 0.2 0.8; 0.3 0.7; 0.4 0.6; 0.5 0.5; 0.6 0.4;0.7, 0.3; 0.8 0.2; 0.9 0.1; 0.99 0.01];%%w=11
    opt.numdir = size(opt.dirs,1);
    opt.methodology = 7;
    objective_data = load(strcat('IGD Calculation/methodology7_',opt.objfunction,'_obj_1.txt'));
    cv_data = load(strcat('IGD Calculation/methodology7_',opt.objfunction,'_cv_1.txt'));
    xvar = load(strcat('IGD Calculation/methodology7_',opt.objfunction,'_var_1.txt'));
    MEMO_data = zeros(size(objective_data,1),1);
    
    
    %-----------------bound------------------------------------------------
    opt.bound = zeros(2,2);
    opt.bound(2,:) = ones(1,2);
    if(strcmpi(opt.objfunction,'zdt4'))
        opt.bound(1,2:end) = opt.bound(1,2:end)+(-5);
        opt.bound(2,2:end) = opt.bound(2,2:end)*5;
    elseif(strcmpi(opt.objfunction,'UF1') ||strcmpi(opt.objfunction,'UF2') || strcmpi(opt.objfunction,'UF5') || strcmpi(opt.objfunction,'UF6') || strcmpi(opt.objfunction,'UF7'))
        opt.bound(1,2:end) = opt.bound(1,2:end)+(-1);
        opt.bound(2,2:end) = opt.bound(2,2:end)*1;
    elseif(strcmpi(opt.objfunction,'UF4') )
        opt.bound(1,2:end) = opt.bound(1,2:end)+(-2);
        opt.bound(2,2:end) = opt.bound(2,2:end)*2;
    elseif(strcmpi(opt.objfunction,'UF8') || strcmpi(opt.objfunction,'UF9')|| strcmpi(opt.objfunction,'UF10'))
        opt.bound(1,3:end) = opt.bound(1,3:end)+(-2);
        opt.bound(2,3:end) = opt.bound(2,3:end)*2;
    elseif(strcmpi(opt.objfunction,'BNH'))
        opt.bound(2,1)=5;
        opt.bound(2,2)=3;
    elseif(strcmpi(opt.objfunction,'OSY'))
        opt.bound(2,1)=10;%x1
        opt.bound(2,2)=10;
        opt.bound(2,6)=10;
        opt.bound(1,3)=1;
        opt.bound(1,5)=1;
        opt.bound(2,3)=5;
        opt.bound(2,5)=5;
        opt.bound(2,4)=6;
    elseif(strcmpi(opt.objfunction,'SRN'))
        opt.bound(1,1:end) = opt.bound(1,1:end)+(-20);
        opt.bound(2,1:end) = opt.bound(2,1:end)*20;
    elseif(strcmpi(opt.objfunction,'TNK'))
        opt.bound(2,1:end) = opt.bound(2,1:end)*pi;
    elseif(strcmpi(opt.objfunction,'WATER'))
        opt.bound(1,1:end) = opt.bound(1,1:end)+0.01;
        opt.bound(2,1) = 0.45;   
        opt.bound(2,2:end) = 0.10; 
    elseif(strcmpi(opt.objfunction,'carside'))
        %opt.bound(1,[1 3 4]) = 0.5; opt.bound(1,2) = 0.45; opt.bound(1,[6 7]) = 0.4; opt.bound(1,5) = 0.875;
        %opt.bound(2,[1 3 4]) = 1.5; opt.bound(2,2) = 1.35; opt.bound(2,[6 7]) = 1.2; opt.bound(2,5) = 2.625;  
        opt.bound(1,1:end) = [0.5 0.45 0.5 0.5 0.875 0.4 0.4];
        opt.bound(2,1:end) = [1.5 1.35 1.5 1.5 2.625 1.2 1.2];
    elseif(strcmpi(opt.objfunction,'welded'))
        opt.bound(1,1:end) = [0.125 0.1 0.1 0.125];
        opt.bound(2,1:end) = [5 10 10 5];
    elseif    (strcmpi(opt.objfunction,'wfg1') || strcmpi(opt.objfunction,'wfg2') || strcmpi(opt.objfunction,'wfg3') ...
            || strcmpi(opt.objfunction,'wfg4') || strcmpi(opt.objfunction,'wfg5') || strcmpi(opt.objfunction,'wfg6')...
            || strcmpi(opt.objfunction,'wfg7') || strcmpi(opt.objfunction,'wfg8') || strcmpi(opt.objfunction,'wfg9'))
        opt.bound(1,:) = zeros(1,opt.V);
        opt.bound(2,:) = ones(1,opt.V);
    end
    
    
    %------------------------Generate samples------------------------------
    
    N = 50;%50x50 grid
    z = [0 0];
    x1 = linspace(0,1,N);
    x2 = linspace(0,1,N);
    
    x1 = 1e-6+opt.bound(1,1)+x1*(opt.bound(2,1)-opt.bound(1,1));
    x2 = 1e-6+opt.bound(1,2)+x2*(opt.bound(2,2)-opt.bound(1,2));
    
    
    [X,Y] = meshgrid(x1,x2);
    sol = zeros(N*N, 2); 
    F = zeros(N*N, 2);  
    MEMO_longarray = zeros(N*N,1);
    
    if opt.C>0 
        G = zeros(N*N, opt.C);
    end
    F1 = zeros(size(X));
    F2 = zeros(size(X));
    MEMO = zeros(size(X));
    
    
    
    %----------------compute objective and ASF for grid--------------------
    l=1;
    k=1;
    for j=1:N*N
        sol(j,:) = [X(l,k) Y(l,k)];
        [F(j,:), u] =  high_fidelity_evaluation(opt, sol(j,:));
        if opt.C>0
           G(j,:) = u; 
        end
        F1(l,k) = F(j,1);
        F2(l,k) = F(j,2);
        arry = [];
        for p=1:size(opt.dirs,1)
            w = opt.dirs(p,:);
            arry(p) = max((F(j,1)-z(1))/w(1),( F(j,2)-z(2))/w(2));% + rho *(((X(j,k)-z(1))/w(1))+((Y(j,k)-z(2))/w(2)));%one ideal
        end
        MEMO(l,k) = min(arry);
        MEMO_longarray(j) =  MEMO(l,k);     
        
        k = k+1;
        if k>N
            k = 1;
            l = l+1;
        end
    end
    
    if opt.C>0
        for i=1:opt.C
           G(G(:,i)<0,i)=0;
        end
        G = sum(G, 2);
    
        index = find(G<=0);
        feasibleASF = MEMO_longarray(index);
        fmax = max(feasibleASF);
        if isempty(fmax)
            fmax = max(G);
        end

        index = find(G>0);
        if ~isempty(index)
            temp_asf = fmax + MEMO_longarray(index);
            MEMO_longarray(index) = temp_asf;
        end

        l = 1;
        for i=1:N
            for j=1:N 
                MEMO(i,j) = MEMO_longarray(l);
                l = l+1;
            end
        end
    end
    
    
    %----------------compute objective and ASF for Data--------------------
    
    for k=1:size(objective_data,1)
        s = xvar(k,1:2);
        [obj,~] =  high_fidelity_evaluation(opt, s);
        arry = [];
        for p=1:size(opt.dirs,1)
            w = opt.dirs(p,:);
            arry(p) = max((obj(1)-z(1))/w(1),( obj(2)-z(2))/w(2));
        end
        MEMO_data(k) = min(arry);
    end
    
    %-------------------sequentially build model space---------------------
    %%{
    for i=20:20:size(objective_data,1)    
        
        %------------------train-------------------------------------------
        asf = MEMO_data(1:i);
        if opt.C>0 
            cv = cv_data(1:i);
            index = find(cv<=0);
            feasibleASF = asf(index);
            fmax = max(feasibleASF);
            if isempty(fmax)
                fmax = max(asf);
            end

            index = find(cv>0);
            if ~isempty(index)
                temp_asf = fmax + cv(index);
                asf(index) = temp_asf;
            end
        end
        net = train_neural_network(xvar(1:i,1:2), asf);
        
        %----------------predicted space-----------------------------------
        Z = net(sol');
        Z = Z';
        Z = reshape(Z, N, N);
        Z = Z';
        
        figure;
        surfc(X,Y,Z);
        zlabel('Selection(x)');
        title(opt.objfunction);
        xlabel('x1');
        ylabel('x2');
        set(gca,'fontsize',18);
        drawnow;


        figure;
        xlim([0 .5]); 
        surfc(F1,F2,Z);
        zlabel('Selection(x)');
        title(opt.objfunction);
        xlabel('F1');
        ylabel('F2');
        set(gca,'fontsize',18);
        drawnow;
    end
    %}
    

    
    %{
    figure;
    surfc(X,Y,MEMO);
    zlabel('Selection(x)');
    title('MEMO');
    xlabel('x1');
    ylabel('x2');
    set(gca,'fontsize',18);
    drawnow;
    
    figure;
    surfc(F1,F2,MEMO);
    zlabel('Selection(x)');
    title('MEMO');
    xlabel('F1');
    ylabel('F2');
    set(gca,'fontsize',18);
    drawnow;
    %}


end
