function [f]=KKT(opt, x)


        nvar = opt.V; % for ZDT1&ZDT2&ZDT3& ZDT4 (30 variables) problem
        
        nobj = opt.M; 
        rho = 0.01;


        y = []; fv=[]; df=[]; g=[]; dg=[]; um=[]; uj=[];v=[];
 
 
        %nvar = 30; ncons = 60; nobj = 2; % for ZDT1&ZDT2&ZDT3& ZDT4 (30 variables) problem
        %nvar = 10; ncons = 20; nobj = 2; %  ZDT6
        %nvar = 2; ncons = 6; nobj = 2; % for BNH , TNK, SRN
        %nvar = 6; ncons = 18; nobj = 2; % for OSY
        %nvar =5; ncons =10; nobj = 3; % for DTLZ 1_ 3 objective
        %nvar =8; ncons =16; nobj = 3; % for DTLZ 2 _3 objective
        %nvar =7; ncons =14; nobj = 5; % for DTLZ 1 5 objective
        %nvar =14; ncons =28; nobj = 5; % for DTLZ 2 5 objective
        %nvar =12; ncons =24; nobj = 10; % for DTLZ 1 10 objective 
        %nvar =19; ncons =38; nobj = 10; % for DTLZ 2 10 objective
        %nvar =12; ncons =24; nobj = 3; % for DTLZ5 objective

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
        problem = lower(opt.objfunction);
        %disp(problem)
        switch(problem)
            case 'zdt1'
                ncons = 2*nvar; 
                numVar = length(x);
                fv(1) = x(1);
                df(:,1)= [1; zeros(nvar-1,1)];
                gg = 1 + 9*sum(x(2:numVar))/(numVar-1);
                fv(2) = (1.0-sqrt(x(1)/gg))*gg;
                df(:,2) = [-0.5*sqrt(gg/x(1));(1-0.5*sqrt(x(1)/gg))*(9/(nvar-1))*ones(nvar-1,1)];
                g = [-x(1:nvar),x(1:nvar)-1];
                dg = [-eye(nvar),eye(nvar)];
                z =  [-0.001; -0.001];
            case 'zdt2' 
                ncons = 2*nvar; 
                fv(1) = x(1);   
                df(:,1)= [1; zeros(nvar-1,1)];
                yy=sum((x(2:nvar)));
                gg = 1.0+(9/(nvar-1))*yy;
                fv(2) = (1.0-(x(1)/gg)^2)*gg;
                df(:,2) = [-2*x(1)/gg;(1+x(1)*x(1)/(gg*gg))*(9/(nvar-1))*ones(nvar-1,1)]; 
                g = [-x(1:nvar),x(1:nvar)-1];
                dg = [-eye(nvar),eye(nvar)];
                z=[-0.0001; -0.0001];
            case 'zdt3'
                ncons = 2*nvar; 
                fv(1) = x(1);
                df(:,1)= [1; zeros(nvar-1,1)];
                yy=sum((x(2:nvar)));
                gg = 1.0+(9/(nvar-1))*yy;
                fv(2) = (1.0-sqrt(x(1)/gg)-(x(1)/gg)*sin(10*pi*x(1)))*gg; 
                df(:,2) = [-0.5*sqrt(gg/x(1))-sin(10*pi*x(1))-10*pi*x(1)*cos(10*pi*x(1)); (1-0.5*sqrt(x(1)/gg))*(9/(nvar-1))*ones(nvar-1,1)]; 
                g = [-x(1:nvar),x(1:nvar)-1];
                dg = [-eye(nvar),eye(nvar)];
                z=[-0.0001;-1.0001];
            case 'zdt4'
                ncons = 2*nvar; 
                fv(1) = x(1); 
                df(:,1)= [1; zeros(nvar-1,1)];
                yy=10*(nvar-1)+ sum((x(2:nvar)).^2 -10* cos(4*pi*(x(2:nvar))));
                gg = 1.0+yy;
                fv(2) = (1.0-sqrt(x(1)/gg))*gg;
                df(:,2)=[-0.5*sqrt(gg/x(1));(1.0 -0.5*x(1)/sqrt(x(1)*gg))*(2*x(2:nvar)+40*pi*sin(4*pi*x(2:nvar)))'.*ones(nvar-1,1)];
                g = [-x(1),-1-x(2:nvar),x(1:nvar)-1];
                dg = [-eye(nvar),eye(nvar)];
                z=[-0.0001;-0.0001];
            case 'zdt6'
                ncons = 2*nvar; 
                fv(1) =1-exp(-4.0*x(1))*(sin(6*pi*x(1)))^6; 
                s=4*exp(-4*x(1))*(sin(6*pi*x(1)))^6-36*pi*exp(-4*x(1))*(sin(6*pi*x(1)))^5*cos(6*pi*x(1));
                df(:,1)= [s; zeros(nvar-1,1)];
                yy=sum((x(2:nvar)));
                gg = 1.0+9.0*(yy/9.0);%^0.25
                fv(2) = (1.0-(fv(1)/gg)^2)*gg; 
                df(:,2) =[-2.0*(fv(1)*s)/gg;(1+(fv(1)*fv(1))/gg^2)*ones(nvar-1,1)];%*0.25*((yy/9.0)^(-0.75))
                g = [-x(1:nvar),x(1:nvar)-1];
                dg = [-eye(nvar),eye(nvar)];
                z=[0.2;-0.01];
            case 'dtlz1'  
                ncons = 2*nvar; 
                xv = x(:);
                l = 100*((nvar-nobj+1) + sum((xv(nobj:nvar)-0.5).^2 - cos(20*pi*(xv(nobj:nvar)-0.5))));
                for jj=1:nobj
                    if (jj==nobj)
                        prodval = 1.0;
                    else
                        prodval = prod(xv(1:nobj-jj));
                    end;

                    fv(jj)=0.5*(1.0+l)*prodval*(1.0-(jj>1)*xv(nobj-jj+1));
                    df(:,jj) = zeros(nvar,1);
            
                    df(nobj:nvar,jj)=0.5*prodval*(1.0-(jj>1)*xv(nobj-jj+1))*100*(2.0*(xv(nobj:nvar)'-0.5)+20*pi*sin(20*pi*(xv(nobj:nvar)'-0.5)));
                    for k=1:(nobj-1-(jj>1)*(jj-2))
                        yv = xv;
                        yv(k) = 1;
                        if (jj==nobj)
                            prodval = 1.0;
                        else
                            prodval = prod(yv(1:nobj-jj));
                        end;
                 
                        df(k,jj)=0.5*(1+l)*prodval*((k==nobj-jj+1)*(-1.0)+(k<nobj-jj+1)*(1.0-(jj>1)*xv(nobj-jj+1)));
                    end;
                end
                g = [-xv(1:nvar);xv(1:nvar)-1];
                g = g';
                dg = [-eye(nvar),eye(nvar)];
                z = -0.05*ones(nobj,1);

            case 'dtlz2'
                ncons = 2*nvar; 
                xv = x(:);
                l = sum((xv(nobj:nvar)-0.5).^2);
                for jj=1:nobj
                    if (jj==nobj)
                        prodval = 1.0;
                    else
                        prodval = prod(cos(0.5*pi*xv(1:nobj-jj)));
                    end;
                    fv(jj)=(1.0+l)*prodval*((jj==1)+(jj>1)*sin(0.5*pi*xv(nobj-jj+1)));
                    df(:,jj) = zeros(nvar,1);
                    df(nobj:nvar,jj)=prodval*((jj==1)+(jj>1)*sin(0.5*pi*xv(nobj-jj+1)))*(2.0*(xv(nobj:nvar)'-0.5));
                    for k=1:(nobj-1-(jj>1)*(jj-2))
                        yv=[];
                        yv = cos(0.5*pi*xv(1:nobj-jj)); % direct values
                        if (k<=nobj-jj)
                            yv(k) = -0.5*pi*sin(0.5*pi*xv(k)); % derivative of cos(x_k) terms
                        else % k=nobj-jj+1
                            yv(k) =  0.5*pi*cos(0.5*pi*xv(k)); % derivative of the extreme sine term
                        end;
                        df(k,jj)=(1.0+l)*prod(yv)*((k==nobj-jj+1)+(k<nobj-jj+1)*((jj==1)+(jj>1)*sin(0.5*pi*xv(nobj-jj+1))));
                    end;
                end
                g = [-xv(1:nvar);xv(1:nvar)-1];
                g = g';
                dg = [-eye(nvar),eye(nvar)];
                z = -0.05*ones(nobj,1);   
            case 'c2dtlz2'
                ncons = 2*nvar+1; 
                xv = x(:);
                l = sum((xv(nobj:nvar)-0.5).^2);
                for jj=1:nobj
                    if (jj==nobj)
                        prodval = 1.0;
                    else
                        prodval = prod(cos(0.5*pi*xv(1:nobj-jj)));
                    end;
                    fv(jj)=(1.0+l)*prodval*((jj==1)+(jj>1)*sin(0.5*pi*xv(nobj-jj+1)));
                    df(:,jj) = zeros(nvar,1);
                    df(nobj:nvar,jj)=prodval*((jj==1)+(jj>1)*sin(0.5*pi*xv(nobj-jj+1)))*(2.0*(xv(nobj:nvar)'-0.5));
                    for k=1:(nobj-1-(jj>1)*(jj-2))
                        yv=[];
                        yv = cos(0.5*pi*xv(1:nobj-jj)); % direct values
                        if (k<=nobj-jj)
                            yv(k) = -0.5*pi*sin(0.5*pi*xv(k)); % derivative of cos(x_k) terms
                        else % k=nobj-jj+1
                            yv(k) =  0.5*pi*cos(0.5*pi*xv(k)); % derivative of the extreme sine term
                        end;
                        df(k,jj)=(1.0+l)*prod(yv)*((k==nobj-jj+1)+(k<nobj-jj+1)*((jj==1)+(jj>1)*sin(0.5*pi*xv(nobj-jj+1))));
                    end;
                end
                [~,G] = high_fidelity_evaluation(opt, xv);
                [~,gradG] = numericGradient(xv,opt);
                g = [-xv(1:nvar);xv(1:nvar)-1];
                g = vertcat(g,G);%added by proteek
                g = g';
                dg = [-eye(nvar),eye(nvar)];
                dg = horzcat(dg,gradG);
                z = -0.05*ones(nobj,1);
                
            case 'dtlz3'
                ncons = 2*nvar; 
            case 'dtlz4'
                ncons = 2*nvar; 
                n = length(x);
                xg=x(3:n);
                gx=sum((xg-0.5).^2);
                fv(1)=(1+gx)*cos((x(1)^100)*0.5*pi)*cos((x(2)^100)*0.5*pi);
                fv(2)=(1+gx)*cos((x(1)^100)*0.5*pi)*sin((x(2)^100)*0.5*pi);
                fv(3)=(1+gx)*sin((x(1)^100)*0.5*pi);
                %%%%%%%%%%%

                df(:,1)= [-50*pi*x(1)^99*(sin(x(1)^100*0.5*pi)*cos(x(2)^100*0.5*pi)*(1+gx));...
                    -50*pi*x(2)^99*(cos(x(1)^100*0.5*pi)*sin(x(2)^100*0.5*pi)*(1+gx));...
                    2*(x(3)-0.5)*cos(x(1)^100*0.5*pi)*cos(x(2)^100*0.5*pi);...
                    2*(x(4)-0.5)*cos(x(1)^100*0.5*pi)*cos(x(2)^100*0.5*pi);...
                    2*(x(5)-0.5)*cos(x(1)^100*0.5*pi)*cos(x(2)^100*0.5*pi);...
                    2*(x(6)-0.5)*cos(x(1)^100*0.5*pi)*cos(x(2)^100*0.5*pi);...
                    2*(x(7)-0.5)*cos(x(1)^100*0.5*pi)*cos(x(2)^100*0.5*pi)];

                df(:,2)= [-50*pi*x(1)^99*(sin(x(1)^100*0.5*pi)*sin(x(2)^100*0.5*pi)*(1+gx));...
                    50*pi*x(2)^99*(cos(x(1)^100*0.5*pi)*cos(x(2)^100*0.5*pi)*(1+gx));...
                    2*(x(3)-0.5)*cos(x(1)^100*0.5*pi)*sin(x(2)^100*0.5*pi);...
                    2*(x(4)-0.5)*cos(x(1)^100*0.5*pi)*sin(x(2)^100*0.5*pi);...
                    2*(x(5)-0.5)*cos(x(1)^100*0.5*pi)*sin(x(2)^100*0.5*pi);...
                    2*(x(6)-0.5)*cos(x(1)^100*0.5*pi)*sin(x(2)^100*0.5*pi);...
                    2*(x(7)-0.5)*cos(x(1)^100*0.5*pi)*sin(x(2)^100*0.5*pi)];

                df(:,3)= [50*pi*x(1)^99*cos(x(1)^100*0.5*pi)*(1+gx);...
                0;...
                2*(x(3)-0.5)*sin((x(1)^100)*0.5*pi);...
                2*(x(4)-0.5)*sin((x(1)^100)*0.5*pi);...
                2*(x(5)-0.5)*sin((x(1)^100)*0.5*pi);...
                2*(x(6)-0.5)*sin((x(1)^100)*0.5*pi);...
                2*(x(7)-0.5)*sin((x(1)^100)*0.5*pi)];
                               
               
                
                xv = x(:);
                g = [-xv(1:nvar);xv(1:nvar)-1];
                g = g';
                dg = [-eye(nvar),eye(nvar)];
                z = [-0.05;-0.05;-0.05];
                %disp(df)
            case 'dtlz5'
                ncons = 2*nvar; 
                %%%%%%%%%%%%%% DTLZ5 -3objective          
                m = 1 + sum((x(3:end)-0.5).^2);
                fv(1)=m*cos(0.5*x(1)*pi)*cos(pi*(1+2*(m-1)*x(2))/(m*4));
                fv(2)=m*cos(0.5*x(1)*pi)*sin(pi*(1+2*(m-1)*x(2))/(m*4));
                fv(3)=m*sin(0.5*x(1)*pi);
                temp_df(:,1) = [-m*0.5*pi*sin(0.5*x(1)*pi)*cos((1/4)*pi*(1+(2*(m-1))*x(2))/m);-(1/4)*cos(0.5*x(1)*pi)*sin((1/4)*pi*(1+(2*(m-1))*x(2))/m)*pi*(2*m-2);(2*x(3)-1.0)*cos(0.5*x(1)*pi)*cos(pi*(1+2*(m-1)*x(2))/(m*4))-m*cos(0.5*x(1)*pi)*((1/2)*pi*(2*x(3)-1.0)*x(2)/m-(1/4)*pi*(1+(2*(m-1))*x(2))*(2*x(3)-1.0)/m^2)*sin(pi*(1+2*(m-1)*x(2))/(m*4));(2*x(4)-1.0)*cos(0.5*x(1)*pi)*cos(pi*(1+2*(m-1)*x(2))/(m*4))-m*cos(0.5*x(1)*pi)*((1/2)*pi*(2*x(4)-1.0)*x(2)/m-(1/4)*pi*(1+(2*(m-1))*x(2))*(2*x(4)-1.0)/m^2)*sin(pi*(1+2*(m-1)*x(2))/(m*4));(2*x(5)-1.0)*cos(0.5*x(1)*pi)*cos(pi*(1+2*(m-1)*x(2))/(m*4))-m*cos(0.5*x(1)*pi)*((1/2)*pi*(2*x(5)-1.0)*x(2)/m-(1/4)*pi*(1+(2*(m-1))*x(2))*(2*x(5)-1.0)/m^2)*sin(pi*(1+2*(m-1)*x(2))/(m*4));(2*x(6)-1.0)*cos(0.5*x(1)*pi)*cos(pi*(1+2*(m-1)*x(2))/(m*4))-m*cos(0.5*x(1)*pi)*((1/2)*pi*(2*x(6)-1.0)*x(2)/m-(1/4)*pi*(1+(2*(m-1))*x(2))*(2*x(6)-1.0)/m^2)*sin(pi*(1+2*(m-1)*x(2))/(m*4));(2*x(7)-1.0)*cos(0.5*x(1)*pi)*cos(pi*(1+2*(m-1)*x(2))/(m*4))-m*cos(0.5*x(1)*pi)*((1/2)*pi*(2*x(7)-1.0)*x(2)/m-(1/4)*pi*(1+(2*(m-1))*x(2))*(2*x(7)-1.0)/m^2)*sin(pi*(1+2*(m-1)*x(2))/(m*4))];
                temp_df(:,2) = [-m*0.5*pi*sin(0.5*x(1)*pi)*sin((1/4)*pi*(1+(2*(m-1))*x(2))/m);(1/4)*cos(0.5*x(1)*pi)*cos((1/4)*pi*(1+(2*(m-1))*x(2))/m)*pi*(2*m-2);(2*x(3)-1.0)*cos(0.5*x(1)*pi)*sin(pi*(1+2*(m-1)*x(2))/(m*4))+m*cos(0.5*x(1)*pi)*((1/2)*pi*(2*x(3)-1.0)*x(2)/m-(1/4)*pi*(1+(2*(m-1))*x(2))*(2*x(3)-1.0)/m^2)*cos(pi*(1+2*(m-1)*x(2))/(m*4));(2*x(4)-1.0)*cos(0.5*x(1)*pi)*sin(pi*(1+2*(m-1)*x(2))/(m*4))+m*cos(0.5*x(1)*pi)*((1/2)*pi*(2*x(4)-1.0)*x(2)/m-(1/4)*pi*(1+(2*(m-1))*x(2))*(2*x(4)-1.0)/m^2)*cos(pi*(1+2*(m-1)*x(2))/(m*4));(2*x(5)-1.0)*cos(0.5*x(1)*pi)*sin(pi*(1+2*(m-1)*x(2))/(m*4))+m*cos(0.5*x(1)*pi)*((1/2)*pi*(2*x(5)-1.0)*x(2)/m-(1/4)*pi*(1+(2*(m-1))*x(2))*(2*x(5)-1.0)/m^2)*cos(pi*(1+2*(m-1)*x(2))/(m*4));(2*x(6)-1.0)*cos(0.5*x(1)*pi)*sin(pi*(1+2*(m-1)*x(2))/(m*4))+m*cos(0.5*x(1)*pi)*((1/2)*pi*(2*x(6)-1.0)*x(2)/m-(1/4)*pi*(1+(2*(m-1))*x(2))*(2*x(6)-1.0)/m^2)*cos(pi*(1+2*(m-1)*x(2))/(m*4));(2*x(7)-1.0)*cos(0.5*x(1)*pi)*sin(pi*(1+2*(m-1)*x(2))/(m*4))+m*cos(0.5*x(1)*pi)*((1/2)*pi*(2*x(7)-1.0)*x(2)/m-(1/4)*pi*(1+(2*(m-1))*x(2))*(2*x(7)-1.0)/m^2)*cos(pi*(1+2*(m-1)*x(2))/(m*4))];
                temp_df(:,3) = [m*0.5*pi*cos(0.5*x(1)*pi);0;(2*x(3)-1.0)*sin(0.5*x(1)*pi);(2*x(4)-1.0)*sin(0.5*x(1)*pi);(2*x(5)-1.0)*sin(0.5*x(1)*pi);(2*x(6)-1.0)*sin(0.5*x(1)*pi);(2*x(7)-1.0)*sin(0.5*x(1)*pi)];
                df(1:7,1) = temp_df(1:7,1);
                df(1:7,2) = temp_df(1:7,2);
                df(1:7,3) = temp_df(1:7,3);
                
                xv = x(:);
                g = [-xv(1:nvar);xv(1:nvar)-1];
                g = g';
                dg = [-eye(nvar),eye(nvar)];
                z = [-0.05;-0.05;-0.05];
            case 'dtlz7'
                ncons = 2*nvar; 
            case 'bnh'  
                ncons = 6; 
                fv(1) = 4*x(1)^2+4*x(2)^2;
                fv(2) = (x(1)-5)^2+(x(2)-5)^2;
                g(1) = (x(1)-5)^2+x(2)^2-25;
                g(2) = -(x(1)-8)^2-(x(2)+3)^2+7.7;
                g(3) =-x(1);
                g(4) =-x(2);
                g(5) =x(1)-5;
                g(6) =x(2)-3;
                df(:,1)=[8*x(1);8*x(2)];
                df(:,2)=[2*(x(1)-5);2*(x(2)-5)];
                dg(:,1)=[2*(x(1)-5);2*x(2)];
                dg(:,2) =[-2*(x(1)-8);-2*(x(2)+3)];
                dg(:,3) = [-1;0];
                dg(:,4) = [0;-1];
                dg(:,5) = [1;0];
                dg(:,6) = [0;1];
                z = [-0.05;-0.05];            
            case 'osy'
                ncons = 18; 
                fv(1) = -(25*(x(1)-2)^2+(x(2)-2)^2+(x(3)-1)^2+(x(4)-4)^2+(x(5)-1)^2);
                fv(2) = x(1)^2+x(2)^2+x(3)^2+x(4)^2+x(5)^2+x(6)^2;
                g(1) = -x(1)-x(2)+2;
                g(2) = x(1)+x(2)-6;
                g(3) = -x(1)+x(2)-2;
                g(4) =x(1)-3*x(2)-2;
                g(5) = -4 + (x(3)-3)^2 +x(4);
                g(6) = -(x(5)-3)^2-x(6)+4;
                g(7) = -x(4);
                g(8) = -x(6);
                g(9) = 1-x(5);
                g(10) = 1-x(3);
                g(11) = -x(1);
                g(12) = -x(2);
                g(13) = x(4)-6;
                g(14) = x(6)-10;
                g(15) = x(5)-5;
                g(16) =x(3)-5;
                g(17) =x(1)-10;
                g(18) =x(2)-10;
                %%%% Derivatives here
                df(:,1) = [-50*(x(1)-2); -2*(x(2)-2); -2*(x(3)-1); -2*(x(4)-4); -2*(x(5)-1);0];
                df(:,2) = [2*x(1); 2*x(2); 2*x(3); 2*x(4); 2*x(5);2*x(6)];
                dg(:,1) = [-1; -1; 0; 0; 0; 0];
                dg(:,2) = [1; 1; 0; 0; 0; 0];
                dg(:,3) = [-1; 1; 0; 0; 0; 0];
                dg(:,4) = [1; -3; 0; 0; 0; 0];
                dg(:,5) = [0; 0; 2*(x(3)-3); 1; 0; 0];
                dg(:,6) = [0; 0; 0; 0; -2*(x(5)-3); -1];
                dg(:,7) = [0; 0; 0; -1; 0; 0];
                dg(:,8) = [0; 0; 0; 0; 0; -1];
                dg(:,9) = [0; 0; 0; 0; -1; 0];
                dg(:,10)= [0; 0; -1; 0; 0; 0];
                dg(:,11)= [-1; 0; 0; 0; 0; 0];
                dg(:,12)= [0; -1; 0; 0; 0; 0];
                dg(:,13)= [0; 0; 0; 1; 0; 0];
                dg(:,14)= [0; 0; 0; 0; 0; 1];
                dg(:,15)= [0; 0; 0; 0; 1; 0];
                dg(:,16)= [0; 0; 1; 0; 0; 0];
                dg(:,17)= [1; 0; 0; 0; 0; 0];
                dg(:,18)= [0; 1; 0; 0; 0; 0];
                z = [-300; 0];
            case 'srn'
                ncons = 6; 
                fv(1) = 2.0+(x(1)-2.0)^2+(x(2)-1.0)^2;
                fv(2) = 9.0*x(1)-(x(2)-1)^2;
                % Constraints
                g(1) = (x(1)^2+x(2)^2-225.0)/225.0;
                g(2) = (x(1)-3.0*x(2)+10.0)/10.0;
                g(3) =-x(1)-20.0;
                g(4) =x(1)-20.0;
                g(5) =-x(2)-20.0;
                g(6) =x(2)-20.0;
                % Derivatives of objectives
                df(:,1)=[2.0*(x(1)-2.0); 2*(x(2)-1.0)];
                df(:,2)=[9.0; -2.0*(x(2)-1.0)];
                % Derivatives of constraints
                dg(:,1)=[2*x(1); 2*x(2)]./225.0;
                dg(:,2)=[1.0; -3.0]./10.0;
                dg(:,3)=[-1.0; 0.0];
                dg(:,4)=[1.0; 0.0];
                dg(:,5)=[0.0; -1.0];
                dg(:,6)=[0.0; 1.0];
                % Ideal Point
                z = [0.0; -300];
            case 'tnk'
                ncons = 6;
                fv(1) = x(1);
                fv(2) = x(2);
                g(1) = -x(1)^2-x(2)^2+1.0+0.1*cos(16*atan(x(1)/x(2)));
                g(2) = (x(1)-0.5)^2+(x(2)-0.5)^2-0.5;
                g(3) =-x(1);
                g(4) =-x(2);
                g(5) = x(1)-3.15;
                g(6) =x(2)-3.15;
                df(:,1)=[1; 0];
                df(:,2)=[0; 1];
                dg(:,1)=[-2*x(1)-1.6*sin(16*atan(x(1)/x(2)))*x(2)/(x(2)^2+x(1)^2); -2*x(2)+1.6*sin(16*atan(x(1)/x(2)))*x(1)/(x(2)^2+x(1)^2)];
                dg(:,2)=[2*(x(1)-0.5);2*(x(2)-0.5)];
                dg(:,3) = [-1; 0];
                dg(:,4) = [0; -1];
                dg(:,5) = [1; 0];
                dg(:,6) = [0; 1];
                z = [-0.01;-0.01];
            case 'carside'%%%%%%%% Car impact 
                ncons = 7*2+10; %number of variable*2 + number of constraints 
                [fv, g2] = high_fidelity_evaluation(opt,x);
                g(1) = g2(1);
                g(2) = g2(2);
                g(3) = g2(3);
                g(4) = g2(4);
                g(5) = g2(5);
                g(6) = g2(6);
                g(7) = g2(7);
                g(8) = g2(8);
                g(9) = g2(9);
                g(10) = g2(10);
                g(11) = -x(1)-0.5;
                g(12) = -x(2)-0.45;
                g(13) = -x(3)-0.5;
                g(14) = -x(4)-0.5;
                g(15) = -x(5)-0.875;
                g(16) = -x(6)-0.4;
                g(17) = -x(7)-0.4;
                g(18) = x(1)-1.5;
                g(19) = x(2)-1.35;
                g(20) = x(3)-1.5;
                g(21) = x(4)-1.5;
                g(22) = x(5)-2.625;
                g(23) = x(6)-1.2;
                g(24) = x(7)-1.2;
                
                yy0 = 0.345;
                yy1 = 0.192;
                df(:,1) = [4.9;6.67;6.98;4.01;1.78;0.00001;2.73];
                df(:,2) = [0;- 0.19 *x(3);- 0.19 *x(2);-0.5;0;0;0];
                df(:,3) = [- 0.674 *0.5*x(2);- 0.674 * x(1)*0.5- 1.95 *0.5* yy0;- 0.489*0.5*x(7);0;- 0.843 *0.5*x(6);- 0.843 *x(5)*0.5;- 0.489 *x(3)*0.5];
                dg(:,1) = [0;- 0.3717*x(4);- 0.0092928;- 0.3717 *x(2);0;0;0];
                dg(:,2) = [(-0.0159*x(2)-0.188*x(1)*yy0)/0.32;(-0.0159*x(1)-0.019*x(7))/0.32;0.0144*x(5)/0.32;0;0.0144*x(3)/0.32;0.08045*yy1/0.32;-0.019*x(2)/0.32];
                dg(:,3) = [(-0.131*yy0-0.0704*yy1)/0.32;(0.03099*x(6)-0.018*x(7)-0.018*x(2)*2)/0.32;(0.0208*yy0+0.121*yy1)/0.32;0;(0.00817-0.00364*x(6))/0.32;(0.03099*x(2)-0.00364*x(5))/0.32;-0.018*x(2)/0.32];
                dg(:,4) = [0;(-0.61+ 0.227 *x(2)*2)/0.32;-0.031296/0.32;0;0;0;-0.166*yy1/0.32];
                dg(:,5) = [- 4.2*x(2)/32.0;-4.2*x(1)/32.0;3.818/32.0;0;0;6.63*yy1/32.0;-7.77*yy0/32.0];
                dg(:,6) = [-5.057*x(2)/32.0;(-5.057*x(1)-11*yy0)/32.0;2.95/32.0;0;0;0;-9.98*yy0/32.0];
                dg(:,7) = [- 12.9 *yy0/32.0;- 9.9/32.0;0;0;0;0;0];
                dg(:,8) = [0;- 0.19*x(3)/4.0;- 0.19 *x(2)/4.0;-0.5/4.0;0;0;0];
                dg(:,9) = [- 0.674*x(2)/9.9;(- 0.674*x(1)-1.95*yy0)/9.9;0;0;0;0;0];
                dg(:,10)= [0;0;-0.489*x(7)/15.7;0;- 0.843*x(6)/15.7;- 0.843 *x(5)/15.7;- 0.489*x(3)/15.7];
                
                dg(:,11)= [-1.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0];
                dg(:,12)= [0.0; -1.0; 0.0; 0.0; 0.0; 0.0; 0.0];
                dg(:,13)= [0.0; 0.0; -1.0; 0.0; 0.0; 0.0; 0.0];
                dg(:,14)= [0.0; 0.0; 0.0; -1.0; 0.0; 0.0; 0.0];
                dg(:,15)= [0.0; 0.0; 0.0; 0.0; -1.0; 0.0; 0.0];
                dg(:,16)= [0.0; 0.0; 0.0; 0.0; 0.0; -1.0; 0.0];
                dg(:,17)= [0.0; 0.0; 0.0; 0.0; 0.0; 0.0; -1.0];
                dg(:,18)= [1.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0];
                dg(:,19)= [0.0; 1.0; 0.0; 0.0; 0.0; 0.0; 0.0];
                dg(:,20)= [0.0; 0.0; 1.0; 0.0; 0.0; 0.0; 0.0];
                dg(:,21)= [0.0; 0.0; 0.0; 1.0; 0.0; 0.0; 0.0];
                dg(:,22)= [0.0; 0.0; 0.0; 0.0; 1.0; 0.0; 0.0];
                dg(:,23)= [0.0; 0.0; 0.0; 0.0; 0.0; 1.0; 0.0];
                dg(:,24)= [0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 1.0];
                
                z = ([24.3180    3.5352   10.5610])';
            case 'welded'%%% Welded beam
                ncons = 4*2+4;%number of variable*2 + number of constraints 
                [fv, g2] = high_fidelity_evaluation(opt,x);
                g(1) = g2(1);
                g(2) = g2(2);
                g(3) = g2(3);
                g(4) = g2(4);
                g(5) = -x(1)-0.125;
                g(6) = -x(2)-0.1;
                g(7) = -x(3)-0.1;
                g(8) = -x(4)-0.125;
                g(9) = x(1)-5;
                g(10) = x(2)-10;
                g(11) = x(3)-10;
                g(12) = x(4)-5;
                
                df(:,1) = [2.20942* x(1) *x(2);1.10471 * x(1)^2+ 0.04811 *x(3)* x(4) ;0.04811*x(4)*(14.0 + x(2));0.04811* x(3) *(14.0 + x(2))];
                df(:,2) = [0;0;-6.585600000/(x(4)*x(3)^4);-2.195200000/(x(4)^2*x(3)^3)];
                dg(:,1) = [(0.3676470588e-4*(-3.600000002*10^7/(x(1)^3*x(2)^2)+.4999999996*(84000.00+3000.000000*x(2))^2*...
                    (.50*x(3)+.50*x(1))/(x(1)^2*x(2)^2*(0.8333333333e-1*x(2)^2+(.5*x(3)+.5*x(1))^2)^2)-.9999999992*(84000.00+3000.000000*x(2))^2*...
                    (.2500000000*x(2)^2+(.5*x(3)+.5*x(1))^2)/(x(1)^3*x(2)^2*(0.8333333333e-1*x(2)^2+(.5*x(3)+.5*x(1))^2)^2)-.9999999992*...
                    (84000.00+3000.000000*x(2))^2*(.2500000000*x(2)^2+(.5*x(3)+.5*x(1))^2)*(.50*x(3)+.50*x(1))/(x(1)^2*x(2)^2*(0.8333333333e-1*x(2)^2+(.5*x(3)+...
                    .5*x(1))^2)^3)-(6000.000000*(84000.00+3000.000000*x(2)))/(x(1)^3*x(2)*(0.8333333333e-1*x(2)^2+(.5*x(3)+.5*x(1))^2))-(3000.000000*(84000.00+ ...
                    3000.000000*x(2)))*(.50*x(3)+.50*x(1))/(x(1)^2*x(2)*(0.8333333333e-1*x(2)^2+(.5*x(3)+.5*x(1))^2)^2)))/sqrt(1.800000001*10^7/(x(1)^2*x(2)^2)+...
                    .4999999996*(84000.00+3000.000000*x(2))^2*(.2500000000*x(2)^2+(.5*x(3)+.5*x(1))^2)/(x(1)^2*x(2)^2*(0.8333333333e-1*x(2)^2+(.5*x(3)+.5*x(1))^2)^2)+...
                    (3000.000000*(84000.00+3000.000000*x(2)))/(x(1)^2*x(2)*(0.8333333333e-1*x(2)^2+(.5*x(3)+.5*x(1))^2)));...
                    (0.3676470588e-4*(-3.600000002*10^7/(x(1)^2*x(2)^3)+(2999.999998*(84000.00+3000.000000*x(2)))*(.2500000000*x(2)^2+(.5*x(3)+.5*x(1))^2)/(x(1)^2*x(2)^2*...
                    (0.8333333333e-1*x(2)^2+(.5*x(3)+.5*x(1))^2)^2)+.2499999998*(84000.00+3000.000000*x(2))^2/(x(2)*x(1)^2*(0.8333333333e-1*x(2)^2+(.5*x(3)+.5*x(1))^2)^2)-.9999999992*...
                    (84000.00+3000.000000*x(2))^2*(.2500000000*x(2)^2+(.5*x(3)+.5*x(1))^2)/(x(1)^2*x(2)^3*(0.8333333333e-1*x(2)^2+(.5*x(3)+.5*x(1))^2)^2)-.1666666666*...
                    (84000.00+3000.000000*x(2))^2*(.2500000000*x(2)^2+(.5*x(3)+.5*x(1))^2)/(x(1)^2*x(2)*(0.8333333333e-1*x(2)^2+(.5*x(3)+.5*x(1))^2)^3)-(3000.000000*...
                    (84000.00+3000.000000*x(2)))/(x(1)^2*x(2)^2*(0.8333333333e-1*x(2)^2+(.5*x(3)+.5*x(1))^2))+9.000000000*10^6/(x(1)^2*x(2)*(0.8333333333e-1*x(2)^2+(.5*x(3)+.5*...
                    x(1))^2))-(500.0000001*(84000.00+3000.000000*x(2)))/(x(1)^2*(0.8333333333e-1*x(2)^2+(.5*x(3)+.5*x(1))^2)^2)))/sqrt(1.800000001*10^7/(x(1)^2*x(2)^2)+.4999999996*(84000.00+...
                    3000.000000*x(2))^2*(.2500000000*x(2)^2+(.5*x(3)+.5*x(1))^2)/(x(1)^2*x(2)^2*(0.8333333333e-1*x(2)^2+(.5*x(3)+.5*x(1))^2)^2)+(3000.000000*(84000.00+3000.000000*x(2)))/(x(1)^2*...
                    x(2)*(0.8333333333e-1*x(2)^2+(.5*x(3)+.5*x(1))^2)));(0.3676470588e-4*(1.799999999*10^7*(14.0+(1/2)*x(2))^2*(.50*x(3)+.50*x(1))/(x(1)^2*x(2)^2*((1/12)*x(2)^2+(.5*x(3)+...
                    .5*x(1))^2)^2)-3.599999998*10^7*(14.0+(1/2)*x(2))^2*(.2500000000*x(2)^2+(.5*x(3)+.5*x(1))^2)*(.50*x(3)+.50*x(1))/(x(1)^2*x(2)^2*((1/12)*x(2)^2+(.5*x(3)+.5*x(1))^2)^3)-...
                    1.272792206*10^7*sqrt(2)*(14.0+(1/2)*x(2))*(.50*x(3)+.50*x(1))/(x(1)^2*x(2)*((1/12)*x(2)^2+(.5*x(3)+.5*x(1))^2)^2)))/sqrt(1.800000000*10^7/(x(1)^2*x(2)^2)+...
                    1.799999999*10^7*(14.0+(1/2)*x(2))^2*(.2500000000*x(2)^2+(.5*x(3)+.5*x(1))^2)/(x(1)^2*x(2)^2*((1/12)*x(2)^2+(.5*x(3)+.5*x(1))^2)^2)+1.272792206*10^7*sqrt(2)*(14.0+...
                    (1/2)*x(2))/(x(1)^2*x(2)*((1/12)*x(2)^2+(.5*x(3)+.5*x(1))^2)));0];
                dg(:,2) = [0;0;-33.60000000/(x(4)*x(3)^3);-16.80000000/(x(4)^2*x(3)^2)];
                dg(:,3) = [1/x(4);0;0;-x(1)/x(4)^2];
                dg(:,4) = [0;0;-(10.79100362*(1.0-0.2823462196e-1*x(3)))*x(3)*x(4)^6/sqrt(x(3)^2*x(4)^6)+.3046799078*sqrt(x(3)^2*x(4)^6);-(32.37301086*(1.0-0.2823462196e-1*x(3)))*x(3)^2*x(4)^5/sqrt(x(3)^2*x(4)^6)];
                dg(:,5) = [-1.0; 0.0; 0.0; 0.0];
                dg(:,6) = [0.0; -1.0; 0.0; 0.0];
                dg(:,7) = [0.0; 0.0; -1.0; 0.0];
                dg(:,8) = [0.0; 0.0; 0.0; -1.0];
                dg(:,9) = [1.0; 0.0; 0.0; 0.0];
                dg(:,10) = [0.0; 1.0; 0.0; 0.0];
                dg(:,11) = [0.0; 0.0; 1.0; 0.0];
                dg(:,12) = [0.0; 0.0; 0.0; 1.0];
                
                z  = ([2.3316   -0.0496])';
            otherwise
                input('Inside KKT.m, Function and Derivative defition is not found');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
     infeas=0;
     % identify if solution is feasible and compute constraint violation
     sumg = 0.0;
     if (ncons>0)
         for n=1:ncons
             if (g(n)>0) 
                 sumg = sumg + g(n)*g(n);
             end
         end
     end;
     if (sumg > 0.0)
         finfeas = 1.0+sumg;
         fvalue = finfeas;
         infeas=1; % it is an infeasible point
         y = [y, fvalue];
         f = fvalue;
     end
     %disp(sumg)
     if (infeas == 0)
         am=[]; aj=[]; gsq=[]; invaj=[]; um=[]; uj=[];
         
         % compute the weight vector
         w = (sqrt((fv'-z)'*(fv'-z)))*ones(nobj,1)./(fv'-z);
         h=w.*(fv'-z);
          A=max(h)+rho*sum(w.*(fv'-z));
         % form the A_m matrix from gradients of objectives
         rhosum = repmat(w',nvar,1).*df;
         for n=1:nobj
             am(n,:) = w(n)*df(:,n)'+rho*sum(rhosum',1);
          
         end;
         % form the A_j matrix from gradients of constraints
         % gsq is a matrix for taking care of penalty term in obj. function
         biga=[]; rb=[]; uvec=[];
         if (ncons>0)
             gsq = zeros(ncons,ncons);
             for n=1:ncons
                 aj(n,:) = dg(:,n)';
                 gsq(n,n) = gsq(n,n) + g(n)*g(n);
             end
             biga = [[ones(nobj,nobj) + am*am',am*aj']; [aj*am', aj*aj'+gsq]];
             rb = [ones(nobj,1);zeros(ncons,1)];
         else
             biga = ones(nobj,nobj) + am*am';
             rb = ones(nobj,1);
         end
         
%       uvecn = quadprog(biga,-rb,[],[],[],[],zeros(nobj+ncons,1),[]);  %  for using quadratic programming 
        
         %disp('before first pinv')
         uvec = pinv(biga)*rb;
         %disp('after first pinv')
         uvecn = uvec;
         cbiga = biga;
         crb = rb;
         negvec = lt(uvec,0); % puts 1 if vec(i) is zero
         count = 0;
         while (max(negvec(1:nobj)) == 1)
             for n=1:nobj
                 if (negvec(n) == 1)
                     cbiga(n,:) = zeros(1,nobj+ncons);
                     cbiga(:,n) = zeros(nobj+ncons,1);
                     cbiga(n,n) = 1;
                     crb(n) = 0;
                     if strcmpi(opt.objfunction,'c2dtlz2')||strcmpi(opt.objfunction,'dtlz4')||strcmpi(opt.objfunction,'dtlz5')||strcmpi(opt.objfunction,'dtlz2')
                        uvecn = inv(cbiga)*crb;
                     else
                        uvecn = pinv(cbiga)*crb;
                     end
                     
                     if uvecn(n) <0 && uvecn(n)> -1e-10
                        uvecn(n) = 0; 
                     end
                     negvec = lt(uvecn,0);
                 end;
             end;
             %if max(negvec(1:nobj))==1 && count>10
             %   disp('1')
             %end
             %disp('inside while')
             %disp(max(negvec(1:nobj)))
             count = count+1;
         end;

%%% uj calculation
         if (ncons>0)
             for n=(nobj+1):(nobj+ncons)
                 if (negvec(n) == 1)
                     cbiga(n,:) = zeros(1,nobj+ncons);
                     cbiga(:,n) = zeros(nobj+ncons,1);
                     cbiga(n,n) = 1;
                     crb(n) = 0;
                     uvecn = pinv(cbiga)*crb;
                     negvec = lt(uvecn,0);
                 end;
             end;
             %disp('inside ncons>0, pinv')
          end;  
         
         um = uvecn(1:nobj);
         if (ncons>0)
             uj = uvecn((nobj+1):(nobj+ncons));
             kktpm = (1-sum(um))^2 + sum(([am;aj]'*uvecn).^2);
            % kktpm = 1 - sum(um) - (g*uj)^2;%%%%
             fvalue = kktpm + sum((uj.*g').^2);
         else
             kktpm = (1-sum(um))^2 + sum((am'*uvecn).^2);
             fvalue = kktpm;
         end;
         %disp('before new way')
         % new way
         ujgj = -g*uj;
         leftside = sum(um)+ ujgj*(1+ujgj);
         
         %B=(kktpm+sum(um.*(w.*(fv'-z)+rho*sum(w.*(fv'-z))))+sum(uj.*g'))/sum(um);
         %if (A<=B)
         if (leftside <= 1)
            %kktpm = kktpm; 
             f=kktpm;
         else
         
            kktpm1 =-g*uj;   % calculates error from Adjusted method
            kktpmnew=(kktpm*g*g'-g*uj)/(1+g*g');  % projection 
            kktpm=(kktpm1+kktpm+kktpmnew)/3.0; %calculates mean error  between projection and the original direct method and the Adjusted method
            f=kktpm;
         end;
      end; % end of infeas==0
      %disp('kkt done')
    end


