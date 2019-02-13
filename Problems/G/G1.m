function varargout = G1(Operation, Global, input)
% Proteek Chandan Roy 2018, email:proteek_buet@yahoo.com

    switch Operation
        case 'init'
            Global.M        = 1;
            Global.D        = 13;
            Global.C        = 9;
            
            bound = zeros(2, Global.D);
            bound(2,1:9) = ones(1,9);
            bound(2,10:12) = 100*ones(1,3);
            bound(2,13) = 1;
                    
            Global.lower    = bound(1,:);
            Global.upper    = bound(2,:);
            Global.operator = @EAreal;
            
            PopDec    = rand(input,Global.D).*repmat(Global.upper-Global.lower,input,1) + repmat(Global.lower,input,1);
            varargout = {PopDec};
        case 'value'
            x = input;
            x1 = x(:,1:4); 
            x2 = x(:,5:13);
            PopObj = 5*sum(x1,2)-5*sum(x1.*x1,2)-sum(x2,2);

            % Constraints
            y(:,1) = 2*x(:,1)+2*x(:,2)+x(:,10)+x(:,11)-10;
            y(:,2) = 2*x(:,1)+2*x(:,3)+x(:,10)+x(:,12)-10;
            y(:,3) = 2*x(:,2)+2*x(:,3)+x(:,11)+x(:,12)-10;
            y(:,4) = -8*x(:,1)+x(:,10);
            y(:,5) = -8*x(:,2)+x(:,11);
            y(:,6) = -8*x(:,3)+x(:,12);
            y(:,7) = -2*x(:,4)-x(:,5)+x(:,10);
            y(:,8) = -2*x(:,6)-x(:,7)+x(:,11);
            y(:,9) = -2*x(:,8)-x(:,9)+x(:,12);
            PopCon = y; 
            varargout = {input,PopObj,PopCon};
        case 'PF'
            f = -15;
            varargout = {f};
        case 'PS'
            varargout  = {[1 1 1 1 1 1 1 1 1 3 3 3 1]}; 
    end
end