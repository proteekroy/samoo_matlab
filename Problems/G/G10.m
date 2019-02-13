function varargout = G10(Operation, Global, input)
% Proteek Chandan Roy 2018, email:proteek_buet@yahoo.com

    switch Operation
        case 'init'
            Global.M        = 1;
            Global.D        = 8;
            Global.C        = 6;
            
            bound = zeros(2, Global.D);
            bound(1,1)  = 100;
            bound(2,1)  = 10000;
            bound(1,2:3)  = 1000;
            bound(2,2:3)  = 10000;
            bound(1,4:8)  = 10;
            bound(2,4:8)  = 1000;
                    
            Global.lower    = bound(1,:);
            Global.upper    = bound(2,:);
            Global.operator = @EAreal;
            
            PopDec    = rand(input,Global.D).*repmat(Global.upper-Global.lower,input,1) + repmat(Global.lower,input,1);
            varargout = {PopDec};
        case 'value'
            x = input;
            PopObj = x(:,1)+x(:,2)+x(:,3);

            % Constraints
            y(:,1) = -1+0.0025*(x(:,4)+x(:,6));
            y(:,2) = -1+0.0025*(-x(:,4)+x(:,5)+x(:,7));
            y(:,3) = -1+0.01*(-x(:,5)+x(:,8));
            y(:,4) = 100*x(:,1)-x(:,1).*x(:,6)+833.33252*x(:,4)-83333.333;
            y(:,5) = x(:,2).*x(:,4)-x(:,2).*x(:,7)-1250.*x(:,4)+1250.*x(:,5);
            y(:,6) = x(:,3).*x(:,5)-x(:,3).*x(:,8)-2500.*x(:,5)+1250000;
            
            PopCon = y; 
            varargout = {input,PopObj,PopCon};
        case 'PF'
            f = 7049.24802052867;
            varargout = {f};
        case 'PS'
            varargout  = {[579.306685017979589, 1359.97067807935605, 5109.97065743133317, 182.01769963061534, 295.601173702746792, 217.982300369384632, 286.41652592786852, 395.601173702746735]}; 
    end
end