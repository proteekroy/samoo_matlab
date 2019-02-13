function varargout = G6(Operation, Global, input)
% Proteek Chandan Roy 2018, email:proteek_buet@yahoo.com

    switch Operation
        case 'init'
            Global.M        = 1;
            Global.D        = 2;
            Global.C        = 2;
            
            bound = zeros(2, Global.D);
            bound(1,1)  = 13;
            bound(2,1)  = 100;
            bound(1,2)  = 0;
            bound(2,2)  = 100;
                    
            Global.lower    = bound(1,:);
            Global.upper    = bound(2,:);
            Global.operator = @EAreal;
            
            PopDec    = rand(input,Global.D).*repmat(Global.upper-Global.lower,input,1) + repmat(Global.lower,input,1);
            varargout = {PopDec};
        case 'value'
            x = input;
            PopObj = (x(:,1)-10).^3+(x(:,2)-20).^3;

            % Constraints
            y(:,1) = -(x(:,1)-5).^2-(x(:,2)-5).^2+100;
            y(:,2) = (x(:,1)-6).^2+(x(:,2)-5).^2-82.81;

            PopCon = y; 
            varargout = {input,PopObj,PopCon};
        case 'PF'
            f = -6961.81387558015;
            varargout = {f};
        case 'PS'
            varargout  = {[14.09500000000000064, 0.8429607892154795668]}; 
    end
end