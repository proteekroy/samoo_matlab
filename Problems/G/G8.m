function varargout = G8(Operation, Global, input)
% Proteek Chandan Roy 2018, email:proteek_buet@yahoo.com

    switch Operation
        case 'init'
            Global.M        = 1;
            Global.D        = 2;
            Global.C        = 2;
            
            bound = zeros(2, Global.D);
            bound(2,1:Global.D) = 10;
                    
            Global.lower    = bound(1,:);
            Global.upper    = bound(2,:);
            Global.operator = @EAreal;
            
            PopDec    = rand(input,Global.D).*repmat(Global.upper-Global.lower,input,1) + repmat(Global.lower,input,1);
            varargout = {PopDec};
        case 'value'
            x = input;
            PopObj = -(sin(2*pi*x(:,1)).^3.*sin(2*pi*x(:,2)))./(x(:,1).^3.*(x(:,1)+x(:,2)));

            % Constraints
            y(:,1) = x(:,1).^2-x(:,2)+1;
            y(:,2) = 1-x(:,1)+(x(:,2)-4).^2;   
            
            PopCon = y; 
            varargout = {input,PopObj,PopCon};
        case 'PF'
            f = -0.0958250414180359;
            varargout = {f};
        case 'PS'
            varargout  = {[1.22797135260752599, 4.24537336612274885]}; 
    end
end