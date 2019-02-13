function varargout = G3(Operation, Global, input)
% Proteek Chandan Roy 2018, email:proteek_buet@yahoo.com

    switch Operation
        case 'init'
            Global.M        = 1;
            Global.D        = 10;
            Global.C        = 1;
            
            bound = zeros(2, Global.D);
            bound(2,1:Global.D) = ones(1,Global.D);
                    
            Global.lower    = bound(1,:);
            Global.upper    = bound(2,:);
            Global.operator = @EAreal;
            
            PopDec    = rand(input,Global.D).*repmat(Global.upper-Global.lower,input,1) + repmat(Global.lower,input,1);
            varargout = {PopDec};
        case 'value'
            x = input;
            n = size(x,2);
            PopObj = -sqrt(n)^n*prod(x,2);

            % Constraints
            y(:,1) = abs(sum(x.^2,2)-1)-1e-4;  
            PopCon = y; 
            varargout = {input,PopObj,PopCon};
        case 'PF'
            f = -1.00050010001000;
            varargout = {f};
        case 'PS'
            varargout  = {[ 0.31624357647283069, ...
                            0.316243577414338339, 0.316243578012345927, 0.316243575664017895, 0.316243578205526066,...
                            0.31624357738855069, 0.316243575472949512, 0.316243577164883938, 0.316243578155920302,...
                            0.316243576147374916]}; 
    end
end