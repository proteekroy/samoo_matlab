function varargout = G24(Operation, Global, input)
% Proteek Chandan Roy 2018, email:proteek_buet@yahoo.com

    switch Operation
        case 'info'
            Global.M        = 1;
            Global.D        = 2;
            Global.C        = 2;
            
            Global.lower    = [0, 0];
            Global.upper    = [3, 4];
            Global.operator = @EAreal;
            
            varargout = {Global};
        case 'evaluate'
            x = input;
            %Objective
            PopObj = -x(:, 1) - x(:, 2);
            % Constraints
            y(:,1) = -2*power(x(:,1),4)+8*power(x(:,1),3)-8*power(x(:,1),2)+x(:,2);
            y(:,1) = -4*power(x(:,1),4)+32*power(x(:,1),3)-88*power(x(:,1),2)+96*x(:,1)+x(:,2)-36;
            
            PopCon = y; 
            varargout = {input,PopObj,PopCon};
        case 'PF'
            f = -5.50801327159536;
            varargout = {f};
        case 'PS'
            varargout  = {[2.32952019747762 3.17849307411774]}; 
    end
end