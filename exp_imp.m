
function [EIF] = exp_imp(popObj, bestpopASF, MSE)

    n = size(popObj,1);
    s = sqrt(abs(MSE) );
    EIF = zeros(n,1);
    for i=1:size(popObj,1)
        yBest = bestpopASF;%ybest
        y_hat = popObj(i,:);%solution obj
        EIFa = ( yBest - y_hat) * (0.5+0.5*erf((1/sqrt(2))*((yBest-y_hat)/s(i))));
        EIFb = s(i)*(1/sqrt(2*pi)) * exp(-0.5*(yBest-y_hat)^2/abs(MSE(i)));
        EIF(i) = EIFa+EIFb;
        EIF(i) = -EIF(i); % MINIMIZE Standard EIF!
    end
end