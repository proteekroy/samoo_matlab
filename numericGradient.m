function [gradF,gradG] = numericGradient(xv,opt)

gradF = [];
gradG = [];

opt.delh = 0.0001;

if(opt.V>0)
    gradF = zeros(opt.M,opt.V);
    gradG = zeros(opt.C,opt.V);
    
    for k=1:opt.V
        delh1 = opt.delh;%can be different 
        delh2 = opt.delh;%can be different
        if (xv(k) + delh1 > opt.bound(2,k))
            xvCopy11 = xv;
            xvCopy22 = xv ; 
            xvCopy22(k) = xv(k) - delh2;%need to check whether it is lower than the lower bound, not done yet
            delh1 = 0;
        elseif (xv(k) - delh2 < 0.0)
            xvCopy11 = xv;
            xvCopy22 = xv ;
            xvCopy11(k) = xv(k) + delh1;%need to check whether it is higher than the upper bound, not done yet
            delh2 = 0;
        else
            xvCopy11 = xv;
            xvCopy22 = xv ;
            xvCopy11(k) = xv(k) + delh1;
            xvCopy22(k) = xv(k) - delh2;
        end
        
        [perturb1,perturb_cons1] = high_fidelity_evaluation(opt, xvCopy11);
        [perturb2,perturb_cons2] = high_fidelity_evaluation(opt, xvCopy22);
        for kk = 1: opt.M %for all objectives
            gradF(kk,k) = (perturb1(kk)-perturb2(kk))/(delh1+delh2);%delh1 and delh2 are same
        end
        if opt.C>0
            for kk = 1: opt.C %for all objectives
                gradG(kk,k) = (perturb_cons1(kk)-perturb_cons2(kk))/(delh1+delh2);%delh1 and delh2 are same
            end     
        end
    end
end
gradF = gradF';
gradG = gradG';
