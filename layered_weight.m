function w =  layered_weight(alpha, H)

        [n,m] = size(H);
        a = (1-alpha)*(ones(1, m)*(1/m));
        w = a + alpha*H;


end