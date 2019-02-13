function Pnew = shiftRefSet(Pold, popBest, M, p)

if ~isnan(Pold)
    sumPold = sum(Pold(1,:));
    sumPnew = sum(popBest);
    Pnew = Pold * (sumPnew/sumPold);
else
    P = initweight( M, nchoosek(M+p-1, p) );
    if p == 2 || p == 3
        p = p-1;
        Padd = 0.5 * initWeight( M, nchoosek(M+p-1, p) ); 
        n = size(Padd,1);
        Pnew = zeros(size(Padd));
        for i=1:n
            t = (1 - sum(Padd(i,:))) / M;
            Pnew(i,:) = t + Padd(i,:);
        end
        P = [P; Pnew];
    end
    Pnew = P * sum(popBest);
end

