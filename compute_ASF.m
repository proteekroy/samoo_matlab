function ASF = compute_ASF(ASF, Yref, row, N, normF)

D =  (normF + repmat(0.01*ones(size(Yref)), N, 1)) ./ (repmat(Yref + 0.01*ones(size(Yref)), N, 1));
maxD = max(D, [], 2);
ASF(row,:) = maxD';

end