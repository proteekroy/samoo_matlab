function ynew = quadratic_predictor(reg, X)

    %reg = opt.reg{i};
    ynew = zeros(size(X,1),1);
    for j=1:size(X,1)
        xnew = X(j,:);
        NewScores=repmat(xnew,[length(reg.PowerMatrix) 1]).^reg.PowerMatrix;
        EvalScores=ones(length(reg.PowerMatrix),1);
        for ii=1:size(reg.PowerMatrix,2)
            EvalScores=EvalScores.*NewScores(:,ii);
        end

        ynew(j) = reg.Coefficients'*EvalScores; % The estimate for the new data point.
    end
end