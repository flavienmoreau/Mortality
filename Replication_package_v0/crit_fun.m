function Q=crit_fun(X,W)
% FM OCT2017
% used to save space, can be improved with safeguards

    % Protection against Inf Feb 6 2018
    f = isfinite(X);
    
    Q = X(f)'*diag(W(f))*X(f);

end
