function [surv]=survivalcurve(X)
% SEPT2017: 10.^ instead of exp, as MR now in log10
% takes a log10 per 100 mortality rate vector 
% and returns a survival curve of same size + 1 
    TT=length(X);
    XX=10.^(X);
    surv=100*ones(TT,1);
    
    surv(1)=100;
    for t=2:TT
        surv(t)=surv(t-1)*(1-XX(t-1));
    end
end
