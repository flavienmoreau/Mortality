function LE=life_exp4(MR,uptoage)
% return the life expectancy associated with a mortality rate vector
% input: X is ln per 100 mortality rate
% EDIT June 2018 to make it size consistent

    % Size Compatibility: ensure vector is vertical
    [n,m] = size(MR);
    if n==1
        MR = MR';
    end   

    % Handle NAN and terminal value
    nans = find(isnan(MR));
    if isempty(nans)
        T_nan = 101;
    else
        T_nan = min(nans);
    end
    
    surv=survivalcurve(MR);    
    T0 = length(surv);
    
    T_MR=min(T_nan-1,T0);

    
    try                         % if no last age given then use 100 by default
        T_MR=min(uptoage,T_MR);
    catch
        T_MR=min(100,T_MR);
    end 
    
    surv=surv(1:T_MR);            % censor at this age
    MR=MR(1:T_MR);
    
    age=1:T_MR;
    age=age-0.5*ones(1,T_MR);       % o/w bias as people die during the year
    
    LE=age*(surv.*(10.^(MR)/100)) + (surv(T_MR)*(T_MR+1))/100;     % censor after max observed age
    
end