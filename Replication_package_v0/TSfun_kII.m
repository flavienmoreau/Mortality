function [res,fit_S,LMR1,S1,LE1,fit0]=TSfun_kII(year,init_guess)

    global mult g_Eps0 g_EpsN g_EpsU maxit prec_pow
    
    disp('')
    titlestring='';

    % Load Mortality q-rates [Need to be in the folder]
    filename=strcat(strcat('MR_FRA',num2str(year),'F_q.csv')) 

    target = csvread(filename,1,0);
    target = log10(target/100);             % log per 100 rate   
    LMR0 = target;
    Tdata =min(100,length(target));
    S0=survivalcurve(LMR0);              % survival curve
    T_S0=min(100,length(S0));

    % Get age at beginning and end of WW2 
    t1=1940-year
    t2=1946-year

    % Flu pandemic
    tww1=1919-year % really it is 1918-year+1

    % Optimization options
    nvars=8;                % number of parameters to optimize upon

    nb_loops = 20 * mult;           % depth of research

    % Use these options for Hoffman [no figures]
    options_fmin=optimset('Display','off','FunValCheck','on','MaxIter',nb_loops*nvars);
    
    % For local run, with Figures showing optimization progress
    %options_fmin=optimset('Display','iter','PlotFcns',@optimplotfval,'FunValCheck','on','MaxIter',nb_loops*nvars);

    weight=ones(Tdata,1); weight(1)=1;


    % Estimation with fminsearch
    %----------------------------

    % Order:  I , delta, sigma, alpha , mu , WWI (I) , Accident (k), WW2(I)

    fmin_init = init_guess
    
    criterion_s = @(Y) crit_fun(trunc_data(S0,T_S0)-trunc_data(survivalcurve(logMR_IIk([Y(1:5),1,0,Y(6),tww1,tww1,Y(8),t1,t2,Y(7)],g_Eps0,g_EpsN,g_EpsU)),T_S0),weight);

    criterion_bo= @(a) criterion_s(a')

    [res,fit_S]=fminsearch(criterion_s,fmin_init,options_fmin);

    [res,fit]=powell(criterion_bo,res',maxit,prec_pow);
    [res,fit]=fminsearch(criterion_s,res',options_fmin);

    Y=res

    LMR1=logMR_IIk([Y(1:5),1,0,Y(6),tww1,tww1,Y(8),t1,t2,Y(7)],g_Eps0,g_EpsN,g_EpsU);      % Note: function #3 has nothing to do with column #
    %LE_cf=life_exp4(logMRadd_bumps5([Y(1:5),1,0,0,t1,t2,Y(7),14]),Tdata);      % 3. Counterfactual  Life expectancy 


    S1=survivalcurve(LMR1);
    LE1=life_exp4(LMR1,Tdata);                          % 2. Estimate Life expectancy            
    fit_S = criterion_s(res);          % what's the fit at the inital value, to Check if making progress

    % Set anterior shocks to 0:
    if t2<0
        res(8) = 0;
    end
    if tww1 <0
        res(6) = 0;
    end
          
    
    fit0 = criterion_s(fmin_init);          % what's the fit at the inital value, to Check if making progress

%         figure()
%         subplot(1,2,1)
%         plot(1:Tdata,LMR0(1:Tdata),'-x',1:Tdata,LMR1(1:Tdata))
%         legend(num2str(year),'simulation','Location','SouthEast')
%         subplot(1,2,2)
%         plot(1:T_S0,S0(1:T_S0),1:T_S0,S1(1:T_S0))
%        title(titlestring)
end
    