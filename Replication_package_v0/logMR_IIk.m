function logMR=logMR_IIk(X,g_Eps0,g_EpsN,g_EpsU)

% - WW1 : I shock between age tp1 and tp2 
% - I shock for second war period allowed, with tp3 and tp4 
% - Accident shock, lifelong or starting at 14

% Additive specification (instead of log additive)
%  syntax: X =  I,delta,sigma_e,alpha,mu_H,sigma_H,Hbar 
%               I1, tp1, tp2, I2, tp3, tp4, kappa, m1 (start of acc. shck)

    %Parameters for simulation
    %N= 100000; % 1m individuals % reduced to 100k for efficiency July29
    global N T ;
    
    I=X(1);
    delta=X(2);
    sigma_e=X(3);
    alpha=X(4);
    mu_H=X(5);
    sigma_H=X(6);
    Hbar=X(7);
    
    % turning points and exogenous killing for WW1
    I1=X(8);
    
    tp1=X(9);
    tp2=X(10);
      
    %Second war period
    I2=X(11);
    tp3=X(12);
    tp4=X(13);

    %Accident rate
    kappa=X(14);                % fraction of exogenous deaths every year, starting at m1
  
    try                         % if no starting point given assume it is throughout
        m1=X(15);
    catch
        m1=14;
    end


    % Initialization
    H=zeros(N,T);
    H(:,1)=mu_H+sigma_H.*g_Eps0;
    deathyear=T*ones(N,1);              %will record death year
    alive=ones(N,T);                    %record survivors
    Nt=N*ones(T,1);                     %count surviving population
    MR=zeros(1,T);

    % Iteration
    justdied=find(H(:,1)<Hbar);         %find those who just died
    mortality=length(justdied);         %count the deaths
    Nt=N-mortality; 
    deathyear(justdied)=1;              %record infant mortality rate
    alive(justdied,1)=0;
    MR(1,1)=mortality/N;                

    t=2;
    go_on=1;
    while t<T && go_on==1
        Eps=sigma_e*g_EpsN(:,t);         % draw shock        
        exodeath=(g_EpsU(:,t)<kappa) & t>=m1;
        
        if t>=tp1 && t<=tp2
            I=I1;                       % WWI: Investment shock
        elseif t>=tp3 && t<=tp4
            I=I2;                       % WW2: I shock
        else
            I = X(1);
        end
        
        H(:,t)=H(:,t-1)+(I-delta*(t^alpha)+Eps);
        alive(:,t)=alive(:,t-1)==1 & H(:,t)>Hbar & exodeath==0;      % take exo death into account here
        justdied=find(alive(:,t-1)~=alive(:,t));
        mortality=length(justdied);
        Nt=Nt-mortality;
        deathyear(justdied)=t;          % record death year
        
        if sum(alive(:,t))>0
            MR(1,t)=mortality/sum(alive(:,t));
        else
            go_on=0;
        end
        t=t+1;
    end

    logMR=log10(MR');
    
    % safety against NaNs
    log_NaN=isnan(logMR);
    logMR(log_NaN)=10^(-9);
    
    %CHECKS AND FIGURES
    % mortality                     OK
    
    %{
    figure()
    plot(log(MR));
    title('log mortality, 1m indivdiuals');
    %}

end
