function [logMR,H,Halive]=logMRadd(X)

% Generates log mortality rates given parameters  
%  syntax: X =  I,delta,sigma_e,alpha,mu_H,sigma_H,Hbar

    % use same seed across all simulations
    global seed
    rng(seed)

    %Parameters for simulation
    N= 200000; % # of individuals simulated, reduce for improved speed
    %global T
    T=100;      % 100 years or pass it globally
    I=X(1);
    delta=X(2);
    sigma_e=X(3);
    alpha=X(4);
    mu_H=X(5);
    sigma_H=X(6);
    Hbar=X(7);

    % Initialization
    H=zeros(N,T);
    H(:,1)=mu_H+sigma_H.*randn(N,1);
    deathyear=T*ones(N,1); %will record death year
    alive=ones(N,T); %record survivors
    Nt=N*ones(T,1); %count surviving population
    MR=ones(1,T);

    % Iteration
    justdied=find(H(:,1)<Hbar); %find those who just died
    mortality=length(justdied); %count the deaths
    Nt=N-mortality; 
    deathyear(justdied)=1; %record infant mortality rate
    alive(justdied,1)=0;
    MR(1,1)=mortality/N; % here maybe 

    t=2;
    go_on=1;
    while t<=T & go_on==1
        Eps=sigma_e*randn(N,1); %draw shock
        
        H(:,t)=H(:,t-1)+(I-delta*(t^alpha)+Eps);
        alive(:,t)=alive(:,t-1)==1 & H(:,t)>Hbar;
        justdied=find(alive(:,t-1)~=alive(:,t));
        mortality=length(justdied);
        Nt=Nt-mortality;
        deathyear(justdied)=t; %record death year
        
        if sum(alive(:,t))>0
            MR(1,t)=mortality/sum(alive(:,t));
        else
            go_on=0;
        end
        t=t+1;
    end

    logMR=log10(MR');
    
    %safety against NaN 
    %log_NaN=isnan(logMR);       % return indices for which logMR is NaN
    %logMR(log_NaN)=10^(-9);
    
    Halive=H.*alive;
    %CHECKS AND FIGURES
    % mortality                     OK
    %{
    plot(log(MR));
    title('log mortality, 1m indivdiuals');
    %}

end
