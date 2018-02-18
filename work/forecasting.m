function forecasting
%LIST = [1,2,6,7,12,14:17,19,24:25,28:32,36:37,39:42,44,51:52,54,55,58,60:64,66,71,75 ...
%    ,77:79,83:88,91,93:97,99:102,104:109]; % patient numbers 
LIST = 1; 
    
% Control Parameters 
nFitting = 3;                       % Number of cycles used for fitting data
nForecast = 2;                      % Number of cycles of data to forecast
total_n = nFitting + 1 + nForecast; % Total number of treatment cycles 
 
 
% Loads parameters for all patients for all models 
load('par_all_two_pop.mat','par_store')
par2 = par_store; 
load('par_all_portz.mat','par_store')
par3 = par_store;

 
% Remove patients where fitting failed in one of the models 
% y = find(0==sum(par2)); 
% z = find(0==sum(par3));
% LIST(unique([y, z]))=[];      
% par2(:,unique([y, z])) = [];  par3(:,unique([y, z])) = [];
%  
%%
fittingErrorPsa = zeros(length(LIST),4); % error for each patient for each model
forecastErrorPsa = zeros(length(LIST),4); % Forecast error for each patient for each model
fittingErrorAndro = zeros(length(LIST),4); % error for each patient for each model
forecastErrorAndro = zeros(length(LIST),4); % Forecast error for each patient for each model
m = 1; 
index = 1; 
for kk = LIST 
    change = [];
    % Loads Patient data for each patient in list
    patient   = strcat('patient',num2str(kk));       % Patient with corresp number
    file      = strcat('Data/',patient,'.txt');      % Complete name of file patient#.txt
    var       = load(file);                          % holds variable just loaded
    patient   = var;
    t         = patient(:,2);                        %time
    psa       = patient(:,3);                        %psa data
    androgen  = patient(:,4);                        %androgen data
    treatment = patient(:,6);                        %1 is on 0 is off
     
    % Creates change in treatment vector 
    jj = 1;
    change(1) = 1;                        % Treatment starts at t = 0
    for a = 1:length(treatment)
        if treatment(a) ~= mod(jj,2) % When treatment change occurs
            jj = jj + 1;
            change(jj) = a;           % Stores time in change vector
        end
    end
    change(jj+1) = length(treatment);     % Last day of treatment
     
    for i = 1:length(change)-1
        a(i) = change(i+1)-change(i);
    end
    if a > 1 
        b = 1; 
    else
        b = 0; 
    end
    
     
    if (length(change) >= total_n) && b
        m = m + 1; 
               modelCounter = 1;
        modelNames = {'two_pop','portz'};
        for i = 1:2
            modelName = modelNames{i};
            
            % Check model and chooses correct parameters and initial conditions
            switch modelName
                case 'two_pop'
                    params = par2(:,index)';
                    x0 = [androgen(1);99;1;psa(1)]; %collects init cond in a column vector
                case 'portz'
                    params = par3(:,index)';
                    x0 = [99;1;.5;.5;psa(1)];
            end
            
            [errp,erra] =  objective(params,psa,androgen,t,change,x0,nFitting,modelName);
            [errpf,erraf] =  future(params,psa,androgen,t,change,x0,nFitting,nForecast,modelName);
            fittingErrorPsa(kk,modelCounter) = errp;
            forecastErrorPsa(kk,modelCounter) = errpf;
            fittingErrorAndro(kk,modelCounter) = erra;
            forecastErrorAndro(kk,modelCounter) = erraf;
            modelCounter = modelCounter + 1;
             
        end
               
              index = index+1;
     
    end
    
end
  
end
 
function [errp,erra]= objective(params,psadata,and_data,tdata,change,x0,n,model) 
%{
 This function serves as the objective function that fmincon uses to find
 optimal parameter values. 
 params  =  vector of parameters to fit
 psadata =  psa data 
 change  =  vecotor that with the time steps at which treatmet is switched
 tdata   =  time data
 n       =  the number of periods of treatment to fit 
%}
psadata = psadata(change(1):change(n+1));    
and_data = and_data(change(1):change(n+1));
[Y1run] = run_model(params,tdata,change,x0,1,n,model,and_data);
 
switch model 
    case 'one_pop'
        PSA = Y1run(:,3);
        AND = Y1run(:,1);
        errp = sum((PSA-psadata).^2/length(PSA));
        erra = sum((AND-and_data).^2/length(AND));
        err = errp + erra;
         
    case 'two_pop'
        PSA = Y1run(:,4);
        AND = Y1run(:,1);
        errp = sum((PSA-psadata).^2/length(PSA));
        erra = sum((AND-and_data).^2/length(AND));
        err = errp + erra;
         
    case 'portz'
        PSA = Y1run(:,5);
        errp = sum((PSA-psadata).^2/length(PSA));
        erra = 0;
         
    case 'hirata'
        PSA = Y1run(:,1) + Y1run(:,2) + Y1run(:,3);
        errp = sum((PSA-psadata).^2/length(PSA));
        erra = 0;
 
end
end
 
function [errp,erra]= future(params,psadata,and_data,tdata,change,x0,nFitting,nForecast,model) 
%{
 This function serves as the objective function that fmincon uses to find
 optimal parameter values. 
 params  =  vector of parameters to fit
 psadata =  psa data 
 change  =  vecotor that with the time steps at which treatmet is switched
 tdata   =  time data
 n       =  the number of periods of treatment to fit 
%}
psadata = psadata(change(nFitting+1):change(nFitting+1+nForecast));    
and_data = and_data(change(nFitting+1):change(nFitting+1+nForecast));
[Y1run] = run_model(params,tdata,change,x0,nFitting+1,nFitting+nForecast,model,and_data);

 
plot(Y1run(:,3));
 
switch model 
      case 'one_pop'
        PSA = Y1run(:,3);
        AND = Y1run(:,1);
        errp = sum((PSA-psadata).^2/length(PSA));
        erra = sum((AND-and_data).^2/length(AND));
        err = errp + erra;
        fprintf('PSA Error = %.4f \t Androgen Error = %.4f \n',errp, erra);
         
    case 'two_pop'
        PSA = Y1run(:,4);
        AND = Y1run(:,1);
        errp = sum((PSA-psadata).^2/length(PSA));
        erra = sum((AND-and_data).^2/length(AND));
        err = errp + erra;
        fprintf('PSA Error = %.4f \t Androgen Error = %.4f \n',errp, erra);
         
    case 'portz'
        PSA = Y1run(:,5);
        errp = sum((PSA-psadata).^2/length(PSA));
        erra = 0;
        fprintf('PSA Error = %.4f \n',errp);
         
    case 'hirata'
        PSA = Y1run(:,1) + Y1run(:,2) + Y1run(:,3);
        errp = sum((PSA-psadata).^2/length(PSA));
        erra = 0;
        fprintf('PSA Error = %.4f \n',errp);
end
end
 
%%%%%%%%%%%%%%%%  Runs the Model %%%%%%%%%%%%
function [y] = run_model(params,tdata,change,x0,initial_n,final_n,model,androgen)   
y = [];                                      % Initial vector for solution
tint = tdata(change(final_n+1));                   % Time interval to run solution 
for k = initial_n:final_n                                  % K number of cycles of treatment
    u = 1 - mod(k,2);    
    switch model 
        case 'one_pop'
            if tint(end) >= tdata(change(k))  
                [~,Yrun]=ode15s(@(t,x) one_pop(t,x,[params(1:end-1),u]),tdata(change(k):change(k+1)),x0);
                x0 = [Yrun(end,1);Yrun(end,2);Yrun(end,3);Yrun(end,4)];
                if k < final_n
                    y = [y; Yrun(1:end-1,:)]; %#ok<AGROW>
                elseif k == final_n
                    y = [y; Yrun(1:end,:)];   %#ok<AGROW>
                end
            end
        case 'two_pop'
            if tint(end) >= tdata(change(k))
                [~,Yrun]=ode15s(@(t,x) two_pop(t,x,[params(1:end-1),u]),tdata(change(k):change(k+1)),x0);
                x0 = [Yrun(end,1);Yrun(end,2);Yrun(end,3);Yrun(end,4)];
                if k < final_n
                    y = [y; Yrun(1:end-1,:)]; %#ok<AGROW>
                elseif k == final_n
                    y = [y; Yrun(1:end,:)];   %#ok<AGROW>
                end
            end
        case 'portz'
            if tint(end) >= tdata(change(k))
                [~,Yrun]=ode15s(@(t,x) portz(t,x,params,androgen,change(k):change(k+1)),tdata(change(k):change(k+1)),x0);
                x0 = [Yrun(end,1);Yrun(end,2);Yrun(end,3);Yrun(end,4);Yrun(end,5)];
                if k < final_n
                    y = [y; Yrun(1:end-1,:)]; %#ok<AGROW>
                elseif k == final_n
                    y = [y; Yrun(1:end,:)];   %#ok<AGROW>
                end
            end 
        case 'hirata'
            if tint(end) >= tdata(change(k))
                [~,Yrun]=ode15s(@(t,x) hirata(t,x,params,u),tdata(change(k):change(k+1)),x0);
                x0 = [Yrun(end,1);Yrun(end,2);Yrun(end,3)];
                if k < final_n
                    y = [y; Yrun(1:end-1,:)]; %#ok<AGROW>
                elseif k == final_n
                    y = [y; Yrun(1:end,:)];   %#ok<AGROW>
                end
            end
    end   
end
end
% model ODE functions
function dxdt = one_pop(~,x,p) 
% collect parameter values to pass to ODE function
mu = p(1);       q = p(2);      R = p(3);  d = p(4); dd = p(5);
gamma1 = p(6);   gamma2 = p(7); Qm = p(8); b = p(9); sigma = p(10);
epsilon = p(11); u = p(12);
%separates solutions 
Q = x(1); X = x(2); P = x(3);V = x(4);
%%%%%%%%%%%%%%%%%%% ODE system %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dX = mu*(1-q/Q)*X - V*R*X/(Q+R)- abs(dd)*X^2;
dA = (gamma1*u +gamma2)*(Qm -Q)- mu*(Q-q);
dP = abs(b)*Q + abs(sigma)*X*Q - abs(epsilon)*P;
dV = -d*V;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dxdt = [dA; dX; dP; dV]; %Puts the ode in a column vector
end
 
function dxdt = two_pop(~,x,p) 
% collect parameter values to pass to ODE function
um = p(1);        q1= p(2);        q2= p(3);      c1= p(4);         c2= p(5);  K1= p(6);  
K2= p(7);         b= p(8);         sigma1= p(9);  epsilon= p(10);   d1= p(11);  d2= p(12); 
R1= p(13);        R2= p(14);       gamma1= p(15); gamma2= p(16);    dd1 = abs(p(17));
dd2 = abs(p(18)); Qm = abs(p(19)); u= p(20);
%separates solutions 
Q = x(1); X = x(2); Y = x(3); P = x(4);
%%%%%%%%%%%%%%%%%%%  Parameters for ODE System %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if q1>Q
    ux = 0;
else
    ux = abs(um)*(1 - abs(q1)/Q);
end
if q2 > Q
    uy = 0;
else
    uy = abs(um)*(1 - abs(q2)/Q);
end
Dx = abs(d1)*R1/(Q+R1);    Dy = abs(d2)*R2/(Q+R2);
mxy = abs(c1)*K1/(Q + K1); myx = abs(c2)*Q/(Q + K2);
%%%%%%%%%%%%%%%%%%% ODE system %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dX = (ux - Dx  - abs(dd1)*X -  mxy)*X + myx*Y;
dY = (uy - Dy  - abs(dd2)*Y - myx)*Y + mxy*X;
dA = (abs(gamma1)*u +abs(gamma2))*(Qm -Q) - (ux*Q*X + uy*Q*Y)/(X+Y);
dP = b*Q + sigma1*(Y*Q + X*Q) - epsilon*P;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dxdt = [dA; dX; dY;dP]; %Puts the ode in a column vector
end
 
function dxdt = hirata(~,x,p,treat) 
% collect parameter values to pass to ODE function
w110 = p(1);   w11f = p(7);
w21  = p(2);   w12f = p(8);
w220 = p(3);   w22f = p(9);
w31 = p(4);    w33f = p(10);
w32 = p(5);
w33 = p(6);
%separates solutions 
X1 = x(1); X2 = x(2); X3 = x(3);
%%%%%%%%%%%%%%%%%%% ODE system %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if treat == 0
%%%%%%%%%%%%%%%%%%% On Treatment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dX1 = w110*X1;
    dX2 = w21*X1 + w220*X2;
    dX3 = w31*X1 + w32*X2 + w33*X3;
%%%%%%%%%%%%%%%%%% Of Treatment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif treat == 1
    dX1 = w11f*X1;
    dX2 = w12f*X1 + w22f*X2;
    dX3 = w33f*X3;
end
dxdt = [dX1; dX2; dX3]; %Puts the ode in a column vector
end
 
function dxdt = portz(t,x,p,androgen,trange)
% Parameter Values 
um = p(1);     qx = p(2);     qy = p(3);      dx = p(4);     dy = p(5);     c1 = p(6); 
Kxyn = p(7);   c2 = p(8);    Kyxn = p(9);     n = p(10);     qm = p(11);    vm = p(12);    
vh = p(13);    b = p(14);    sigmax = p(15);  sigmay = p(16);rhoxm = p(17); rhoym = p(18);
m = p(19);   sigma0 = p(20); delta  = p(21);
%separates solutions 
X = x(1); Y = x(2); Qx = x(3); Qy = x(4); P = x(5);
% Androgen function constructed from data
A = androgen(end) + (androgen(1)-androgen(end))*exp(-(t-trange(1)));
% Combines parameters to form expressions used in model 
Qxn = Qx^n; Qyn = Qy^n;
Qxm = Qx^m; Qym = Qx^m;
ux = um*(1 - qx/Qx);
uy = um*(1 - qy/Qy);
mxy = c1*Kxyn/(Qxn + Kxyn);
myx = c2*Qyn/(Qyn + Kyxn);
vx = ((qm - Qx)/(qm - qx))*vm*(A/(A + vh));
vy = ((qm - Qy)/(qm - qy))*vm*(A/(A + vh));
% Portz et. al. Model 
dX = (ux - dx - mxy)*X + myx*Y;
dY = (uy - dy - myx)*Y + mxy*X;
dQx = vx - ux*Qx - b*Qx;
dQy = vy - uy*Qy - b*Qy;
dP = sigmax*X*Qxm/(Qxm + rhoxm) + sigmay*Y*Qym/(Qym + rhoym) + sigma0*(X + Y) - delta*P;
%Puts the ode in a column vector
dxdt = [dX; dY; dQx; dQy; dP]; 
end