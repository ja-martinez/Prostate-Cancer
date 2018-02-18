function tin_forecasting
% LIST = [1,2,6,7,12,14:17,19,24:25,28:32,36:37,39:42,44,51:52,54,55,58,60:64,66,71,75 ...
%     ,77:79,83:88,91,93:97,99:102,104:109]; % patient numbers 
LIST = 17; % pick one patient

% Control Parameters 
nFitting = 3;                       % Number of cycles used for fitting data
nForecast = 2;                      % Number of cycles of data to forecast
total_n = nFitting + 1 + nForecast; % Total number of treatment cycles 


% Loads parameters for all patients for two_pop model 
load('par_all_two_pop.mat','par_store')
par2 = par_store; 

% % Remove patients where fitting failed in one of the models 
% y = find(0==sum(par2)); 
% LIST(unique(y))=[];      
% par2(:,unique(y)) = [];
% length(par2)

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
        disp(m)
        modelCounter = 1;
        modelName = 'two_pop';
            params = par2(:,index)';
            x0 = [androgen(1);par2(20,index);par2(21,index);psa(1)]; %collects init cond in a column vector
          
            %Part of plotting
            figure(1)
            hold on;
            
            [errp,erra,PSA0_new,X10_new,X20_new,Q0_new] =  objective(params,psa,androgen,t,change,x0,nFitting,modelName);
            
            %Change initial condition
            x0 = [PSA0_new;X10_new;X20_new;Q0_new];
            
            [errpf,erraf] = future(params,psa,androgen,t,change,x0,nFitting,nForecast,modelName);
            fittingErrorPsa(kk,modelCounter) = errp;
            forecastErrorPsa(kk,modelCounter) = errpf;
            fittingErrorAndro(kk,modelCounter) = erra;
            forecastErrorAndro(kk,modelCounter) = erraf;
            
           %Part of plotting
            scatter(t(1:change(6)),psa(1:change(6)));
            line([t(change(4)),t(change(4))],ylim);
           
        end
              
              index = index+1;
    
    end
           
        
end
 


function [errp,erra,PSA0_new,X10_new,X20_new,Q0_new]= objective(params,psadata,and_data,tdata,change,x0,n,model) 
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
    
        %'two_pop' 
        PSA = Y1run(:,4);
        AND = Y1run(:,1);

        % The plotting is done in this supporting function (1) TWO PARTS.
        plot(tdata(1:length(PSA)),PSA,'b');
        legend('fitting');
        
        errp = sum((PSA-psadata).^2/length(PSA));
        erra = sum((AND-and_data).^2/length(AND));
        err = errp + erra;
        fprintf('PSA Error = %.4f \t Androgen Error = %.4f \n',errp, erra);
        
        %Update X10_new, X20_new.
        X10_new = Y1run(end,2);
        X20_new = Y1run(end,3);
        PSA0_new = Y1run(end,4);
        Q0_new = Y1run(end,1);
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
size(Y1run)
size(psadata)

% plot(Y1run(:,3))
        %'two_pop' 
        PSA = Y1run(:,4);
        AND = Y1run(:,1);
        
        % The plotting is done in this supporting function (2) TWO PARTS
        plot(tdata(change(nFitting+1):change(nFitting+1+nForecast)),PSA,'k');
        legend('forecasting');
        
        errp = sum((PSA-psadata).^2/length(PSA));
        erra = sum((AND-and_data).^2/length(AND));
        err = errp + erra;
        fprintf('PSA Error = %.4f \t Androgen Error = %.4f \n',errp, erra);
end

%%%%%%%%%%%%%%%%  Runs the Model %%%%%%%%%%%%

function [y] = run_model(params,tdata,change,x0,initial_n,final_n,model,androgen)   
y = [];                                      % Initial vector for solution
tint = tdata(change(final_n+1));                   % Time interval to run solution 
for k = initial_n:final_n                                  % K number of cycles of treatment
    u = 1 - mod(k,2);    

            %two_pop
            if tint(end) >= tdata(change(k))
                [~,Yrun]=ode15s(@(t,x) two_pop(t,x,[params(1:end-1),u]),tdata(change(k):change(k+1)),x0);
                x0 = [Yrun(end,1);Yrun(end,2);Yrun(end,3);Yrun(end,4)];
                if k < final_n
                    y = [y; Yrun(1:end-1,:)]; %#ok<AGROW>
                elseif k == final_n
                    y = [y; Yrun(1:end,:)];   %#ok<AGROW>
                end
            end

end
end

% model ODE functions

function dxdt = two_pop(~,x,p) 
% collect parameter values to pass to ODE function
um = p(1);        q1= p(2);        q2= p(3);      c1= p(4);         c2= p(5);  K1= p(6);  
K2= p(7);         b= p(8);         sigma1= p(9);  epsilon= p(10);   d1= p(11);  d2= p(12); 
R1= p(13);        R2= p(14);       gamma1= p(15); gamma2= p(16);    dd1 = abs(p(17));
dd2 = abs(p(18)); Qm = abs(p(19)); u= p(22);
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


