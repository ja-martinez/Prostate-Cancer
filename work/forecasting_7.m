function forecasting_7
%{
    This function saves the fitting and forecasting errors for all models
    and saves the plots forecasting plots. There is one plot per (usable) patient
    and the plots include all three model forecasting.

    To save errors, creat folders 'two_pop', 'three_pop', and 'portz'.
    TO save plots,  create folder 'Plots'
    Choose patients in ~line 79
%}

%%%%%%%%%%%%%%%%%%%%%%%% Model Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose whether you want to fit all patients (q=1) or not (q=0)
q = 0;
if q == 1
    LIST = [1,2,6,7,12,14:17,19,24:25,28,29,31,32,36:37,39:40,42,44,51:52,54,55,58,60:64,66,75 ...
        ,77:79,83:85,87,88,91,93:97,99:102,104:109]; % patient numbers 
elseif q == 0
    LIST = [1]; % patient numbers 
else
    warning('Unexpected q value, choose a different one')
    return;
end

% Control Parameters 
nFitting = 3;                       % Number of cycles used for fitting data
nForecast = 2;                      % Number of cycles of data to forecast
total_n = nFitting + 1 + nForecast; % Total number of treatment cycles 

fittingErrorPsa_two_pop = zeros(length(LIST));  
forecastErrorPsa_two_pop = zeros(length(LIST)); 
fittingErrorAnd_two_pop = zeros(length(LIST)); 
forecastErrorAnd_two_pop = zeros(length(LIST)); 

fittingErrorPsa_three_pop = zeros(length(LIST));  
forecastErrorPsa_three_pop = zeros(length(LIST)); 
fittingErrorAnd_three_pop = zeros(length(LIST)); 
forecastErrorAnd_three_pop = zeros(length(LIST)); 

fittingErrorPsa_portz = zeros(length(LIST));  
forecastErrorPsa_portz = zeros(length(LIST)); 

fpath = 'C:\Users\Alejandro Martinez\Documents\Mathematics\ProstateCancer\ModifiedCode\Plots'; % Plot destination

% Run through models
for yy = 1;
    if yy == 1
        model = 'two_pop';
    elseif yy == 2
        model = 'three_pop';
    elseif yy == 3
        model = 'portz';
    end

    switch model    % Load parameters    
        case 'two_pop'
            load('par_all_two_pop.mat','par_store')
             par = par_store;
        case 'three_pop'
            load('par_all_three_pop.mat','par_store')
            par = par_store;
        case 'portz' 
            load('par_all_portz.mat','par_store')
            par = par_store;
        otherwise 
            warning('Unexpected model, choose a different one')
            return;
    end

    % % Remove patients where fitting failed in one of the models 
    % y = find(0==sum(par)); 
    % LIST(unique(y))=[];      
    % par(:,unique(y)) = [];
    % length(par)
    
    patientIndex = 1; % initialize index to store patients who are forcasted
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Run through all patients   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for kk = LIST
        %change = zeros(1, total_n);                      % Initialize change vector that saves time array index of when treatment changes
        % Loads Patient data for each patient in list
        file      = strcat('Data_7\case',num2str(kk),'_data.mat');      % Complete name of file patient#.txt
        load(file);                          % holds variable just loaded
        patient   = data;
        t         = patient(:,1);                        % time
        psa       = patient(:,3);                        % psa data
        androgen  = patient(:,2);                        % androgen data
        treatment = patient(:,4);                        % 1 is on 0 is off
        index = find(LIST==kk);
        
        
        % Assign values to change vector
        jj = 1;
        change(1) = 1;                        % Treatment starts at t = 0
        for a = 1:length(treatment)
            if treatment(a) ~= mod(jj,2) % When treatment change occurs
                jj = jj + 1;
                change(jj) = a;           % Stores time in change vector
            end
        end
        change(jj+1) = length(treatment);     % Last day of treatment
        
        % Determine if there is only one data point for each treatment cycle. b=1 if more than one, b=0 if not
        b = 1;                  % initialize b value
        for i = 1:length(change)-1
            a = change(i+1)-change(i);       % a(i) tells you lngth of the ith treatment
            if a == 1
                b = 0;
                break;
            end
        end   
        % Do forecasting if there are enough cycles for forecasting and if each cycle has more than one data point
        if (length(change) >= total_n) && b           
            patientNumber(patientIndex) = kk;       % Saves the patients who are forcasted
            patientIndex = patientIndex + 1;
            %Part of plotting
            figure(kk);
            title(strcat('patient',num2str(kk)), 'fontsize',26);
            hold on;

            switch model
                case 'two_pop'
                    params = par(:,index)';
                    x0 = [androgen(1);par(20,index);par(21,index);psa(1)]; %collects init cond in a column vector
                
                    [errors,x0_new] =  objective(params,psa,androgen,t,change,x0,nFitting,model);     % get new initial conditions, fitting errors, and plot fitting
                    errp = errors(1);
                    erra = errors(2);
                    
                    [errorsf] = future(params,psa,androgen,t,change,x0_new,nFitting,nForecast,model);     % get forecasting errors and plot the forecasting
                    errpf = errorsf(1);
                    erraf = errorsf(2);

                    fittingErrorPsa_two_pop(index) = errp;
                    forecastErrorPsa_two_pop(index) = errpf;
                    fittingErrorAnd_two_pop(index) = erra;
                    forecastErrorAnd_two_pop(index) = erraf;    

                case 'three_pop'
                    params = par(:,index)';
                    x0 = [androgen(1);par(22,index);par(23,index);par(24,index);psa(1)]; %collects init cond in a column vector
                    
                    [errors,x0_new] =  objective(params,psa,androgen,t,change,x0,nFitting,model);     % get new initial conditions, fitting errors, and plot fitting
                    errp = errors(1);
                    erra = errors(2);
                    
                    [errorsf] = future(params,psa,androgen,t,change,x0_new,nFitting,nForecast,model);     % get forecasting errors and plot the forecasting
                    errpf = errorsf(1);
                    erraf = errorsf(2);

                    fittingErrorPsa_three_pop(index) = errp;
                    forecastErrorPsa_three_pop(index) = errpf;
                    fittingErrorAnd_three_pop(index) = erra;
                    forecastErrorAnd_three_pop(index) = erraf;     

                case 'portz'
                    params = par(:, index)';
                    x0 = [par(22,index);par(23,index);par(24,index);par(25,index);psa(1)];

                    [errors,x0_new] =  objective(params,psa,androgen,t,change,x0,nFitting,model);     % get new initial conditions, fitting errors, and plot fitting
                    
                    [errorsf] = future(params,psa,androgen,t,change,x0_new,nFitting,nForecast,model);     % get forecasting errors and plot the forecasting                    

                    fittingErrorPsa_portz(index) = errors;
                    forecastErrorPsa_portz(index) = errorsf;
                    
                           
            end
            %Part of plotting
            scatter(t(1:change(total_n)),psa(1:change(total_n)),'k','linewidth',2);        % plot trial data
            hold on;
            line([t(change(nFitting+1)),t(change(nFitting+1))],ylim);         % plot line that separates the fitting and forecasting sections of the plot
            xlabel('Time (Days)', 'fontsize', 26);
            ylabel('PSA (\mug/L)', 'fontsize', 26);
            p1 = plot(nan,nan,'m');
            p2 = plot(nan,nan,'g');
            p3 = plot(nan,nan,'b');
            legend([p1 p2 p3],'Model 2','Model 3','PKN')
            hold off;
        end
           

        
    end
    
    for w = patientIndex
        figure(w)
        fileNameAndrogen = strcat('Patient',num2str(kk),'Androgen');   % Create string for androgen plot file name
        saveas(gcf, fullfile(fpath, fileNameAndrogen), 'jpeg');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Save Error Arrays %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch model
    
        case 'two_pop'
            save('two_pop/fittingErrorPsa_two_pop.mat','fittingErrorPsa_two_pop') 
            save('two_pop/fittingErrorAnd_two_pop.mat','fittingErrorAnd_two_pop')
            save('two_pop/forecastErrorPsa_two_pop.mat','forecastErrorPsa_two_pop')
            save('two_pop/forecastErrorAnd_two_pop.mat','forecastErrorAnd_two_pop')
        case 'three_pop'
            save('three_pop/fittingErrorPsa_three_pop.mat','fittingErrorPsa_three_pop') 
            save('three_pop/fittingErrorAnd_three_pop.mat','fittingErrorAnd_three_pop')
            save('three_pop/forecastErrorPsa_three_pop.mat','forecastErrorPsa_three_pop')
            save('three_pop/forecastErrorAnd_three_pop.mat','forecastErrorAnd_three_pop')
        case 'portz'
            save('portz/fittingErrorPsa_portz.mat','fittingErrorPsa_portz') 
            save('portz/forecastErrorPsa_portz.mat','forecastErrorPsa_portz')
    end
end
end

function [errors,x0_new] = objective(params,psadata,and_data,tdata,change,x0,n,model) 
%{
 This function takes the parameters from the fitting and the initial
 conditions from the trial data.

 This function plots the fitting model data and returns the fitting 
 errors and the values of the dependent variables at the end of the
 fitting cycle.
%}

psadata = psadata(change(1):change(n+1));    
and_data = and_data(change(1):change(n+1));
[Y1run] = run_model(params,tdata,change,x0,1,n,model,and_data);
    
switch model
    case 'two_pop'
        PSA = Y1run(:,4);
        AND = Y1run(:,1);
        errp = sum((PSA-psadata).^2/length(PSA));
        erra = sum((AND-and_data).^2/length(AND));
        errors = [errp, erra];
        
        %Update initial conditions
        X10_new = Y1run(end,2);
        X20_new = Y1run(end,3);
        PSA0_new = Y1run(end,4);
        Q0_new = Y1run(end,1);
        
        x0_new = [Q0_new, X10_new, X20_new, PSA0_new];
        
        plot(tdata(1:length(PSA)),PSA,'m','linewidth',2);
        hold on;

    case 'three_pop'
        PSA = Y1run(:,5);
        AND = Y1run(:,1);
        errp = sum((PSA-psadata).^2/length(PSA));
        erra = sum((AND-and_data).^2/length(AND));
        errors = [errp, erra];

        %Update initial conditions
        X10_new = Y1run(end,2);
        X20_new = Y1run(end,3);
        X30_new = Y1run(end,4);
        PSA0_new = Y1run(end,5);
        Q0_new = Y1run(end,1);
        x0_new = [Q0_new, X10_new, X20_new, X30_new, PSA0_new];
        
        plot(tdata(1:length(PSA)),PSA,'g','linewidth',2);
        hold on;

    case 'portz'
        PSA = Y1run(:,5);
        err = sum((PSA-psadata).^2/length(PSA));
        errors = err;

        %Update initial conditions
        X10_new = Y1run(end,1);
        X20_new = Y1run(end,2);
        Q10_new = Y1run(end,3);
        Q20_new = Y1run(end,4);
        PSA0_new = Y1run(end,5);
        x0_new = [X10_new, X20_new, Q10_new, Q20_new, PSA0_new];
        
        plot(tdata(1:length(PSA)),PSA,'b','linewidth',2);
        hold on;
end


    
end

function [errors]= future(params,psadata,and_data,tdata,change,x0,nFitting,nForecast,model) 
%{
 This function takes the parameters from the fitting and the new initial conditions
 from where the fitting model left off to run the model again.

 It adds a forecasting plot and returns the forecasting errors.
%}
psadata = psadata(change(nFitting+1):change(nFitting+1+nForecast));                % retrieve psa data from after fitting cycles
androgen = and_data;    % Androgen used for portz
and_data = and_data(change(nFitting+1):change(nFitting+1+nForecast));       %androgen used for errors                          % retrieve androgen data from after the fitting cycles
[Y1run] = run_model(params,tdata,change,x0,nFitting+1,nFitting+nForecast,model,androgen);             % run model with the initial conditions set to where the fitting left off

switch model
    case 'two_pop'
        PSA = Y1run(:,4);
        AND = Y1run(:,1);        
        errp = sum((PSA-psadata).^2/length(PSA));                        % forecasting psa error
        erra = sum((AND-and_data).^2/length(AND));                       % forecasting androgen error
        errors = [errp, erra];
        
        plot(tdata(change(nFitting+1):change(nFitting+1+nForecast)),PSA,'m','linewidth',2);            % plot the forecasting model data
        hold on;

    case 'three_pop'
        PSA = Y1run(:,5);
        AND = Y1run(:,1);
        errp = sum((PSA-psadata).^2/length(PSA));
        erra = sum((AND-and_data).^2/length(AND));
        errors = [errp, erra];
        
        plot(tdata(change(nFitting+1):change(nFitting+1+nForecast)),PSA,'g','linewidth',2);            % plot the forecasting model data
        hold on;
    case 'portz'
        PSA = Y1run(:,5);
        err = sum((PSA-psadata).^2/length(PSA));
        errors = err;
        plot(tdata(change(nFitting+1):change(nFitting+1+nForecast)),PSA,'b','linewidth',2);            % plot the forecasting model data
        hold on;
end

end

%%%%%%%%%%%%%%%%  Runs the Model %%%%%%%%%%%%

function [y] = run_model(params,tdata,change,x0,initial_n, final_n, model,androgen)   
y = [];                                      % Initial vector for solution
tint = tdata(change(final_n+1));             % Time interval to run solution 
for k = initial_n:final_n                    % K number of cycles of treatment
    u = 1 - mod(k,2);    
    switch model             
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
%            
        case 'three_pop'
            if tint(end) >= tdata(change(k))
                [~,Yrun]=ode15s(@(t,x) three_pop(t,x,[params(1:end-1),u]),tdata(change(k):change(k+1)),x0);
                x0 = [Yrun(end,1);Yrun(end,2);Yrun(end,3);Yrun(end,4);Yrun(end,5)];
                if k < final_n
                    y = [y; Yrun(1:end-1,:)]; %#ok<AGROW>
                elseif k == final_n
                    y = [y; Yrun(1:end,:)];   %#ok<AGROW>
                end
            end
%            
        case 'portz'
            if tint(end) >= tdata(change(k))
%                options = odeset('AbsTol',1e-14,'RelTol',1e-14);
                [~,Yrun]=ode23tb(@(t,x) portz(t,x,params,androgen,change(k):change(k+1)),tdata(change(k):change(k+1)),x0);
                x0 = [Yrun(end,1);Yrun(end,2);Yrun(end,3);Yrun(end,4);Yrun(end,5)];
                if k < final_n
                    y = [y; Yrun(1:end-1,:)]; %#ok<AGROW>
                elseif k == final_n
                    y = [y; Yrun(1:end,:)];   %#ok<AGROW>
                end
            end 
%            
    end   
end
end

%% Model ODE functions

% Two population model
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


% Three population Model
function dxdt = three_pop(~,x,p) 
% collect parameter values to pass to ODE function
um = p(1);         q1= p(2);          q2= p(3);      q3 = p(4);
c= p(5);           K= p(6);           b= p(7);       sigma= p(8);      
epsilon= p(9);     d1= p(10);         d2= p(11);     d3= p(12);
R1= p(13);         R2= p(14);         R3 = p(15);
gamma1= p(16);     gamma2= p(17);
dd1 = abs(p(18));  dd2 = abs(p(19));  dd3 = abs(p(20)); 
Qm = abs(p(21));   u=p(25);
%separates solutions 
Q = x(1); X = x(2); Y = x(3); Z = x(4); P = x(5);
%%%%%%%%%%%%%%%%%%%  Parameters for ODE System %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if q1>Q
    ux = 0;
    else
    ux = abs(um)*(1 - abs(q1)/Q);
    end
%
    if q2 > Q
    uy = 0;
    else
    uy = abs(um)*(1 - abs(q2)/Q);
    end
%
    if q3 > Q
    uz = 0;
    else
    uz = abs(um)*(1 - abs(q3)/Q);
    end
%
Dx = abs(d1)*R1/(Q+R1);    
Dy = abs(d2)*R2/(Q+R2);
Dz = abs(d3)*R3/(Q+R3);
Lm = abs(c)*K/(Q + K); % "Mutation" rate
gamma = abs(gamma1)*u +abs(gamma2);
%%%%%%%%%%%%%%%%%%% ODE system %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dX = (ux - Dx  - abs(dd1)*X)*X - 2*Lm*X;
dY = (uy - Dy  - abs(dd2)*Y)*Y + Lm*X - 2*Lm*Y ;
dZ = (uz - Dz - abs(dd3)*Z)*Z + Lm*(X+Y) ;
dA = gamma*(Qm -Q) - (ux*Q*X + uy*Q*Y +uz*Q*Z)/(X+Y+Z);
dP = b*Q + sigma*(Y*Q + X*Q + Z*Q) - epsilon*P;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dxdt = [dA; dX; dY; dZ; dP]; %Puts the ode in a column vector
end

% PKN model
function dxdt = portz(t,x,p,androgen,trange)
% Parameter Values 
um = p(1);     qx = p(2);    qy = p(3);       dx = p(4);     dy = p(5);     c1 = p(6); 
Kxyn = p(7);   c2 = p(8);    Kyxn = p(9);     n = p(10);     qm = p(11);    vm = p(12);    
vh = p(13);    b = p(14);    sigmax = p(15);  sigmay = p(16);rhoxm = p(17); rhoym = p(18);
m = p(19);   sigma0 = p(20); delta  = p(21);
%separates solutions 
X = x(1); Y = x(2); Qx = x(3); Qy = x(4); P = x(5);
% Androgen function constructed from data
A = androgen(trange(end)) + (androgen(trange(1))-androgen(trange(end)))*exp(-(t-trange(1)));
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