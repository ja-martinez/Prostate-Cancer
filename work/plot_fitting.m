function [androgenError, psaError] =  plot_fitting
%{
This files takes the patients numbers from variable LIST and performs parameter 
estimation using fmincon. Parameters are stored for each patient in p#.mat, and for 
every patient in par_all.mat.

Parameter estimation can be done for two_pop,three_pop, portz,and hirata models. 

Created by Javier Baez Sep 2016.
Modified by Tin May 2017.

%} 

n = 6;                                       % number of periods to fit 
    
%%%%%%%%%%%%%%%%%%%%%%%% Initialization of arrays and variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose whether you want to fit all patients (q=1) or not (q=2)
q = 2;
if q == 1
    LIST = [1,2,6,7,12,14:17,19,24:25,28,29,31,32,36:37,39:40,42,44,51:52,54,55,58,60:64,66,75 ...
        ,77:79,83:85,87,88,91,93:97,99:102,104:109]; % patient numbers 
elseif q == 2
    %LIST = [30, 41, 71, 86]; % patient numbers 
    LIST = [15];
else
    warning('Unexpected q value, choose a different one')
    return;
end

totalPatients = length(LIST);              % total number of patients  
patients = cell(1,totalPatients);          % creates a cell array to hold patient data. 
androgenError = zeros([totalPatients, 2]); % creates array to output androgen error
psaError = zeros([totalPatients, 3]);      % creates array to output androgen error
counter = 1;                               % counter for patients and change array
ind = cell(1,totalPatients);               % initializes the index for the S structure
change = ones(1,n+1);                      % vector stores times when treatment changed

%%%%%%% arrays to save data for model   %%%%%%
% Model data
model_2_and = ones([82,totalPatients]);
model_2_psa = ones([82,totalPatients]);
model_3_and = ones([82,totalPatients]);
model_3_psa = ones([82,totalPatients]);
model_portz_psa = ones([82,totalPatients]);
% Trial data
trial_psa = ones([82,totalPatients]);
trial_and = ones([82,totalPatients]);

options = optimset('Algorithm','interior-point','TolX',1e-13,'TolFun',1e-13,'TolCon'...
    ,1e-13,'MaxIter',1000);                % Optimizer Options
    
%%%%%%%%%%%%%%%%%%%% Creates Structure to store all patient data %%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = LIST;  
   patient = strcat('patient',num2str(ii));   % Patient with corresp number  
   file = strcat('Data/',patient,'.txt');     % Complete name of file patient#.txt
   var = load(file);                          % holds variable just loaded 
   patients(counter) = {var};                 % stores patient data into cell patients
   x = strcat('a',num2str(counter));          % Creates index for structure S
   ind{counter} = x;                          % puts all indexed in ind cell array
   counter = counter+ 1;                      % increases counter 
end
S = cell2struct(patients,ind,2);              % cell to structure

for yy = 1:3 %% Uncomment this to recover the full code for all 4 models.
    % Don't forget to uncomment the "end" also (line ~185).
%yy = 1; %% Comment this out if line 12 (for loop) is uncommented. This runs specifically for three population model.
if yy == 1
    model = 'two_pop';
elseif yy == 2
    model = 'three_pop';
elseif yy == 3
    model = 'portz';
elseif yy == 4
    model = 'hirata';
end

switch model    % Number of parameters, initial cond (xi0) are parameters.      
    case 'two_pop'
        nParams = 22;
    case 'three_pop'
        nParams = 25;
    case 'hirata' 
        nParams = 13; 
    case 'portz' 
        nParams = 25; 
    otherwise 
        warning('Unexpected model, choose a different one')
        return;
end
par_store = zeros(nParams,totalPatients);  % used to store parameter values for each fit 


%%%%%%%%%%%%%%%%%%%%%%%%%%% Runs the fitting for the selected patients %%%%%%%%%%%%%%%%%%%%%

for i = 1:totalPatients                     
    try                               % In case of error for loop will move to next iteration  
    name = ['p',num2str(LIST(i))];    % where parameters will be saved   
    index = char(ind(i));             % Index to calll specific patients
    patient = S.(index);              % calls specific patient from list S
    t = patient(:,2);                 % time vector in days 
    tsize = length(t);                % how many data points there are for this patient
    psa = patient(:,3);               % psa vector of values
    trial_psa(1:tsize,i) = psa';              % save trial psa data for plotting
    androgen = patient(:,4);          % androgen vector of value 
    trial_and(1:tsize,i) = androgen';         % save trial androgen data for plotting
    treatment = patient(:,6);         % treatment vector of values
%%%%%%%%%%%%%%%% Finds periods of on and off treatment %%%%%%%%%%%%%%%%%%%% 
jj = 1; 
change(1) = 1;                        % Treatment starts at t = 0  
for a = 1:length(treatment) 
         if treatment(a) ~= mod(jj,2) % When treatment change occurs    
            jj = jj + 1;
            change(jj) = a;           % Stores time index in change vector 
         end
end
change(jj+1) = length(treatment);     % Last day of treatment 
%%%%%%%%%%%%%%%%%%%%%% Bounds For Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch model 
    case 'two_pop' 
        a = max(androgen);
        b = min(androgen(change(1):change(4)));
         %  um            % q1             % q2              % c1
        LB(1) = 0.01;    LB(2) = b+.1;    LB(3) = 0;        LB(4) = 0.00001;
        UB(1) = 0.1;     UB(2) = b+.5;    UB(3) = b +.1;    UB(4) = .0001;
        %  c2            %  K1            % K2              % b
        LB(5) = 0.00001; LB(6) = 0;       LB(7) = 0;        LB(8) = 0;
        UB(5) = 0.0001;  UB(6) = 1;       UB(7) = 1;        UB(8) = 0.0025;
        %  sigma1        % epsilon        % d1              % d2 
        LB(9) = 0;       LB(10) = 0.01;   LB(11) = 0.002;   LB(12) = 0;
        UB(9) = 1;       UB(10) = 1;      UB(11) = .09;     UB(12) = .001;
        %  R1            %  R2             % gamma1         %  gamma2
        LB(13) = 0;      LB(14) = 0;       LB(15) = 20;     LB(16) = 0;
        UB(13) = 3;      UB(14) = 3;       UB(15) = 20;     UB(16) = .001;
        % dd1               % dd2                 % Qm            
        LB(17) = 0.000001;  LB(18) = 0.000001;    LB(19) = a - 4;  
        UB(17) = .00009;    UB(18) = .00009;      UB(19) = a;      
        % X10              % X20              % u 
        LB(20) = 90;       LB(21) = 0;        LB(22) = 0;
        UB(20) = 100;      UB(21) = 10;       UB(22) = 0;

        x0 = [androgen(1);99;1;psa(1)]; %collects init cond in a column vector 
    case 'three_pop'
        a = max(androgen);
        b = min(androgen(change(1):change(4)));
        %  um            % q1             % q2              % q3
        LB(1) = 0.01;    LB(2) = b+.1;    LB(3) = 0;        LB(4) = 0;
        UB(1) = 0.1;     UB(2) = b+.5;    UB(3) = b +.1;    UB(4) = b +.1;
        %  c             % K              % b               % sigma
        LB(5) = 0.00001; LB(6) = 0;       LB(7) = 0;        LB(8) = 0;
        UB(5) = 0.0001;  UB(6) = 1;       UB(7) = 0.0025;   UB(8) = 1;
        %  epsilon       % d1             % d2              % d3 
        LB(9) = 0.01;    LB(10) = 0.002;  LB(11) = 0;       LB(12) = 0;
        UB(9) = 1;       UB(10) = 0.09;   UB(11) = .001;    UB(12) = .001;
        %  R1            %  R2            % R3              % gamma1
        LB(13) = 0;      LB(14) = 0;      LB(15) = 0;       LB(16) = 18;
        UB(13) = 3;      UB(14) = 3;      UB(15) = 3;       UB(16) = 23;
        %  gamma2        %  u - The way the program is structured, u is the last parameter.         
        LB(17) = 0;      LB(25) = 0;      
        UB(17) = 0.001;  UB(25) = 0;      
        % dd1              % dd2              % dd3              % Qm 
        LB(18) = 0.000001; LB(19) = 0.000001; LB(20) = 0.000001; LB(21) = a - 4;
        UB(18) = .00009;   UB(19) = .00009;   UB(20) = .00009;   UB(21) = a;
        % X10              % X20              % X30 
        LB(22) = 90;       LB(23) = 0;        LB(24) = 0;
        UB(22) = 100;      UB(23) = 5;        UB(24) = 2;
%       % A0                                  % P0
%       LB(26) = 0.5*androgen(1);             LB(27) = 0.5*psa(1);
%       UB(26) = 1.5*androgen(1);             UB(27) = 1.5*psa(1);
        
        x0 = [androgen(1);99;0.9;0.1;psa(1)]; %collects init cond in a column vector 
        %This initial condition vector needs to be updated in the objective
        %function. (Careful since in the end, we'll have 5 models instead
        %so the update process an "if", or similar, statement).
    case 'hirata' 
        % Taken from Everett et. al. 
        % w11o          % w21          % w22o         % w31          % w32      w33o
        LB(1) = -.15;  LB(2) = .0006; LB(3) = -.015; LB(4) = .0003; LB(5) = 0; LB(6) = 0.002;
        UB(1) = -.015; UB(2) = .002;  UB(3) = .0009; UB(4) = .001;  UB(5) = 0; UB(6) = 0.003;
        % w11f           % w12f        % w22f             % w33f
        LB(7) = .001;   LB(8) = .049;  LB(9) = 0.002;    LB(10) = -.13;
        UB(7) = .003;   UB(8) = .18;   UB(9) = 0.008;    UB(10) = -0.0044;
        % X10                   % X20                     % X30 
        LB(11) = psa(1)*.9;     LB(12) = 0;               LB(13) = 0;
        UB(11) = psa(1);        UB(12) = psa(1)*.1;       UB(13) = psa(1)*.05;
        %Possibly need a contraints on initial condition (to add up to P)?
        x0 = [psa(1)*.95;psa(1)*.049;psa(1)*.001]; %collects init cond in a column vector 
    case 'portz' 
        % um           % qx           % qy          % dx           % dy             % c1
        LB(1) = 0.01;  LB(2) = .175;  LB(3) = .1;   LB(4) = 0.15;  LB(5) = 0.215;   LB(6) = .01;
        UB(1) = .1;    UB(2) = .29;   UB(3) = .21;  UB(4) = .4;    UB(5) = .4;      UB(6) = .015;
        % Kxyn         % c2           % Kyxn        % n          % qm          % vm
        LB(7) = .05;   LB(8) = .01;   LB(9) = 1.2;  LB(10) = 1;  LB(11) = 2;   LB(12) = .075;
        UB(7) = .08;   UB(8) = .015;  UB(9) = 1.7;  UB(10) = 1;  UB(11) = 5;   UB(12) = .275;
        % vh          % b             % sigmax      % sigmay       %  rhoxm         % rhoym
        LB(13) = 2;   LB(14) = 0.02;  LB(15) = 0 ;  LB(16) = 0;    LB(17) = .3;     LB(18) = 1;
        UB(13) = 4;   UB(14) = 0.09;  UB(15) = .4;  UB(16) = .4;   UB(17) = 1.3;    UB(18) = 1.3;
        % m           % sigma0        % delta
        LB(19) = 1;   LB(20) = 0;     LB(21) = 0.008;
        UB(19) = 1;   UB(20) = .04;   UB(21) = 0.08;
        % X10          % X20          % Q10         % Q20
        LB(22) = 95;   LB(23) = 0;    LB(24) = 0;   LB(25) = 0;
        UB(22) = 100;  UB(23) = 5;    UB(24) = 1;   UB(25) = 1;
        x0 = [99;1;.5;.5;psa(1)];
end
%%%%%%%%%%%%%%%%%%%%%% Optimization Step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
IC = UB; % initial parameter values  
if yy == 1
    tempCell = struct2cell(load('par_all_two_pop.mat'));
    params = tempCell{1};
elseif yy == 2
    tempCell = struct2cell(load('par_all_three_pop.mat'));
    params = tempCell{1};
elseif yy == 3
    tempCell = struct2cell(load('par_all_portz.mat'));
    params = tempCell{1};
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%% Save model data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[errorAND, errorPSA, Y] = plotData(params,psa,androgen,t,change,n,model);
switch model
    case 'two_pop'
        model_2_psa(1:length(Y),i) = Y(:,4);
        model_2_and(1:length(Y),i) = Y(:,1);
    case 'three_pop'
        model_3_psa(1:length(Y),i) = Y(:,5);
        model_3_and(1:length(Y),i) = Y(:,1);
    case 'portz'
        model_portz_psa(1:length(Y),i) = Y(:,5);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% Save androgen and PSA errors %%%%%%%%%%%%%%%%%%%%%%%%%
if yy~=3
    androgenError(i,yy) = errorAND;
end
psaError(i,yy) = errorPSA;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

par_store(:,i) = params; % stores parameters in a matrix to be used later for predictions and errors
cd(strcat('parameters_',model))
save(name,'params')      % saves parameters for individual patients in file p(patient#).mat
cd ..
    catch ME
         disp(ME)        % If an error occurs it will be displayed

    end

end
%%
save(strcat('par_all_',model,'.mat'),'par_store') % saves the matrix of all patients parameter values
end  % end for the loop going over all the models 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tIndex = change(n+1);     % Index at which the n cycles end
fpath = 'C:\Users\Alejandro Martinez\Documents\Mathematics\ProstateCancer\ModifiedCode\Plots'; % Plot destination

% Plot androgen for model 2 and 3 for each patient
for j= 1:totalPatients
    patientName = strcat('Patient ',num2str(LIST(i)));                  % Create string for patient number used in plot
    fileNameAndrogen = strcat('Patient',num2str(LIST(i)),'Androgen');   % Create string for androgen plot file name
    fileNamePSA = strcat('Patient',num2str(LIST(i)),'PSA');             % Create string for PSA plot file name
    
    scatter(t(1:tIndex), trial_and(1:tIndex,j));                        % Plot trial data
    hold on;
    plot(t(1:tIndex), model_2_and(1:tIndex,j));                         % Plot model 2 data
    hold on;
    plot(t(1:tIndex), model_3_and(1:tIndex,j));                         % Plot model 3 data
    hold off;
    title(patientName, 'fontsize', 26);
    xlabel('Time (Days)', 'fontsize', 26);
    ylabel('Androgen (nM)', 'fontsize', 26);
    legend('Data','Model2','Model3');
    saveas(gcf, fullfile(fpath, fileNameAndrogen), 'jpeg');
    
    scatter(t(1:tIndex), trial_psa(1:tIndex,j));                        % Plot trial data
    hold on;
    plot(t(1:tIndex), model_2_psa(1:tIndex,j));                         % Plot model 2 data
    hold on;
    plot(t(1:tIndex), model_3_psa(1:tIndex,j));                         % Plot model 3 data
    hold on;
    plot(t(1:tIndex), model_portz_psa(1:tIndex,j));                     % Plot PKN Data
    hold off;
    title(patientName, 'fontsize', 26);
    xlabel('Time (Days)', 'fontsize', 26);
    ylabel('PSA (\mug/L)', 'fontsize', 26);
    legend('Data','Model2','Model3','PKN');
    saveas(gcf, fullfile(fpath, fileNamePSA), 'jpeg');
end

end

%%%%%%%%%%%%%%%%%  Retrieve error and model data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [erra, errp, Y1run]= plotData(params,psadata,and_data,tdata,change,n,model) 
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
% x0 = [and_data(1);params(22);params(23);params(24);psadata(1)];
% [Y1run] = run_model(params,tdata,change,x0,n,model,and_data);
erra = 0;      %initialize androgen error so that there is some value for PKN model

switch model        
    case 'two_pop' 
        x0 = [and_data(1);params(20);params(21);psadata(1)];
        [Y1run] = run_model(params,tdata,change,x0,n,model,and_data);        
        %
        PSA = Y1run(:,4);
        AND = Y1run(:,1);
        errp = sum((PSA-psadata).^2/length(PSA));
        erra = sum((AND-and_data).^2/length(AND));
%        
    case 'three_pop' 
        x0 = [and_data(1);params(22);params(23);params(24);psadata(1)];
        [Y1run] = run_model(params,tdata,change,x0,n,model,and_data);
        %
        PSA = Y1run(:,5);
        AND = Y1run(:,1);
        errp = sum((PSA-psadata).^2/length(PSA));
        erra = sum((AND-and_data).^2/length(AND));       
%        
    case 'portz'
        x0 = [params(22);params(23);params(24);params(25);psadata(1)];
        [Y1run] = run_model(params,tdata,change,x0,n,model,and_data);
        PSA = Y1run(:,5);
        errp = sum((PSA-psadata).^2/length(PSA));
%        
    case 'hirata' 
        x0 = [params(11);params(12);params(13)];
        [Y1run] = run_model(params,tdata,change,x0,n,model,and_data);
        %
        PSA = Y1run(:,1) + Y1run(:,2) + Y1run(:,3);
        errp = sum((PSA-psadata).^2/length(PSA));

end
end

%%%%%%%%%%%%%%%%%  Minimizes the Error %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [err]= objective(params,psadata,and_data,tdata,change,~,n,model) 
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
% x0 = [and_data(1);params(22);params(23);params(24);psadata(1)];
% [Y1run] = run_model(params,tdata,change,x0,n,model,and_data);

switch model        
    case 'two_pop' 
        x0 = [and_data(1);params(20);params(21);psadata(1)];
        [Y1run] = run_model(params,tdata,change,x0,n,model,and_data);        
        %
        PSA = Y1run(:,4);
        AND = Y1run(:,1);
        errp = sum((PSA-psadata).^2/length(PSA));
        erra = sum((AND-and_data).^2/length(AND));
        err = errp + erra;
        fprintf('PSA Error = %.4f \t Androgen Error = %.4f \n',errp, erra);
%        
    case 'three_pop' 
        x0 = [and_data(1);params(22);params(23);params(24);psadata(1)];
        [Y1run] = run_model(params,tdata,change,x0,n,model,and_data);
        %
        PSA = Y1run(:,5);
        AND = Y1run(:,1);
        errp = sum((PSA-psadata).^2/length(PSA));
        erra = sum((AND-and_data).^2/length(AND));
        err = errp + erra;
        fprintf('PSA Error = %.4f \t Androgen Error = %.4f \n',errp, erra);        
%        
    case 'portz'
        x0 = [params(22);params(23);params(24);params(25);psadata(1)];
        [Y1run] = run_model(params,tdata,change,x0,n,model,and_data);
        %
        PSA = Y1run(:,5);
        err = sum((PSA-psadata).^2/length(PSA));
        fprintf('PSA Error = %.4f \n',err);
%        
    case 'hirata' 
        x0 = [params(11);params(12);params(13)];
        [Y1run] = run_model(params,tdata,change,x0,n,model,and_data);
        %
        PSA = Y1run(:,1) + Y1run(:,2) + Y1run(:,3);
        err = sum((PSA-psadata).^2/length(PSA));
        fprintf('PSA Error = %.4f \n',err);

end
end

%%
%%%%%%%%%%%%%%%%  Runs the Model and Generates synthetic data %%%%%%%%%%%%
function [y] = run_model(params,tdata,change,x0,n,model,androgen)   
y = [];                                      % Initial vector for solution
tint = tdata(change(n+1));                   % Time interval to run solution 
for k = 1:n                                  % done to cycle though u values
    u = 1 - mod(k,2);    
    switch model             
        case 'two_pop'
            if tint(end) >= tdata(change(k))
                [~,Yrun]=ode15s(@(t,x) two_pop(t,x,[params(1:end-1),u]),tdata(change(k):change(k+1)),x0);
                x0 = [Yrun(end,1);Yrun(end,2);Yrun(end,3);Yrun(end,4)];
                if k < n
                    y = [y; Yrun(1:end-1,:)]; %#ok<AGROW>
                elseif k == n
                    y = [y; Yrun(1:end,:)];   %#ok<AGROW>
                end
            end
%            
        case 'three_pop'
            if tint(end) >= tdata(change(k))
                [~,Yrun]=ode15s(@(t,x) three_pop(t,x,[params(1:end-1),u]),tdata(change(k):change(k+1)),x0);
                x0 = [Yrun(end,1);Yrun(end,2);Yrun(end,3);Yrun(end,4);Yrun(end,5)];
                if k < n
                    y = [y; Yrun(1:end-1,:)]; %#ok<AGROW>
                elseif k == n
                    y = [y; Yrun(1:end,:)];   %#ok<AGROW>
                end
            end
%            
        case 'portz'
            if tint(end) >= tdata(change(k))
%                 options = odeset('AbsTol',1e-14,'RelTol',1e-14);
                [~,Yrun]=ode23tb(@(t,x) portz(t,x,params,androgen,change(k):change(k+1)),tdata(change(k):change(k+1)),x0);
                x0 = [Yrun(end,1);Yrun(end,2);Yrun(end,3);Yrun(end,4);Yrun(end,5)];
                if k < n
                    y = [y; Yrun(1:end-1,:)]; %#ok<AGROW>
                elseif k == n
                    y = [y; Yrun(1:end,:)];   %#ok<AGROW>
                end
            end 
%            
        case 'hirata' 
            if tint(end) >= tdata(change(k))
                [~,Yrun]=ode15s(@(t,x) hirata(t,x,params,u),tdata(change(k):change(k+1)),x0);
                x0 = [Yrun(end,1);Yrun(end,2);Yrun(end,3)];
                if k < n
                    y = [y; Yrun(1:end-1,:)]; %#ok<AGROW>
                elseif k == n
                    y = [y; Yrun(1:end,:)];   %#ok<AGROW>
                end
            end
    end   
end
end

%% model ODE functions

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
%%
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
%%

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