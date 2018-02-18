function bifurcation
    % Creates bifurcation plot 
    % to change patient, chage patient and patientString variable
    % to change parameter, change the upper and lower bounds, and the
    % parameters index inside loop

    LIST = [1,2,6,7,12,14:17,19,24:25,28,29,31,32,36:37,39:40,42,44,51:52,54,55,58,60:64,66,75 ...
        ,77:79,83:85,87,88,91,93:97,99:102,104:109];    % load patient array to know parameter indexing
    load('par_all_three_pop.mat','par_store');          % load parameters
    patient = 15;
    patientString = '15';
    patientIndex = find(LIST==patient);                      % find index for patient_
    parameters = par_store(:,patientIndex);             % store patient_ parameters                        
    file = strcat('Data_62/patient',patientString,'.txt');    % define file name
    data = load(file);                                 % load patient_ data 
    initialPSA = data(1,3);                             % retrieve initial PSA
    initialAnd = data(1,4);                             % retrieve initial androgen
    
    tLength = 5*1008;                  % length of treatment
    t = 0:1:tLength;                 % create time array
    tCheckup = 14;                   % number of days between doctor visits
    numVisits = tLength/tCheckup;    % total number of visits to solve in parts
    tStart = 4.5*1008;                    % when to start values to include in bifurcation diagram
    upperPSA = 10;                   % upper PSA threshold
    lowerPSA = 1;                    % lower PSA threshold
    
    % I will look at PSA values with respect to changes in eps
    lowerBound = 0.001;                                  % lower bound for chosen variable to vary
    upperBound = 0.01;                                   % upper bound for chosen variable to vary
    numValues = 1000;                                    % total number of values to evaluate
    valInterval = (upperBound - lowerBound)/numValues;    % length between interval between values
    values = lowerBound:valInterval:upperBound;          % vector holding all values to be evaluated
    
    x0 = [initialAnd,parameters(22),parameters(23),parameters(24),initialPSA];    % initial conditions
    PSA = zeros(tLength,length(values));          % initialize array to store PSA values
    maxVal = zeros(1,length(values));           % initialize array that stores max val for each value
    minVal = zeros(1,length(maxVal));           % initialize array that stores min val for each value
    
    % loop for all chosen values of the parameter you will vary
    for index = 1:numValues+1
        value = values(index);      % retrieve value of parameter
        parameters(9) = value;     % evaluate at specific variable value
        u = 0;                      % u=0 on treatment and u=1 off treatment
    
        for visit = 1:numVisits     % solve in parts to change 
            
            % solve model for specific set of parameters
            [~, Y]=ode15s(@(t,x) three_pop(t,x,[parameters(1:end-1)',u]),t(14*visit-13):t(14*visit),x0);
            
            weekPSA = Y(:,5);                       % extract PSA values
            PSA(14*visit-13:14*visit,index) = weekPSA;    % save PSA values
            
            % if PSA goes out of bounds, change treatment (u)
            if weekPSA(end) >= upperPSA
                u = 0;
            elseif weekPSA(end) <= lowerPSA
                u = 1;
            end
            
            x0 = Y(end,:);      % change initial conditions for the next iteration to follow this one
        end 
        
        maxVal(index) = max(PSA(tStart:end,index));
        minVal(index) = min(PSA(tStart:end,index));
        
    end
    
    % plot
    figure
    plot(values,maxVal);
    hold on;
    plot(values,minVal);
    legend('max values','min values');
    xlabel('\epsilon (day^-^1)');
    ylabel('PSA \mu g/L');

    
end

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
%
    ux = abs(um)*(1 - abs(q1)/Q);
%
    uy = abs(um)*(1 - abs(q2)/Q);
%
    uz = abs(um)*(1 - abs(q3)/Q);
%
Dx = abs(d1)*R1/(Q+R1);    
Dy = abs(d2)*R2/(Q+R2);
Dz = abs(d3)*R3/(Q+R3);
Lm = abs(c)*K/(Q + K); % "Mutation" rate
gamma = abs(gamma1)*u +abs(gamma2);
%%%%%%%%%%%%%%%%%%% ODE system %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dX = (ux - Dx  - abs(dd1)*X)*X - 2*Lm*X;
dY = (uy - Dy  - abs(dd2)*Y)*Y + Lm*X - 2*Lm*Y;
dZ = (uz - Dz - abs(dd3)*Z)*Z + Lm*(X+Y) ;
dA = gamma*(Qm -Q) - (ux*Q*X + uy*Q*Y +uz*Q*Z)/(X+Y+Z);
dP = b*Q + sigma*(Y*Q + X*Q + Z*Q) - epsilon*P;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dxdt = [dA; dX; dY; dZ; dP]; %Puts the ode in a column vector
end