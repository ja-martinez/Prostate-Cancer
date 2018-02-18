classdef Patients < handle
    %%UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    %%
    properties
        time_data
        psa_data
        androgen_data
        treatment_data
        x0_2
        x0_3
        LB_2
        UB_2
        LB_3
        UB_3
        change                  % array stores index of where treatment changes
        parameters_2
        parameters_3
        nFitting = 3
        nForecasting = 2
    end
    %%
    methods
        function obj = Patients(patientNumber)
            % Class constructor
            % Loads patient clinical data and parameter ranges into object properties
            
            
            %% load data
            patient = strcat('patient',num2str(patientNumber));   % Patient with corresp number  
            file = strcat('Data_62/',patient,'.txt');             % Complete name of file patient#.txt
            data = load(file);                                    % loads patient data
            
            %% put data in object properties 
            obj.time_data = data(:,2);                            
            obj.psa_data = data(:,3);                             
            obj.androgen_data = data(:,4);     
            obj.treatment_data = data(:,6);        
       
            %% Assign values to change vector
            jj = 1;
            obj.change = 1;                        % Treatment starts at t = 0                        
            for i = 2:length(obj.treatment_data)
                if obj.treatment_data(i) ~= mod(jj,2)         % When treatment change occurs
                    jj = jj + 1;
                    obj.change(jj) = i;                  % Stores time in change vector
                end
            end

            if length(obj.treatment_data) ~= obj.change(jj) % add this in case the last treatment is only one measurement
                obj.change(jj+1) = length(obj.treatment_data);     % Last day of treatment
            end
            
            
            %% two pop parameter ranges
            a = max(obj.androgen_data);
            b = min(obj.androgen_data(obj.change(1):obj.change(obj.nFitting)));
             %  um                    % q1                     % q2                      % c1
            obj.LB_2(1) = 0.01;    obj.LB_2(2) = b+.1;    obj.LB_2(3) = 0;        obj.LB_2(4) = 0.00001;
            obj.UB_2(1) = 0.1;     obj.UB_2(2) = b+.5;    obj.UB_2(3) = b +.1;    obj.UB_2(4) = .0001;
            %  c2                     %  K1                    % K2                      % b
            obj.LB_2(5) = 0.00001; obj.LB_2(6) = 0;       obj.LB_2(7) = 0;        obj.LB_2(8) = 0;
            obj.UB_2(5) = 0.0001;  obj.UB_2(6) = 1;       obj.UB_2(7) = 1;        obj.UB_2(8) = 0.0025;
            %  sigma1                 % epsilon                % d1                      % d2 
            obj.LB_2(9) = 0;       obj.LB_2(10) = 0.01;   obj.LB_2(11) = 0.002;   obj.LB_2(12) = 0;
            obj.UB_2(9) = 1;       obj.UB_2(10) = 1;      obj.UB_2(11) = .09;     obj.UB_2(12) = .001;
            %  R1                     %  R2                     % gamma1                 %  gamma2
            obj.LB_2(13) = 0;      obj.LB_2(14) = 0;       obj.LB_2(15) = 20;     obj.LB_2(16) = 0;
            obj.UB_2(13) = 3;      obj.UB_2(14) = 3;       obj.UB_2(15) = 20;     obj.UB_2(16) = .001;
            % dd1                        % dd2                         % Qm            
            obj.LB_2(17) = 0.000001;  obj.LB_2(18) = 0.000001;    obj.LB_2(19) = a - 4;  
            obj.UB_2(17) = .00009;    obj.UB_2(18) = .00009;      obj.UB_2(19) = a;      
            % X10                       % X20                      % u 
            obj.LB_2(20) = 90;       obj.LB_2(21) = 0;        obj.LB_2(22) = 0;
            obj.UB_2(20) = 100;      obj.UB_2(21) = 10;       obj.UB_2(22) = 0;

            %% three pop parameter ranges
            a = max(obj.androgen_data);
            b = min(obj.androgen_data(obj.change(1):obj.change(obj.nFitting)));
             %  um                       % q1                       % q2                        % q3
            obj.LB_3(1) = 0.01;    obj.LB_3(2) = b+.1;    obj.LB_3(3) = 0;        obj.LB_3(4) = 0;
            obj.UB_3(1) = 0.1;     obj.UB_3(2) = b+.5;    obj.UB_3(3) = b +.1;    obj.UB_3(4) = b +.1;
             %  c                        % K                        % b                         % sigma
            obj.LB_3(5) = 0.00001; obj.LB_3(6) = 0;       obj.LB_3(7) = 0;        obj.LB_3(8) = 0;
            obj.UB_3(5) = 0.0001;  obj.UB_3(6) = 1;       obj.UB_3(7) = 0.0025;   obj.UB_3(8) = 1;
             %  epsilon                  % d1                       % d2                        % d3 
            obj.LB_3(9) = 0.01;    obj.LB_3(10) = 0.002;  obj.LB_3(11) = 0;       obj.LB_3(12) = 0;
            obj.UB_3(9) = 1;       obj.UB_3(10) = 0.09;   obj.UB_3(11) = .001;    obj.UB_3(12) = .001;
             %  R1                       %  R2                      % R3                        % gamma1
            obj.LB_3(13) = 0;      obj.LB_3(14) = 0;      obj.LB_3(15) = 0;       obj.LB_3(16) = 18;
            obj.UB_3(13) = 3;      obj.UB_3(14) = 3;      obj.UB_3(15) = 3;       obj.UB_3(16) = 23;
             %  gamma2                   %  u - The way the program is structured, u is the last parameter.         
            obj.LB_3(17) = 0;      obj.LB_3(25) = 0;      
            obj.UB_3(17) = 0.001;  obj.UB_3(25) = 0;      
             % dd1                         % dd2                        % dd3                        % Qm 
            obj.LB_3(18) = 0.000001; obj.LB_3(19) = 0.000001; obj.LB_3(20) = 0.000001; obj.LB_3(21) = a - 4;
            obj.UB_3(18) = .00009;   obj.UB_3(19) = .00009;   obj.UB_3(20) = .00009;   obj.UB_3(21) = a;
             % X10                         % X20                        % X30 
            obj.LB_3(22) = 90;       obj.LB_3(23) = 0;        obj.LB_3(24) = 0;
            obj.UB_3(22) = 100;      obj.UB_3(23) = 5;        obj.UB_3(24) = 2;
   %         % A0                                            % P0
   %        obj.LB_3(26) = 0.5*androgen(1);             obj.LB_3(27) = 0.5*psa(1);
   %        obj.UB_3(26) = 1.5*androgen(1);             obj.UB_3(27) = 1.5*psa(1);
        end
        function fit(obj,cycleStart,cycleFinish,errType,model,accuracy)
%%            does fitting for a specific interval and object changes parameters
%             cycleStart = cycle in which you wish to start fitting
%             cycleFinish = cycle in which you wish to end fitting (it will finish this particular cycle)
%             errType = choose either 'MSE' or 'relative'
%             model = choose either ''two_pop' or 'three_pop'
%             accuracy = accuracy of optimization (good results happen ~e-13

            %% optimizer options
            options = optimset('Algorithm','interior-point','TolX',accuracy,'TolFun',accuracy,'TolCon'...
                ,accuracy,'MaxIter',1000);                %  orig: 1e-13
            
            switch model
                case 'two_pop'
                    %% find optiized parameters for two pop model
                    IC = obj.UB_2; % initial parameter values  
                    [obj.parameters_2,~] = fmincon(@(params)objective(obj,params,cycleStart,cycleFinish....
                        ,model,errType),IC,[],[],[],[],obj.LB_2,obj.UB_2,[],options);
                
                case 'three_pop'
                    %% find optiized parameters for three pop model
                    IC = obj.UB_3; % initial parameter values  
                    [obj.parameters_3,~] = fmincon(@(params)objective(obj,params,cycleStart,cycleFinish...
                        ,model,errType),IC,[],[],[],[],obj.LB_3,obj.UB_3,[],options);
            end
        end
        function err = objective(obj,params,cycleStart,cycleFinish,model,errType)
        %% This is the objective function that gives out the corresponding 
        % error given a set of parameters
        switch model
            case 'two_pop'
                [Y] = run_model(obj,params,cycleStart,cycleFinish,model);
                
                PSA = Y(:,4);
                AND = Y(:,5);
                
                switch errType
                    case 'MSE'
                        err = mseError(obj,PSA,AND);
                        
                    case 'relative'
                        err = relativeError(obj,PSA,AND);                        
                end
                
            case 'three_pop'
                [Y] = run_model(obj,params,cycleStart,cycleFinish,model);
                
                PSA = Y(:,5);
                AND = Y(:,1);
                
                switch errType
                    case 'MSE'
                        err = mseError(obj,PSA,AND);
                        
                    case 'relative'
                        err = relativeError(obj,PSA,AND);                        
                end
        end
        end
        function y = run_model(obj,params,cycleStart,cycleFinish,model)
%%          solves the model equations from given parameters and time interval
%           outputs vector of variable values            
            % get last time value in this interval
            tint = tdata(obj.change(cycleFinish+1));
            
            % loop through cycles 
            for k = cycleStart:cycleFinish
                u = 1 - mod(k,2);
                
                switch model
                    case 'two_pop'
                            [~,Yrun]=ode15s(@(t,x) two_pop(t,x,[params(1:end-1),u]),tdata(obj.change(k):obj.change(k+1)),obj.x0);
                            x0 = [Yrun(end,1);Yrun(end,2);Yrun(end,3);Yrun(end,4)];
                            if k < n
                                y = [y; Yrun(1:end-1,:)]; %#ok<AGROW>
                            elseif k == n
                                y = [y; Yrun(1:end,:)];   %#ok<AGROW>
                            end
                        
                    case 'three_pop'
                        
                        
                end
            end
            
        end
        function err = mseError(obj,PSA,AND)
            %% function calculates the MSE error given PSA and AND values
            
            errp = sum((PSA-obj.psa_data).^2/length(PSA));          % psa error
            erra = sum((AND-obj.androgen_data).^2/length(AND));     % androgen error
            err = errp + erra;                                      % total error
        end
        function err = relativeError(obj,PSA,AND)
            %% function calculates the relative error given PSA and AND values
            
            errp = sum((PSA-obj.psa_data)/obj.psa_data);            % psa error
            erra = sum((AND-obj.androgen_data)/obj.androgen_data);  % androgen error
            err = errp + erra;                                      % total error
            
        end
        function x0 = change.x0_2(obj)
            %% set initial conditions for two pop
            obj.x0_2 = [obj.androgen_data(1);99;1;obj.psa_data(1)];         
            
        end   
        function x0 = change.x0_3(obj,)
            %% set initial conditions for three pop
            x0 = [obj.androgen_data(1);99;0.9;0.1;obj.psa_data(1)];  
        end
        
    end
    
end

