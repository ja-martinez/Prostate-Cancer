function graph

LIST = [1,2,6,7,12,14:17,19,24:25,28,29,31,32,36:37,39:40,42,44,51:52,54,55,58,60:64,66,75 ...
        ,77:79,83:85,87,88,91,93:97,99:102,104:109];    % load patient array to know parameter indexing
    load('par_all_three_pop.mat','par_store');          % load parameters
    patient = 15;
    patientString = '15';
    patientIndex = find(LIST==patient);                      % find index for patient_
    parameters = par_store(:,patientIndex);             % store patient_ parameters                        
    file = strcat('Data_62/patient',patientString,'.txt');    % define file name
    data = load(file);                                 % load patient_ data 
    psa_data = data(:,3);                             % retrieve initial PSA
    and_data = data(:,4);                             % retrieve initial androgen
    
    
    
end