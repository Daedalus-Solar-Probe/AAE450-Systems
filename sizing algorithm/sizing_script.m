%% Clear

clear
clc

%% Inputs

% system design variable
design.launch_vehicle.name = "Falcon Heavy (Expendable)";
design.launch_vehicle.capacity = [15010 12345 10115 8225 6640 5280 4100 3080 2195 1425 770]; % [kg]
design.launch_vehicle.C3 = 0:10:100; % [km^2/s^2]