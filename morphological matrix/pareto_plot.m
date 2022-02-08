%% Clear

clear
clc;

%% load payload package data

% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 11);

% Specify sheet and range
opts.Sheet = "Payload Package";
opts.DataRange = "A3:K6";

% Specify column names and types
opts.VariableNames = ["ID", "Option", "Reliability1", "Reliability2", "Reliability3", "Cost1", "Cost2", "Cost3", "Value1", "Value2", "Value3"];
opts.VariableTypes = ["double", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify variable properties
opts = setvaropts(opts, "Option", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Option", "EmptyFieldRule", "auto");

% Import the data
payload_package = readtable("Master Morphological Matrix.xlsx", opts, "UseExcel", false);

% Clear temporary variables
clear opts

%% load second stage propulsion data

% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 11);

% Specify sheet and range
opts.Sheet = "Second Stage Propulsion";
opts.DataRange = "A3:K6";

% Specify column names and types
opts.VariableNames = ["ID", "Option", "Reliability1", "Reliability2", "Reliability3", "Cost1", "Cost2", "Cost3", "Value1", "Value2", "Value3"];
opts.VariableTypes = ["double", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify variable properties
opts = setvaropts(opts, "Option", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Option", "EmptyFieldRule", "auto");

% Import the data
second_stage = readtable("Master Morphological Matrix.xlsx", opts, "UseExcel", false);

% Clear temporary variables
clear opts

%% load lift-off stage data

% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 11);

% Specify sheet and range
opts.Sheet = "Lift-off Stage";
opts.DataRange = "A3:K9";

% Specify column names and types
opts.VariableNames = ["ID", "Option", "Reliability1", "Reliability2", "Reliability3", "Cost1", "Cost2", "Cost3", "Value1", "Value2", "Value3"];
opts.VariableTypes = ["double", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify variable properties
opts = setvaropts(opts, "Option", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Option", "EmptyFieldRule", "auto");

% Import the data
lift_off_stage = readtable("Master Morphological Matrix.xlsx", opts, "UseExcel", false);

% Clear temporary variables
clear opts

%% load inclination change data

% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 11);

% Specify sheet and range
opts.Sheet = "Inclination Change";
opts.DataRange = "A3:K8";

% Specify column names and types
opts.VariableNames = ["ID", "Option", "Reliability1", "Reliability2", "Reliability3", "Cost1", "Cost2", "Cost3", "Value1", "Value2", "Value3"];
opts.VariableTypes = ["double", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify variable properties
opts = setvaropts(opts, "Option", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Option", "EmptyFieldRule", "auto");

% Import the data
inclination_change = readtable("Master Morphological Matrix.xlsx", opts, "UseExcel", false);

% Clear temporary variables
clear opts

%% load final operating orbit data

% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 11);

% Specify sheet and range
opts.Sheet = "Final Operating Orbit";
opts.DataRange = "A3:K9";

% Specify column names and types
opts.VariableNames = ["ID", "Option", "Reliability1", "Reliability2", "Reliability3", "Cost1", "Cost2", "Cost3", "Value1", "Value2", "Value3"];
opts.VariableTypes = ["double", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify variable properties
opts = setvaropts(opts, "Option", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Option", "EmptyFieldRule", "auto");

% Import the data
final_orbit = readtable("Master Morphological Matrix.xlsx", opts, "UseExcel", false);

% Clear temporary variables
clear opts

%% create permuation matrix

% number of options of each category
N_pay = length(payload_package.Option);
N_sec = length(second_stage.Option);
N_lift = length(lift_off_stage.Option);
N_incl = length(inclination_change.Option);
N_orbit = length(final_orbit.Option);

% total number of possible configurations
N_total = N_pay*N_sec*N_lift*N_incl*N_orbit;

% permuation structure
design_space(N_total) = struct;

% permutation matrix
dsi = zeros(N_total,5);

% figure of merit vectors
reliabilities = zeros(N_total,3);
costs = zeros(N_total,3);
values = zeros(N_total,3);

% loop over
n = 1;
for i = 1:N_pay
    for j = 1:N_sec
        for k = 1:N_lift
            for l = 1:N_incl
                for m = 1:N_orbit

                    % design
                    design_space(n).payload_package = payload_package.Option(j);
                    design_space(n).second_stage = second_stage.Option(j);
                    design_space(n).lift_off_stage = lift_off_stage.Option(k);
                    design_space(n).inclination_change = inclination_change.Option(l);
                    design_space(n).final_orbit = final_orbit.Option(m);

                    dsi(n,:) = [i j k l m];

                    % reliabilities
                    reliabilities(n,1) = (payload_package.Reliability1(j) + ...
                        second_stage.Reliability1(j) + ...
                        lift_off_stage.Reliability1(k) + ...
                        inclination_change.Reliability1(l) + ...
                        final_orbit.Reliability1(m))/5;

                    reliabilities(n,2) = (payload_package.Reliability2(j) + ...
                        second_stage.Reliability2(j) + ...
                        lift_off_stage.Reliability2(k) + ...
                        inclination_change.Reliability2(l) + ...
                        final_orbit.Reliability2(m))/5;

                    reliabilities(n,3) = (payload_package.Reliability3(j) + ...
                        second_stage.Reliability3(j) + ...
                        lift_off_stage.Reliability3(k) + ...
                        inclination_change.Reliability3(l) + ...
                        final_orbit.Reliability3(m))/5;

                    % costs
                    costs(n,1) = (payload_package.Cost1(j) + ...
                        second_stage.Cost1(j) + ...
                        lift_off_stage.Cost1(k) + ...
                        inclination_change.Cost1(l) + ...
                        final_orbit.Cost1(m))/5;

                    costs(n,2) = (payload_package.Cost2(j) + ...
                        second_stage.Cost2(j) + ...
                        lift_off_stage.Cost2(k) + ...
                        inclination_change.Cost2(l) + ...
                        final_orbit.Cost2(m))/5;

                    costs(n,3) = (payload_package.Cost3(j) + ...
                        second_stage.Cost3(j) + ...
                        lift_off_stage.Cost3(k) + ...
                        inclination_change.Cost3(l) + ...
                        final_orbit.Cost3(m))/5;

                    % value
                    values(n,1) = (payload_package.Value1(j) + ...
                        second_stage.Value1(j) + ...
                        lift_off_stage.Value1(k) + ...
                        inclination_change.Value1(l) + ...
                        final_orbit.Value1(m))/5;

                    values(n,2) = (payload_package.Value2(j) + ...
                        second_stage.Value2(j) + ...
                        lift_off_stage.Value2(k) + ...
                        inclination_change.Value2(l) + ...
                        final_orbit.Value2(m))/5;

                    values(n,3) = (payload_package.Value3(j) + ...
                        second_stage.Value3(j) + ...
                        lift_off_stage.Value3(k) + ...
                        inclination_change.Value3(l) + ...
                        final_orbit.Value3(m))/5;


                    n = n + 1;
                end
            end
        end
    end
end

%% find infeasible designs

% solar sails conflict with:
%   inclination change:
%       direct transfer
%       bielliptic transfer
%   
tmp = dsi(:,2) == 1 & (dsi(:,4) == 1 | dsi(:,4) == 2);

% chemical, ion, and nuclear conflict with:
%   inclination change:
%       solar sail
%       venus GA + solar sail
%   final orbit:
%       circular orbit (cranking)
%       non-keplarian
%       polaris final
%

tmp = dsi(:,)
%% evaluate the merits of the design

% % data to pass into function
% data.payload_package = payload_package;
% data.second_stage = second_stage;
% data.lift_off_stage = lift_off_stage;
% data.inclination_change = inclination_change;
% data.final_orbit = final_orbit;
% 
% % initialize arrays
% reliabilities = zeros(N_total,1);
% costs = zeros(N_total,1);
% values = zeros(N_total,1);
% 
% % iterate over all designs
% for i = 1:N_total
%     [reliabilities(i), costs(i), values(i)] = calc_merit(design_space(i),data);
% end
% 
% function [reliability, cost, value] = calc_merit(design,data)
% 
% 
% reliability = nan;
% cost = nan;
% value = nan;
% 
% if design.inclination_change == "Direct Transfer from 1 AU" || ...
%         design.inclination_change == "Bielliptic Transfer" || ...
%         design.inclination_change == "Gravity Assist Inner Planets" || ...
%         design.inclination_change == "Gravity Assist Outer Planets"
% 
%     if design.second_stage == "Solar Sail" || design.second_stage == "Ion"
%         return
%     end
% end
% 
% if design.inclination_change == "Solar Sail" || ...
%         design.inclination_change == "Venus GA + Solar Sail"
% 
%     if design.second_stage ~= "Solar Sail"
%         return
%     end
% end
% 
% if design.final_orbit == "Elliptical Orbit (Venus flyby)"
% 
%     if design.inclination_change ~= "Gravity Assist Inner Planets" && ...
%             design.inclination_change ~= "Venus GA + Solar Sail"
%         return
%     end
% end
% 
% if design.final_orbit == "Elliptical Orbit (Earth flyby)" || ...
%         "Elliptical Orbit (Mars flyby)"
% 
%     if design.inclination_change ~= "Gravity Assist Inner Planets"
%         return
%     end
% end
% 
% if design.final_orbit == "Modified Ulysses orbit (Jupiter flyby)"
%     if design.inclination_change ~= "Gravity Assist Outer Planets"
%         return
%     end
% end
% 
% if design.final_orbit == "Circular orbit (Spiral + Orbit Cranking)" || ...
%         design.final_orbit == "Non-Keplarian Orbit (beta=0.8)" || ...
%         design.final_orbit == "POLARIS Final Orbit"
%     if design.second_stage ~= "Solar Sail"
%         return
%     end
% end
% 
% reliability = payload_package(design.payload_package == data.payload_package.Option,'Reliability1') + ...
%     second_stage(design.second_stage == data.second_stage.Option,'Reliability1') + ...
%     lift_off_stage(design.lif)
% 
% 
% end
