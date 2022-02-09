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

% % modify data
% lift_off_stage.Cost1 = log(lift_off_stage.Cost1) + 1;
% lift_off_stage.Cost2 = log(lift_off_stage.Cost2) + 1;
% lift_off_stage.Cost3 = log(lift_off_stage.Cost3) + 1;

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

% % modify data
% inclination_change.Cost1 = log(inclination_change.Cost1) + 1;
% inclination_change.Cost2 = log(inclination_change.Cost2) + 1;
% inclination_change.Cost3 = log(inclination_change.Cost3) + 1;


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

%% Inputs

%%%%%%%%%%%% Weights %%%%%%%%%%%%%

% payload
cw1 = 1/5; % cost
vw1 = 1/5; % mission value
rw1 = 1/5; % success/reliability

% second stage propulsion system cost weight
cw2 = 1/5; % cost
vw2 = 1/5; % mission value
rw2 = 1/5; % success/reliability

% lift off stage cost weight
cw3 = 1/5; % cost
vw3 = 1/5; % mission value
rw3 = 1/5; % success/reliability

% inclination change method cost weight
cw4 = 1/5; % cost
vw4 = 1/5; % mission value
rw4 = 1/5; % success/reliability

% final orbit cost weight
cw5 = 1/5; % cost
vw5 = 1/5; % mission value
rw5 = 1/5; % success/reliability

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
                    design_space(n).payload_package = payload_package.Option(i);
                    design_space(n).second_stage = second_stage.Option(j);
                    design_space(n).lift_off_stage = lift_off_stage.Option(k);
                    design_space(n).inclination_change = inclination_change.Option(l);
                    design_space(n).final_orbit = final_orbit.Option(m);

                    dsi(n,:) = [i j k l m];

                    % reliabilities
                    reliabilities(n,1) = rw1*payload_package.Reliability1(i) + ...
                        rw2*second_stage.Reliability1(j) + ...
                        rw3*lift_off_stage.Reliability1(k) + ...
                        rw4*inclination_change.Reliability1(l) + ...
                        rw5*final_orbit.Reliability1(m);

                    reliabilities(n,2) = rw1*payload_package.Reliability2(i) + ...
                        rw2*second_stage.Reliability2(j) + ...
                        rw3*lift_off_stage.Reliability2(k) + ...
                        rw4*inclination_change.Reliability2(l) + ...
                        rw5*final_orbit.Reliability2(m);

                    reliabilities(n,3) = rw1*payload_package.Reliability3(i) + ...
                        rw2*second_stage.Reliability3(j) + ...
                        rw3*lift_off_stage.Reliability3(k) + ...
                        rw4*inclination_change.Reliability3(l) + ...
                        rw5*final_orbit.Reliability3(m);

                    % costs
                    costs(n,1) = cw1*payload_package.Cost1(i) + ...
                        cw2*second_stage.Cost1(j) + ...
                        cw3*lift_off_stage.Cost1(k) + ...
                        cw4*inclination_change.Cost1(l) + ...
                        cw5*final_orbit.Cost1(m);

                    costs(n,2) = cw1*payload_package.Cost2(i) + ...
                        cw2*second_stage.Cost2(j) + ...
                        cw3*lift_off_stage.Cost2(k) + ...
                        cw4*inclination_change.Cost2(l) + ...
                        cw5*final_orbit.Cost2(m);

                    costs(n,3) = cw1*payload_package.Cost3(i) + ...
                        cw2*second_stage.Cost3(j) + ...
                        cw3*lift_off_stage.Cost3(k) + ...
                        cw4*inclination_change.Cost3(l) + ...
                        cw5*final_orbit.Cost3(m);

                    % value
                    values(n,1) = vw1*payload_package.Value1(i) + ...
                        vw2*second_stage.Value1(j) + ...
                        vw3*lift_off_stage.Value1(k) + ...
                        vw4*inclination_change.Value1(l) + ...
                        vw5*final_orbit.Value1(m);

                    values(n,2) = vw1*payload_package.Value2(i) + ...
                        vw2*second_stage.Value2(j) + ...
                        vw3*lift_off_stage.Value2(k) + ...
                        vw4*inclination_change.Value2(l) + ...
                        vw5*final_orbit.Value2(m);

                    values(n,3) = vw1*payload_package.Value3(i) + ...
                        vw2*second_stage.Value3(j) + ...
                        vw3*lift_off_stage.Value3(k) + ...
                        vw4*inclination_change.Value3(l) + ...
                        vw5*final_orbit.Value3(m);


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

tmp = tmp | (dsi(:,2) ~= 1 & (dsi(:,5) == 5 | dsi(:,5) == 6 | dsi(:,5) == 7));

% costs(tmp) = nan;
% values(tmp) = nan;
% reliabilities(tmp) = nan;

%% plot value/cost

x = costs(:,2);
x(~tmp) = nan;

y = values(:,2);
y(~tmp) = nan;

figure(1)

plot(x,y,'k.','MarkerSize',15)
xlabel("Normalized Success Factor",'FontSize',16)
ylabel("Normalized Mission Value",'FontSize',16)
title("Value/Reliability Trade",'FontSize',18)

hold on

%%% look in areas of the plot %%%

% ii = find(y > 1.7 & x < 0.7);
% i = 4;
% plot(x(ii(i)),y(ii(i)),'ro','MarkerSize',10,'LineWidth',2)
% design_space(ii(i))

% xlim([0.85 1.1])
% ylim([0.75 2])

%% look at payloads %%%

% ii = find(dsi(:,1) == 1); % payload == all
% plot(x(ii),y(ii),'r.','MarkerSize',15)
% 
% ii = find(dsi(:,1) == 2); % payload == remote + mag
% plot(x(ii),y(ii),'g.','MarkerSize',15)
% 
% ii = find(dsi(:,1) == 3); % payload == in situ
% plot(x(ii),y(ii),'b.','MarkerSize',15)
% 
% ii = find(dsi(:,1) == 4); % payload == DSI-EUVI-MAG
% plot(x(ii),y(ii),'c.','MarkerSize',15)
% 
% h(1) = plot(nan,nan,'r.');
% h(2) = plot(nan,nan,'g.');
% h(3) = plot(nan,nan,'b.');
% h(4) = plot(nan,nan,'c.');
% legend(h,["All" "Remote + Mag" "In Situ" "DSI-UEVI-MAG"],'Location','best')

%%% look at second stage %%%

ii = find(dsi(:,2) == 1); % second stage == solar sail
plot(x(ii),y(ii),'ro','MarkerSize',6)

ii = find(dsi(:,2) == 2); % second stage == chemical
plot(x(ii),y(ii),'go','MarkerSize',6)

ii = find(dsi(:,2) == 3); % second stage == ion
plot(x(ii),y(ii),'bo','MarkerSize',6)

ii = find(dsi(:,2) == 4); % second stage == nuclear
plot(x(ii),y(ii),'co','MarkerSize',6)

h(1) = plot(nan,nan,'ro');
h(2) = plot(nan,nan,'go');
h(3) = plot(nan,nan,'bo');
h(4) = plot(nan,nan,'co');
legend(h,["Solar Sail" "Chemical" "Ion" "Nuclear"],'Location','best')

%%% look at lift off stage %%%

% ii = find(dsi(:,3) == 1); % first stage = SLS
% plot(x(ii),y(ii),'r^','MarkerSize',8)
% 
% ii = find(dsi(:,3) == 2); % first stage = falcon heavy
% plot(x(ii),y(ii),'g^','MarkerSize',8)
% 
% ii = find(dsi(:,3) == 3); % first stage = sharship
% plot(x(ii),y(ii),'b^','MarkerSize',8)
% 
% ii = find(dsi(:,3) == 4); % first stage = atlas V
% plot(x(ii),y(ii),'c^','MarkerSize',8)
% 
% ii = find(dsi(:,3) == 5); % first stage = delta IV
% plot(x(ii),y(ii),'m^','MarkerSize',8)
% 
% ii = find(dsi(:,3) == 6); % first stage = Vulcan centaur heavy
% plot(x(ii),y(ii),'y^','MarkerSize',8)
% 
% ii = find(dsi(:,3) == 7); % first stage = soyuz 2.1b
% plot(x(ii),y(ii),'k^','MarkerSize',8)
% 
% h(1) = plot(nan,nan,'r^');
% h(2) = plot(nan,nan,'g^');
% h(3) = plot(nan,nan,'b^');
% h(4) = plot(nan,nan,'c^');
% h(5) = plot(nan,nan,'m^');
% h(6) = plot(nan,nan,'y^');
% h(7) = plot(nan,nan,'k^');
% legend(h,["NASA SLS" "SpaceX F-Heavy" "SpaceX Starship" "ULA Atlas V 551" "ULA Delta IV Heavy" "ULA Vulcan Centaur Heavy" "Soyuz 2.1b"],'Location','best')

%%% look at inclination change %%%

% ii = find(dsi(:,4) == 1); % inclination change == Direct
% plot(x(ii),y(ii),'rs','MarkerSize',10,'LineWidth',1.1)
% 
% ii = find(dsi(:,4) == 2); % inclination change == bielliptic
% plot(x(ii),y(ii),'gs','MarkerSize',10,'LineWidth',1.1)
% 
% ii = find(dsi(:,4) == 3); % inclination change == solar sail
% plot(x(ii),y(ii),'bs','MarkerSize',10,'LineWidth',1.1)
% 
% ii = find(dsi(:,4) == 4); % inclination change == gravity inner
% plot(x(ii),y(ii),'cs','MarkerSize',10,'LineWidth',1.1)
% 
% ii = find(dsi(:,4) == 5); % inclination change == gravity outer
% plot(x(ii),y(ii),'ms','MarkerSize',10,'LineWidth',1.1)
% 
% ii = find(dsi(:,4) == 6); % inclination change == venus GA + solar sail
% plot(x(ii),y(ii),'ys','MarkerSize',10,'LineWidth',1.1)
% 
% h(1) = plot(nan,nan,'rs');
% h(2) = plot(nan,nan,'gs');
% h(3) = plot(nan,nan,'bs');
% h(4) = plot(nan,nan,'cs');
% h(5) = plot(nan,nan,'ms');
% h(6) = plot(nan,nan,'ys');
% legend(h,["Direct Transfer from 1 AU" "Bielliptic Transfer" "Solar Sail" "Gravity Assist Inner Planets" "Gravity Assist Outer Planets" "Venus GA + Solar Sail"],'Location','best')

%% look at final orbit %%%

% ii = find(dsi(:,5) == 1); % inclination change == elliptical w/ venus
% plot(x(ii),y(ii),'rp','MarkerSize',8)
% 
% ii = find(dsi(:,5) == 2); % inclination change == elliptical w/ earth
% plot(x(ii),y(ii),'gp','MarkerSize',8)
% 
% ii = find(dsi(:,5) == 3); % inclination change == elliptical w/ mars
% plot(x(ii),y(ii),'bp','MarkerSize',8)
% 
% ii = find(dsi(:,5) == 4); % inclination change == ulysses (jupiter)
% plot(x(ii),y(ii),'cp','MarkerSize',8)
% 
% ii = find(dsi(:,5) == 5); % inclination change == orbit cranking
% plot(x(ii),y(ii),'mp','MarkerSize',8)
% 
% ii = find(dsi(:,5) == 6); % inclination change == non-keplarian
% plot(x(ii),y(ii),'yp','MarkerSize',8)
% 
% ii = find(dsi(:,5) == 7); % inclination change == polaris
% plot(x(ii),y(ii),'kp','MarkerSize',8)
% 
% h(1) = plot(nan,nan,'rp');
% h(2) = plot(nan,nan,'gp');
% h(3) = plot(nan,nan,'bp');
% h(4) = plot(nan,nan,'cp');
% h(5) = plot(nan,nan,'mp');
% h(6) = plot(nan,nan,'yp');
% h(7) = plot(nan,nan,'kp');
% legend(h,["Elliptical Orbit (Venus flyby)" "Elliptical Orbit (Earth flyby)" "Elliptical Orbit (Mars flyby)" "Modified Ulysses orbit (Jupiter flyby)" "Circular orbit (Spiral + Orbit Cranking)" "Non-Keplarian Orbit (beta=0.8)" "POLARIS Final Orbit"],'Location','best')


% plot the reference mission
plot(1,1,'rx','MarkerSize',15,'LineWidth',2)

hold off

grid on
