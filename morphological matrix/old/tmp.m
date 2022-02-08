%% Clear
clear
clc

%% load

load data.mat

%% attitude_actuators

opts = attitude_actuators.Option;
cost = attitude_actuators.Cost;
val = attitude_actuators.Value;

figure(1)
for i = 1:length(opts)
    plot(cost(i),val(i),'.','MarkerSize',10)
    hold on
end
hold off
grid on
xlabel("Cost Function")
ylabel("Value Function")
title("Attitude Actuation System")

xlim([5 85])
ylim([10 70])
text(cost+1,val+1,opts)


%% attitude control

opts = attitude_control.Option;
cost = attitude_control.Cost;
val = attitude_control.Value;

figure(2)
for i = 1:length(opts)
    plot(cost(i),val(i),'.','MarkerSize',10)
    hold on
end
hold off
grid on
xlabel("Cost Function")
ylabel("Value Function")
title("Attitude Control System")

xlim([20 75])
ylim([10 40])
text(cost+1,val+1,opts)

%% science package

opts = science_package.Option1;
cost = science_package.Cost;
val = science_package.Value;

figure(3)
for i = 1:length(opts)
    plot(cost(i),val(i),'.','MarkerSize',10)
    hold on
end
hold off
grid on
xlabel("Cost Function")
ylabel("Value Function")
title("Science Instrament Package")

xlim([0.75 2.5])
ylim([35 250])
text(cost+0.02,val+5,opts)

%% second stage

opts = second_stage.Option;
cost = second_stage.Cost;
val = second_stage.Value;

figure(4)
for i = 1:length(opts)
    plot(cost(i),val(i),'.','MarkerSize',10)
    hold on
end
hold off
grid on
xlabel("Cost Function")
ylabel("Value Function")
title("Second Stage Propulsion")

xlim([2.5 8])
ylim([0 13])
text(cost+0.05,val+0.5,opts)

%% final orbit

opts = final_orbit.Option;
cost = final_orbit.Cost;
val = final_orbit.Value;

figure(5)
for i = 1:length(opts)
    plot(cost(i),val(i),'.','MarkerSize',10)
    hold on
end
hold off
grid on
xlabel("Cost Function")
ylabel("Value Function")
title("Final Operating Orbit")

xlim([0.5 3])
ylim([1 4])
text(cost(1:5)+0.05,val(1:5)+0.08,opts(1:5))
text(cost(end)+0.05,val(end)+0.02,opts(end))

%% flybys

opts = flybys.Option;
cost = flybys.Cost;
val = flybys.Value;

figure(6)
for i = 1:length(opts)
    plot(cost(i),val(i),'.','MarkerSize',10)
    hold on
end
hold off
grid on
xlabel("Cost Function")
ylabel("Value Function")
title("Planetary Flybys")

xlim([0 250])
ylim([0 1100])
text(cost+1,val+50,opts)

%% inclination change

opts = inclin_change.Option;
cost = inclin_change.Cost;
val = inclin_change.Value;

figure(7)
for i = 1:length(opts)
    plot(cost(i),val(i),'.','MarkerSize',10)
    hold on
end
hold off
grid on
xlabel("Cost Function")
ylabel("Value Function")
title("Inclination Change Method")

xlim([20 115])
ylim([5 80])
text(cost+1,val+1,opts)

%% shitty pareto plot


result = zeros(length(inclin_change.Option)*length(flybys.Option)*length(final_orbit.Option)*length(attitude_control.Option)*length(attitude_actuators.Option)*length(second_stage.Option)*length(science_package.Option1),9);

tic
ind = 1;
for i = 1:length(inclin_change.Option)
    for j = 1:length(flybys.Option)
        for k = 1:length(final_orbit.Option)
            for l = 1:length(attitude_control.Option)
                for m = 1:length(attitude_actuators.Option)
                    for n = 1:length(second_stage.Option)
                        for o = 1:length(science_package.Option1)
                            result(ind,1) = ...
                                inclin_change.Cost(i)+...
                                flybys.Cost(j)+...
                                final_orbit.Cost(k)+...
                                attitude_control.Cost(l)+...
                                attitude_actuators.Cost(m)+...
                                second_stage.Cost(n)+...
                                science_package.Cost(o);
                            result(ind,2) = ...
                                inclin_change.Value(i)+...
                                flybys.Value(j)+...
                                final_orbit.Value(k)+...
                                attitude_control.Value(l)+...
                                attitude_actuators.Value(m)+...
                                second_stage.Value(n)+...
                                science_package.Value(o);
                            result(ind,3:end) = [i j k l m n o];
                            ind = ind + 1;
                        end
                    end
                end
            end
        end
    end
end
toc

%% plot 

figure(8)
plot(result(:,1),result(:,2),'x')
grid on
xlabel("Sum of Cost Functions")
ylabel("Sum of Value Functions")
title("Design Space Plot")
