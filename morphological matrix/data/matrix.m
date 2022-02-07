%% clear

clear
clc

%% data

load data.mat

%% create permutation matrix

A = length(final_orbit.opts);
B = length(first.opts);
C = length(flybys.opts);
D = length(inclin.opts);
E = length(kicker.opts);
F = length(payload.opts);

N = A*B*C*D*E*F;

tmp(N) = struct();

%% do the loop

tic
n = 1;
for a = 1:A
    for b = 1:B
        for c = 1:C
            for d = 1:D
                for e = 1:E
                    for f = 1:F

                        costs = (final_orbit.cost(a) + ...
                            first.cost(b) + ...
                            flybys.cost(c) + ...
                            inclin.cost(d) + ...
                            kicker.cost(e) + ...
                            payload.cost(f))/6;

                        values = (final_orbit.value(a) + ...
                            first.value(b) + ...
                            flybys.value(c) + ...
                            inclin.value(d) + ...
                            kicker.value(e) + ...
                            payload.value(f))/6;

                        reliabilities = (final_orbit.reliability(a) + ...
                            first.reliability(b) + ...
                            flybys.reliability(c) + ...
                            inclin.reliability(d) + ...
                            kicker.reliability(e) + ...
                            payload.reliability(f))/6;

                        tmp(n).final_orbit = final_orbit.opts(a);
                        tmp(n).first = first.opts(b);
                        tmp(n).flybys = flybys.opts(c);
                        tmp(n).inclin = inclin.opts(d);
                        tmp(n).kicker = kicker.opts(e);
                        tmp(n).payload = payload.opts(f);

                        tmp(n).cost = costs;
                        tmp(n).value = values;
                        tmp(n).reliability = reliabilities;

                        n = n + 1;
                    end
                end
            end
        end
    end
end
toc

%% clean ups

costs = zeros(N,1);
values = zeros(N,1);
reliabilities = zeros(N,1);

tic
for n = 1:N
    [costs(n), values(n), reliabilities(n)] = calc_results(tmp(n));
end
toc

%% plot

figure(1)
plot(costs,values,'.')
grid on
xlabel("Relative Cost",'FontSize',16)
ylabel("Relative Mission Value",'FontSize',16)
title("Mission Value/Cost Trade",'FontSize',20)

figure(2)
plot(costs,reliabilities,'.')
grid on
xlabel("Relative Cost",'FontSize',16)
ylabel("Relative Reliability",'FontSize',16)
title("Reliability/Cost Trade",'FontSize',20)

figure(3)
plot(reliabilities,values,'.')
grid on
xlabel("Relative Reliability",'FontSize',16)
ylabel("Relative Mission Value",'FontSize',16)
title("Mission Value/Reliability Trade",'FontSize',20)


% %%
%
% [v,i] = max(values./reliabilities,[],'all');
%
% get_opts(i,selections,final_orbit.opts,first.opts,flybys.opts,inclin.opts,kicker.opts,payload.opts)
%
% hold on
% plot(reliabilities(i),values(i),'rx','MarkerSize',10)


%% get options function
function [value,cost,reliability] = calc_results(data)

value = nan;
cost = nan;
reliability = nan;

if data.final_orbit == "Elliptical Orbit (Venus flyby)" && ...
        data.flybys ~= "Venus"
    return
end
if data.final_orbit == "Elliptical Orbit (Earth flyby)" && ...
        data.flybys ~= "Earth"
    return
end
if data.final_orbit == "Elliptical Orbit (Mars flyby)" && ...
        data.flybys ~= "Mars"
    return
end
if data.final_orbit == "Modified Ulysses orbit (Jupiter flyby)" && ...
        data.flybys ~= "Jupiter"
    return
end
if data.final_orbit == "Circular orbit (Spiral + Orbit Cranking)" && ...
        data.kicker ~= "Solar Sail"
    return
end
if data.final_orbit == "Non-Keplarian Orbit (beta=0.8)" && ...
        data.kicker ~= "Solar Sail"
    return
end
if data.inclin == "Solar Sail" && data.kicker ~= "Solar Sail"
    return
end
if data.inclin == "Direct Transfer from 1 AU" && data.kicker == "Solar Sail"
    return
end
if data.inclin == "Bielliptic Transfer" && data.kicker == "Solar Sail"
    return
end
if data.inclin == "Gravity Assist Inner Planets" && data.flybys == "Jupiter"
    return
end
if data.inclin == "Gravity Assist Outer Planets" && data.flybys ~= "Jupiter"
    return
end
if data.inclin == "Venus GA + Solar Sail" && (data.flybys ~= "Venus" || data.kicker ~= "Solar Sail")
    return
end

value = data.value;
cost = data.cost;
reliability = data.reliability;

end