%% clear
clear
clc

%% data

load tmp.mat

%% no-no's

% "Elliptical Orbit (Venus flyby)" needs "Venus" flyby
selections(selections(:,1) == 1 & selections(:,3) ~= 1) = nan;

% "Elliptical Orbit (Earth flyby)" needs "Earth" flyby
selections(selections(:,1) == 2 & selections(:,3) ~= 4) = nan;

% "Elliptical Orbit (Mars flyby)" needs "Mars" flyby
selections(selections(:,1) == 3 & selections(:,3) ~= 3) = nan;

% "Modified Ulysses orbit (Jupiter flyby)" needs "Jupiter" flyby
selections(selections(:,1) == 4 & selections(:,3) ~= 2) = nan;

% "Circular orbit (Spiral + Orbit Cranking)" needs "solar sail"
selections(selections(:,1) == 5 & selections(:,5) ~= 1) = nan;

% "Non-Keplarian Orbit (beta=0.8)" needs "solar sail"
selections(selections(:,1) == 6 & selections(:,5) ~= 1) = nan;

