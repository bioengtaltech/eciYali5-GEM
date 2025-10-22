function [kcat] = calcKcatFromActivity(ecModel,enzymeUNIPROT,specificActivity)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
enzMW = ecModel.ec.mw(strcmp(ecModel.ec.enzymes,enzymeUNIPROT)); % Get MW of the enzyme.
kcat = specificActivity; % umol/min/mg protein, same as mmol/min/g protein.
kcat = kcat / 1000; % mol/min/g protein.
kcat = kcat / 60; % mol/sec/g protein.
kcat = kcat * enzMW;
end