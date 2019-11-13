% spike train properties
howManyTrains = 1;
howManySpikes = 100;
lambda = 1;
duration = 99; % T = 100ms;
% latency values
latencyMin = -duration/3;
latencyMax = duration/3;
step = 1; % 1 = 1ms;
%or
latency = [latencyMin : step : latencyMax];
% binning values
howManyBins = 100; % can be a vector also
edges = [latencyMax : (duration + latencyMin - latencyMax)/howManyBins : duration + latencyMin]; %fix this
pearCorr = zeros(size(latency));
spkDist = zeros(size(latency));
vicPurDist = zeros(size(latency));
%some figure properties
% Color %standard
% Linewidth 6
% Font... Lato