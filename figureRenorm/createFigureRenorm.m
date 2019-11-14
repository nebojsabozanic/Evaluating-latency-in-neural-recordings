% initialize
ini
% generate two Poisson spike trains
trains = generatePoisson(howManyTrains, howManySpikes, lambda);
% in this figure we are comparing two identical Poisson spike trains
trains{2} = trains{1};
for indx = 1 : numel(latency)
    trains = shiftTrains(trains, latency(indx));
    % bin it in order to calculate correlation % or % kernelize it + smooth
    binnedTrains = binning(trains, edges);
    % calculate Pearson correlation
    pearCorr_temp = corrcoef(binnedTrains{1}, binnedTrains{2}); 
    pearCorr(indx) = pearCorr_temp(1,2); % since corrcoef returns [2x2]...
    % calculate SPIKE-distance
    spkDist(indx) = spkDist_mex(trains, max(0, latency(indx)), min(duration, duration+latency(indx)), [0, latency(indx)], [duration, duration+latency(indx)]);
    % calculate Victor-Purpura distance
    vicPurDist(indx) = SPIKY_Victor_MEX(trains{1}, trains{2}, 10);
end
% plot original
figure, plot(latency, pearCorr, latency, spkDist, latency, vicPurDist, 'LineWidth', 4); %6
legend('Pearson correlation', 'SPIKE-distance', 'Victor-Purpura metric')
% plot normalized
figure, plot(latency, normalize(pearCorr), latency, normalize(spkDist), latency, normalize(vicPurDist), 'LineWidth', 4); %6
legend('Pearson correlation', 'SPIKE-distance', 'Victor-Purpura metric')
