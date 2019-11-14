function binnedTrains = binning(trains, edges)
for cnt = 1 : numel(trains)
    binnedTrains{cnt} = histcounts(trains{cnt}, edges); 
end
end