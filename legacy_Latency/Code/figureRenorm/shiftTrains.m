function trains = shiftTrains(trains, latency)

    trains{2} = trains{1} + latency;

end