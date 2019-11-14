function poisSpikesCell = generatePoisson(howManyTrains, howManySpikes, rate)
% generate a sequence of uniformly distributed isi's
unifSpikes = rand(howManyTrains, howManySpikes);
% rate = 1; %function's input
refracTime = 0;
% create a log distribution (of isi's) out of the uniform one
poisSpikes = refracTime - log(1 - unifSpikes)/rate;
% from isi's to spikes (integrate)
poisSpikes = cumsum(poisSpikes, 2);
% from matrix to cell
poisSpikesCell = num2cell(poisSpikes',1)
end