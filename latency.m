function b = latency()
% Determing the latency using spike train metrics
% April Fool 2015
% @Nebojsa @Mario @Thomas d_KreuzLab_b

    nStats = 10;
    nSpikes = 1000;

    T = 100;

    %rate is defined as nSpikes/T

    %modify to a point Spike_distance in order to cover the edge
%     Tau = [-T+0.01 : .01 : T-0.01];
    Tau = [-T/2 : .01 : T/2];

    T_d = zeros(1, numel(Tau));

    cases = ['same', 'different', 'same with some noise'];

    doit = cases(3);

    for i = 1 : nStats
        i
       % spikes{1} = poissonSpikeGen(nSpikes/T, T); 
        spikes{1} = poissonGenKreuz(nSpikes, T);

            switch doit

            case cases(1)

                referentTrain = spikes{1};

            case cases(2)

%                referentTrain = poissonSpikeGen(nSpikes/T, T); 
                referentTrain = poissonGenKreuz(nSpikes, T); %

            case cases(3)
                %jitter
                jitter = sort(randn(size(spikes{1})));
                jitter = 0.01*T/nSpikes.*jitter./max(jitter); %30% of fR
                jitter(1) = abs(jitter(1)); jitter(end) = -abs(jitter(end));
                referentTrain = spikes{1} + jitter;
                %false spike detection - deleting
                indx = rand(1, .1*nSpikes);
                indx = indx./max(indx)*nSpikes;
                indx = unique(floor(indx));
                indx = indx(2:end);
                
                referentTrain(indx) = [];
                % parallel coding neuron - adding
                referentTrain = sort([referentTrain, rand(1, .1*nSpikes)]);
                referentTrain = referentTrain - min(referentTrain);
                referentTrain = T.*referentTrain./max(referentTrain);

        end

        S_d = zeros(1, numel(Tau));
        for count = 1 : numel(Tau)
            spikes{2} = referentTrain + Tau(count);
    %         clc 
            (i-1)/nStats*100 + (count-1)/numel(Tau)*100/nStats

            %save test.mat spikes Tau count
            temp =  spkDistance2(spikes, max(0, Tau(count)), min(T, T+Tau(count)), [0, Tau(count)], [T, T+Tau(count)]); 
            %if isnan(temp) || temp > 1 || temp < 0
             %   disp('error') %fix for one point%
            %else
                S_d(count) = temp;
            %end
        end

        T_d = T_d + S_d;
    end

    % hold on,
    figure,
    plot(Tau, T_d./nStats, 'k')
    size(T_d)
    b = mean(T_d./nStats)

end

function spikes = poissonGenKreuz(nSpikes, T)

    rng('shuffle');
    uniform = rand(1,nSpikes);
    spikes = -log(1-uniform)/(nSpikes/T);
    spikes = cumsum(spikes);
    spikes = spikes - min(spikes);
    spikes = T.*spikes./max(spikes);

end

function spikes = poissonSpikeGen(fr, T)

    dt = 1/1000000; % s
    nBins = floor(T/dt);
    spikeMat = rand(1, nBins) < fr*dt;
    tVec = 0:dt:T-dt;

    spikes = tVec(find(spikeMat == 1));

end