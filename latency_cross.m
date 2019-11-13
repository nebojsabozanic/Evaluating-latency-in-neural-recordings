function b = latency_cross(nBin, noiseLevelj)
% jitter should be from 0 to 1; fool proof later
% Determing the latency using spike train metrics
% April Fool 2015
% @Nebojsa @Mario @Thomas d_KreuzLab_b

    nStats = 10;
    nSpikes = 100;

    T = 100;
    
%     noiseLevelj = 0.1;
    noiseLeveld = 0; noiseLevela = 0;

    %rate is defined as nSpikes/T

    %modify to a point Spike_distance in order to cover the edge
%     Tau = [-T+0.01 : .01 : T-0.01];
    tauStep = T/nBin;
    Tau = [-T/2 : tauStep : T/2];

    T_d = zeros(1, numel(Tau));
    Tc_d = zeros(1, numel(Tau));
    Tv_d = zeros(1, numel(Tau));
    
    cases = ['same', 'different', 'same with some noise'];

    doit = cases(3);

    for i = 1 : nStats
%         i
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
                if (noiseLevelj)
                    jitter = noiseLevelj*T/nSpikes*jitter./max(jitter); %30% of fR
                else
                    jitter = zeros(size(jitter));
                end
                jitter(1) = abs(jitter(1)); jitter(end) = -abs(jitter(end));
                referentTrain = spikes{1} + jitter;
                %false spike detection - deleting
%                 indx = rand(1, noiseLeveld*nSpikes);
%                 indx = indx./max(indx)*nSpikes;
%                 indx = unique(floor(indx));
%                 indx = indx(2:end);
%                 
%                 referentTrain(indx) = [];
                % parallel coding neuron - adding
%                 add2referentTrain = rand(1, noiseLevela*nSpikes);
%                 add2referentTrain = poissonGenKreuz(noiseLevela*nSpikes, T);
% %                 add2referentTrain = add2referentTrain - min(add2referentTrain);%
% %                 add2referentTrain = add2referentTrain/max(add2referentTrain)*T;
%                 add2referentTrain = add2referentTrain(2:end-1);
%                 referentTrain = sort([referentTrain, add2referentTrain]);
%                 

%                 referentTrain = sort([referentTrain, rand(1, noiseLevela*nSpikes)]);
%                 referentTrain = referentTrain - min(referentTrain);
%                 referentTrain = T.*referentTrain./max(referentTrain);
        end

        C_c = zeros(1, numel(Tau));
        S_d = zeros(1, numel(Tau));
        V_c = zeros(1, numel(Tau));
        for count = 1 : numel(Tau)
            spikes{2} = referentTrain + Tau(count);
          %   clc 
            %(i-1)/nStats*100 + (count-1)/numel(Tau)*100/nStats

            %save test.mat spikes Tau count
            temp =  spkDistance2(spikes, max(0, Tau(count)), min(T, T+Tau(count)), [0, Tau(count)], [T, T+Tau(count)]);
            tauStep = (min(T, T+Tau(count)) - max(0, Tau(count)))/nBin;
            psth{1} = histc(spikes{1}, max(0, Tau(count)) : tauStep : min(T, T+Tau(count)));
            psth{2} = histc(spikes{2}, max(0, Tau(count)) : tauStep : min(T, T+Tau(count))); %tauStep
            R = corrcoef(psth{1}(1:end-1), psth{2}(1:end-1));
            vc = 0; %SPIKY_Victor_MEX(single(spikes{1}), single(spikes{2}), single(10)); %put as an argument 100
            
            
            C_c(count) = R(1,2);
            S_d(count) = temp;
            V_c(count) = vc;
            
        end
    
        T_d = T_d + S_d;
        Tc_d = Tc_d + C_c;
        Tv_d = Tv_d + V_c;
    end

    % hold on,
    %figure,
    normT_d = (T_d./nStats - mean(T_d./nStats))./std(T_d./nStats);
    normTc_d = -(Tc_d./nStats - mean(Tc_d./nStats))./std(Tc_d./nStats);
    normTv_d = (Tv_d./nStats - mean(Tv_d./nStats))./std(Tv_d./nStats);
    %plot(Tau, normT_d, 'k', Tau, normTc_d, 'm')   
%     title(str2num(nSpikes/tauStep));
    %size(T_d);
    %b = [mean(T_d./nStats) mean(Tc_d./nStats); std(T_d./nStats) std(Tc_d./nStats)] 
    b(1) = min(normT_d);
    b(2) = min(normTc_d); %change 501 into the half of all
    b(3) = min(normTv_d);
    
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