function latency_batch

% nBin = [50:50:500];%10.^[1:3];
nBin = 100;
jitter = [0, 0.05, 0.1, 0.2, 0.4, 0.8];
spkDistv = zeros(size(jitter));
crosscorv = zeros(size(jitter));
vp = zeros(size(jitter));

for count = 1 : numel(jitter);
    count
    
    b = latency_cross(nBin, jitter(count));
    spkDistv(count) = b(1);
    crosscorv(count) = b(2);
    vp(count) = b(3);
end

plot(jitter, spkDistv, jitter, crosscorv, jitter, vp, 'LineWidth',6)
legend('SPIKE-distance amplitude of the peak', 'Cross-correlation amplitude of the peak', 'Victor purpura aplitude of the peak');