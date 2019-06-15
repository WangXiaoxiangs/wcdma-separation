% GDE algorithm and AIC\MDL algorithm
%% parameters
N = 4;          % number of sources
M = 12;          % received attennas
num = 1e4;      % singal samples

mimochannel = comm.MIMOChannel(...
    'SampleRate', 1000, ...
    'PathDelays', zeros(1, N), ...
    'AveragePathGains', zeros(1, N), ...
    'MaximumDopplerShift', 0, ...
    'SpatialCorrelationSpecification', 'None', ...
    'NumTransmitAntennas', N, ...
    'NumReceiveAntennas', M, ...
    'PathGainsOutputPort', true);

%% signal generate
data = randi([0 3], num, N);        % suppose QPSK modulation
datamod = pskmod(data, 4, pi/4);    % QPSK modulation
txsig = awgn(datamod, 20, 'measured');   % add noise  
[rxsig, H] = mimochannel(txsig);    % after channel

%% GDE algorithm
Rx = rxsig' * rxsig / num;          % autocorrlation
Rx_part = Rx((1: M-1), (1: M-1));   % M - 1\ order
[V, D] = eig(Rx_part);              % eigenvalue decomp.
[~, idx] = sort(diag(D), 'descend');% sort eigenvalue
V1 = V(:, idx);

GammaMartix = [V1, zeros(M - 1, 1); ...
    zeros(1, M - 1), 1];            % Unity martix
T = GammaMartix' * Rx * GammaMartix;% after transform
% GerschgorinPlot(T);                 % Gerschgorin graph

GDE = zeros(M - 2, 1);
temp = sum(abs(T((1: M - 1), end))) * 0.5 / (M - 1);
for k = 1: 1: M - 2
    GDE(k) = abs(T(k, end)) - temp;
end
number_source_GDE = find(GDE < 0, 1) - 1;

%% AIC and MDL algorithm
ev = sort(eig(Rx), 'descend');      % front part is signal subspace
aic_rx = zeros(M, 1);               % AIC criterion
mdl_rx = zeros(M, 1);               % MDL criterion
for k = 1: M                        % AIC and MDL calc
    nec = mean(ev(k: M));           % get noise subspace, calc noise power
    tmp_snr_est = -num * sum(log(ev(k: M))) + (M - k + 1) * num * log(nec);
    aic_rx(k) = tmp_snr_est + 2 * (k - 1) * (2 * M - k + 1);
    mdl_rx(k) = tmp_snr_est + k / 2 * (2 * M - k + 1) * log(num);
end
[~, K_idx_aic] = min(real(aic_rx));
N_est_aic = K_idx_aic - 1;
[~, K_idx_mdl] = min(real(mdl_rx));
N_est_mdl = K_idx_mdl - 1;