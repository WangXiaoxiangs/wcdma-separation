% GDE algorithm and AIC\MDL algorithm
%% parameters
N = 4;          % number of sources
M = 12;          % received attennas
num = 1e4;      % singal samples
upsamples = 8;  % upsample times
mimochannel = comm.MIMOChannel(...
    'SampleRate', 1000, ...
    'PathDelays', zeros(1, N), ...
    'AveragePathGains', zeros(1, N), ...
    'MaximumDopplerShift', 0, ...
    'SpatialCorrelationSpecification', 'None', ...
    'NumTransmitAntennas', N, ...
    'NumReceiveAntennas', M, ...
    'PathGainsOutputPort', true);
sRRC_filt = comm.RaisedCosineTransmitFilter(...
    'Shape', 'Square root', ...
    'RolloffFactor', 0.3, ...
    'FilterSpanInSymbols', 10, ...
    'OutputSamplesPerSymbol', upsamples);

%% signal generate
data = randi([0 3], num, N);        % suppose QPSK modulation
datamod = pskmod(data, 4, pi/4);    % QPSK modulation
dataRRc = sRRC_filt(datamod);       % sRRC filter
[txsig, H] = mimochannel(dataRRc);  % after MIMO channel
noise_Var = 10 + 10 * rand(M, 1);   % add different power noise
rxsig = zeros(upsamples * num, M);  % received signal
for k = 1: 1: M
    rxsig(:, k) = awgn(txsig(:, k), noise_Var(k), 'measured');
end

%% IDGE, detect more sources
% step 1, calc different Rx(\tau)
tau = 10: 2: 40;
len = length(tau);
Zx  = zeros(M ^ 2, len);
Rx_temp = zeros(M, M, len);
for k = 1: 1: len
    % calc a.c.fcn
    Rx_temp(:, :, k) = rxsig' * rxsig / (num * upsamples);
    % make a vector
    Zx(:, k) = reshape(Rx_temp(:, :, k), M ^ 2, 1);
end
% step 2, vector and calc its a.c.fcn
Rz = Zx * Zx' / len;
% step 3, eigenvalue decomp.
gz = 2 * M;
Rz_part = Rz((1: gz), (1: gz));
[V, D] = eig(Rz_part);

ksai = V' * Rz((1: gz), (gz + 1: end));
% step 4, calc radius