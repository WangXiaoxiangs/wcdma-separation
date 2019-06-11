% 28kbps simple WCDMA uplink
% include:
%   CRC, padding, convolution encode, spread spectrum with OVSF,
%   scrambling with Gold, upsampling for 8 times
% paramter:
%   framelen:   280
%   CRC:        CRC - 16, default
%   padding:    4 bits
%   convcode:   len = 9, [561, 753];
%   Gold:       2, with different shift
%   ovsf:       idx = 3
%   spreading factor: 64
%   sRRC:       8 times, rolloff = 0.22
%   one frame, 15 slots, 10 ms, single user, no channel
%   DPCCH * 1 + DPDCH * 1
% input: 
%   user:       1 * 1 double, integer number, for different users
% output: 307200 * 1 complex double
%   (280 + 16 + 4) * 2 * 64 * 8 = 307200;

function [y, data, scrambling_code] = simpleWCDMAgenerate(user)
    % CRC parameter
    crcGen = comm.CRCGenerator('z^16 + z^15 + z^2 + 1');
    % convlution code
    trellis = poly2trellis(9, [561 753]);
    % ovsf code
    ovsf_code = comm.OVSFCode('SpreadingFactor', 64, 'Index', 3, ...
        'SamplesPerFrame', 2560);
    % Gold code
    gold_i = comm.GoldSequence('FirstPolynomial', [25,3,0], ...
        'FirstInitialConditions', [zeros(1, 24), 1], ...
        'SecondPolynomial', [25,3,2,1,0], ...
        'SecondInitialConditions', [zeros(1, 24), 1], ...
        'Index', user, 'Shift', 0, 'SamplesPerFrame', 2560);
    gold_q = comm.GoldSequence('FirstPolynomial', [25,3,0], ...
        'FirstInitialConditions', [zeros(1, 24), 1], ...
        'SecondPolynomial', [25,3,2,1,0], ...
        'SecondInitialConditions', [zeros(1, 24), 1], ...
        'Index', user, 'Shift', 16777232, 'SamplesPerFrame', 2560);
    % channel
    sRRC_uplink = comm.RaisedCosineTransmitFilter(...
        'Shape', 'Normal', 'RolloffFactor', 0.22, ...
        'FilterSpanInSymbols', 10, 'OutputSamplesPerSymbol', ...
        8, 'Gain', 1/sqrt(8));

    % DPDCH data
    data = randi([0 1], 280, 1);            % in 10 ms, 15 slots
    data(1: 16) = ones(16, 1);
    data_crc = [crcGen(data); zeros(4, 1)]; % apply CRC
    data_enc = convenc(data_crc, trellis);  % apply convlution code
    ovsf_sf = ovsf_code();                  % ovsf code
    dpdch_I = zeros(2560, 15);              % dpdch data
    % split into 15 slots
    for k = 1: 1: 15
        dpdch_code = data_enc((k - 1) * 40 + (1: 40));  % each slot
        dpdch_bipolar = pskmod(dpdch_code, 2);          % to Bipolar
        dpdch_sf = reshape(repmat(dpdch_bipolar, 1, 64).', 2560, 1);
        dpdch_I(:, k) = dpdch_sf .* ovsf_sf;            % generate dpdch
    end
    dpdch_I = dpdch_I(:);
    % DPCCH data
    pilot = [1; 1; 1; 1; 1; 1; 0; 0; 0; 0]';            % pilot, random
    sf_code = ones(2560, 1);
    pilot_bipolar = pskmod(pilot, 2);                   % to bipolar
    pilot_sf = repmat(pilot_bipolar, 256, 1);           % repmat the pilot
    dpcch_Q = repmat(pilot_sf(:) .* sf_code, 15, 1);    % dpcch data
    % scrambling code
    scrambling_code = zeros(2560, 15);
    for k = 1: 1: 15
        gold_I = gold_i();
        gold_Q = gold_q();
        I_ch = pskmod(gold_I, 2);
        Q_ch = pskmod(gold_Q, 2);
        scrambling_code(:, k) = I_ch + 1i * Q_ch;
    end
    scrambling_code = scrambling_code(:) / sqrt(2);
    % DPCH generate
    dpch_iq = (dpdch_I + 1i * dpcch_Q) .* scrambling_code / sqrt(2);
    % apply sRRC and channel
    y = sRRC_uplink(dpch_iq);
end