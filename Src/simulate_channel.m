function [y_noisy, h] = simulate_channel(x, N_sc, L, SNR_dB)
% SIMULATE_CHANNEL
% Description: Simulates a multipath fading channel with additive white 
% Gaussian noise (AWGN). Computes the channel impulse response, generates 
% noisy channel output, and calculates noise power.
% Input:
%   - x: Input signal vector.
%   - N_sc: Number of subcarriers in each frame.
%   - L: Length of the channel impulse response (number of taps).
%   - SNR_dB: Signal-to-Noise Ratio in decibels.
% Output:
%   - y_noisy: Output signal after passing through channel and adding noise.
%   - h: Channel impulse response.    
    
    % Simulate channel with random taps
    h = (randn(L, 1) + 1i .* randn(L, 1)) .* sqrt(1/(2*L));        
    y = convolution(h, x);

    % Get Expected Symbol Energy (approx. 1)
    Es = 0;
    for frame = 1:3
        start_idx = (N_sc + L - 1)*(frame - 1) + (L - 1);
        end_idx = frame *(N_sc + L - 1);
        Es = Es + sum(abs(x(start_idx:end_idx)) .^ 2) / N_sc; % Symbol energy
    end
    Es = Es / 3;

    % Noise power
    SNR_linear = 10^(SNR_dB/10);
    N0 = Es / SNR_linear;

    % Generate Gaussian noise
    noise = sqrt(N0/2) .* (randn(size(y)) + 1i .* randn(size(y)));
    
    % Add noise to simulated channel output
    y_noisy = y + noise;
end

