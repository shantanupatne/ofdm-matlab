function [X_tx, bitstring1, bitstring2] = gen_multi_user_repetition(N_sc, L, num_frames, syms, user_indices, pilot_indices, IDFT_matrix)
% GEN_MULTI_USER_REPETITION
% Generates an OFDM signal for multiple users with pilot symbols 
% inserted in the first frame. Uses repetition coding to increase
% diversity.
% Frame 1: 3 Symbols each of UE1 and UE2 with Pilots in between
%           UE1 symbols repeat twice alternating with UE2 symbols 
%           [P, U1S1, U1S2, U1S3, P, U2S2, U2S2, U2S3, P, U1S1, U1S2, U1S3,
%           P, U2S1, U2S2, U3S3, P, U1S1, U1S2, U1S3]
% Frame 2: 8 symbols of UE1 and 4 symbols of UE2. ( no pilots )
%          First 4 symbols of each user repeat twice, alternating between each user
%          Last 4 symbols of UE1 repeat in Frame 3.
% Frame 3: 4 symbols of UE1 and 8 symbols of UE2 (no pilots)
%
% INPUTS:
%   N_sc        - Number of subcarriers in an OFDM frame
%   L           - Number of channel taps
%   num_frames  - Number of OFDM frames
%   syms        - Modulation symbol mapping
%   num_users   - Number of users
%   IDFT_matrix - IDFT matrix for transformation to time domain
%
% OUTPUTS:
%   X_tx         - Serialized OFDM signal with all frames
%   bitstring    - Random bitstring transmitted for each user
%   pilot        - Pilot symbol used
%   user_indices - Cell array containing subcarrier indices allocated to
%                  each user in each frame
%
    num_bits = log2(length(syms));

    cp_length = L - 1; % Cyclic prefix length
    pilot = 1; % Pilot symbol (usually average energy of constellation)

    % Generate random bit sequence for all users
    data_bits = num_bits * 15;
    bitstring1 = randi([0, 1], data_bits, 1);
    bitstring2 = randi([0, 1], data_bits, 1);

    rep_bitstring1 = [repmat(bitstring1(1:3*num_bits), 3, 1); repmat(bitstring1(3*num_bits + 1 : 7*num_bits), 2, 1); repmat(bitstring1(7*num_bits + 1:11*num_bits), 2, 1); bitstring1(11*num_bits + 1:end)];
    rep_bitstring2 = [repmat(bitstring2(1:3*num_bits), 2, 1); repmat(bitstring2(3*num_bits + 1 : 7*num_bits), 2, 1); repmat(bitstring2(7*num_bits + 1:11*num_bits), 2, 1); bitstring2(11*num_bits + 1:end)];
    
    % Initialize output signal
    X_tx = zeros((N_sc + cp_length) * num_frames, 1);
    
    bit_offset = [0, 0];
    for j = 1:num_frames
        X_freq = zeros(N_sc, 1);

        if j == 1
            X_freq(pilot_indices) = pilot;
        end

        [user1_idx, user2_idx] = user_indices{j, :};

        user1_bits = rep_bitstring1(bit_offset(1) + 1: bit_offset(1) + length(user1_idx) * num_bits);
        user1_bits = reshape(user1_bits, [num_bits, length(user1_idx)]);

        user2_bits = rep_bitstring2(bit_offset(2) + 1: bit_offset(2) + length(user2_idx) * num_bits);
        user2_bits = reshape(user2_bits, [num_bits, length(user2_idx)]);

        u1_syms = syms(bit2int(user1_bits, num_bits) + 1);
        u2_syms = syms(bit2int(user2_bits, num_bits) + 1);

        X_freq(user1_idx) = u1_syms;
        X_freq(user2_idx) = u2_syms;

        bit_offset(1) = bit_offset(1) + length(user1_idx)*num_bits;
        bit_offset(2) = bit_offset(2) + length(user2_idx)*num_bits;

        X_time = IDFT_matrix * X_freq;

        X_cp = [X_time(end - cp_length + 1:end); X_time];
        X_tx((N_sc + cp_length) * (j - 1) + 1 : (N_sc + cp_length) * j) = X_cp;

    end
end

