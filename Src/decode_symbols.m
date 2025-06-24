function decoded_bits = decode_symbols(symbols, num_bits)
% Description: Decodes complex symbols into binary bits for different 
% modulation schemes based on the decoding regions for each scheme
% Input:
%   - symbols: Received complex symbol array.
%   - num_bits: Number of bits per symbol (1, 2, or 4).
% Output:
%   - decoded_bits: Binary bit array corresponding to input symbols.
    
    switch num_bits
        case 1
            decoded_bits = real(symbols) > 0;
        case 2
            re_bits = real(symbols) > 0;
            im_bits = imag(symbols) < 0;
            decoded_bits = [im_bits, re_bits];
        case 4
            % [bit1 bit2 bit3 bit4]
            bit1 = real(symbols) > 0; % 0xxx or 1xxx
            bit2 = abs(real(symbols)) < 2/sqrt(10); % x0xx or x1xx
            bit3 = imag(symbols) < 0; % xx0x or xx1x
            bit4 = abs(imag(symbols)) < 2/sqrt(10); % xxx0 or xxx1

            decoded_bits = [bit1, bit2, bit3, bit4];
        otherwise
            disp("Only 1, 2, or 4 bits accepted");
    end

end