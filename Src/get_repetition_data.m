function rep_data = get_repetition_data()
% [
%       repeated sym in prev frame, 
%       num repeated symbols in this frame,
%       num unique symbols
%       repetition_offset
% ]
    rep_data = cell(3, 2);
    rep_data{1, 1} = [0, 3, 3, 0];
    rep_data{2, 1} = [0, 2, 8, 4];
    rep_data{3, 1} = [1, 1, 8, 4];

    rep_data{1, 2} = [0, 2, 3, 0];
    rep_data{2, 2} = [0, 2, 4, 0];
    rep_data{3, 2} = [0, 2, 8, 4];
end

