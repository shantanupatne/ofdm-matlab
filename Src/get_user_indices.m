function indices = get_user_indices(N_sc, num_frames, pilot, repition)
    indices = cell(num_frames, 2);
    idx1 = 1:N_sc/2;
    idx2 = N_sc/2 + 1:N_sc;
    if ~pilot && ~repition
        indices{1, 1} = idx1;
        indices{1, 2} = idx2;
        indices{2, 1} = idx2;
        indices{2, 2} = idx1;
        indices{3, 1} = idx1;
        indices{3, 2} = idx2;
    elseif ~repition
        indices{1, 1} = [2, 3, 4, 6, 7, 8, 10];
        indices{1, 2} = [11, 12, 14, 15, 16, 18, 19];
        indices{2, 1} = idx2;
        indices{2, 2} = idx1;
        indices{3, 1} = idx1;
        indices{3, 2} = idx2;
    else
        indices{1, 1} = [2, 3, 4, 10, 11, 12, 18, 19, 20];
        indices{1, 2} = [6, 7, 8, 14, 15, 16];

        indices{2, 1} = [1, 2, 3, 4, 9, 10, 11, 12, 17, 18, 19, 20];
        indices{2, 2} = setdiff(1:20, indices{2, 1});

        indices{3, 1} = [5, 6, 7, 8, 13, 14, 15, 16];
        indices{3, 2} = setdiff(1:20, indices{3, 1});
    end
end

