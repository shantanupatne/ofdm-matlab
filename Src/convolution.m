function y = convolution(h, X)
    y = zeros(size(X));    
    for n = 1:length(y)
        for p = 1:length(h)
            if (n - p + 1 > 0) && (n - p + 1 <= length(X))
                y(n) = y(n) + X(n - p + 1) * h(p);
            end
        end
    end  
end