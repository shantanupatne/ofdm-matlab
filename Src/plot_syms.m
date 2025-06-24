function plot_syms(syms, plot_title)
    scatterplot(syms, 1, 0, 'b*')
    for k = 1:length(syms)
        text(real(syms(k)) - 0.1,imag(syms(k)), ...
            dec2base(k - 1,2,4),'Color',[0 1 0]);
    end
    title(plot_title)
    axis([-1 1 -1 1])
end
