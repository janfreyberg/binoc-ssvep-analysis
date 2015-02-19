function powspctrm = maxtransform(powspctrm)

for electrode = 1:size(powspctrm, 1)
    for frequency = 1:size(powspctrm, 2)
        
        powspctrm(electrode, frequency, :) = powspctrm(electrode, frequency, :)/max(powspctrm(electrode, frequency, :));
        
    end
end