function new_powspctrm = zscore_transform(powspctrm)

new_powspctrm = zeros(size(powspctrm));

for electrode = 1:size(powspctrm, 1)
    for frequency = 1:size(powspctrm, 2)
        
        mu = nanmean( powspctrm(electrode, frequency, : ) );
        sigma = nanstd( powspctrm(electrode, frequency, : ) );
        new_powspctrm(electrode, frequency, :) = (powspctrm(electrode, frequency, :) -mu) /sigma;
        
    end
end