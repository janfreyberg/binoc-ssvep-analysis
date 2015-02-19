%% make spectrum
cfg.output       = 'pow';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;
freq = ft_freqanalysis(cfg, dataSSVEP);

figure;
plot(freq.freq, freq.powspctrm);
legend(freq.label)
ylim([0 10]);