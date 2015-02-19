%% Time-Freq Analysis
cfg.output       = 'pow';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;
freq = ft_freqanalysis(cfg, dataSSVEP);

%% Multiplot
cfg.baseline     = [-1 0]; 
cfg.baselinetype = 'relative'; 
cfg.zlim         = [0 15];	        
cfg.channel = freq.label(1:64);
cfg.showlabels   = 'yes';
cfg.layout       = 'biosemi64.lay';
figure;
ft_multiplotTFR(cfg, freq);