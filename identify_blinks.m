cfg = [];


% basic preprocessing
cfg.continuous = 'yes';

cfg.demean    = 'no';
cfg.detrend = 'no';

cfg.reref = '1:64';
cfg.refchannel = 'all';

cfg.dataset = fullfile(pwd, 'EEG Feb 2015', 'EG-CTR-0008-BinSSVEP.bdf');
cfg.channel = 1:72;

loaded_data = ft_preprocessing(cfg);

% resample
cfg.resamplefs = 300;
