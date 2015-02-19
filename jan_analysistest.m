clear all

cfg = [];

% preprocessing params
cfg.continuous = 'yes';
cfg.demean    = 'yes';
cfg.detrend = 'yes';
cfg.lpfreq = 49;
cfg.reref = 'yes';
cfg.refchannel = 'all';

% Trial definition
cfg.dataset = 'E:\Dropbox\Documents\University\Matlab Code\EEG Study 2015\Jan_MSSVEP.bdf';
cfg.trialdef.eventtype = 'STATUS';
cfg.trialdef.prestim = 2;
cfg.trialdef.poststim = 10;
cfg.channel = 1:64;

cfg.toi          = -0.5:0.05:10;

cfg1 = cfg;
cfg2 = cfg;
cfg1.trialdef.eventvalue = [11 12];
cfg2.trialdef.eventvalue = [21 22];

% TFR params
cfg1.foi = 20:0.4:45; % frequencies of interest
cfg2.foi = 20:0.4:45;

cfg1 = ft_definetrial(cfg1);
dataSSVEP = ft_preprocessing(cfg1);
cfg = cfg1;

jan_tfr;
cfg.toi = 0:0.05:8;
jan_fourier_transformdata;

cfg2 = ft_definetrial(cfg2);
dataSSVEP = ft_preprocessing(cfg2);
cfg = cfg2;

jan_tfr;
cfg.toi = 0:0.05:8;
jan_fourier_transformdata;