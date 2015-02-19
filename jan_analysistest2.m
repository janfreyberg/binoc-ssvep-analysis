clear all

cfg = [];
cfg.dataset = fullfile(pwd, 'EEG Feb 2015', 'EG-CTR-0005-MSSVEP.bdf');
cfg.trialdef.eventtype = 'STATUS';

% cfg.trialdef.eventvalue = [21, 22];
cfg.trialdef.eventvalue = 1:255;

cfg.trialdef.prestim = 1;
cfg.trialdef.poststim = 9;
cfg.channel = 1:64;
cfg.continuous = 'yes';
cfg.demean    = 'yes';
cfg.detrend = 'yes';
cfg.lpfreq = 49;
cfg.hpfreq = 20;
cfg.reref = 'yes';
cfg.refchannel = 'all';


cfg = ft_definetrial(cfg);

% remove any mistaken trials
kickOut = zeros(size(cfg.trl));
if size(cfg.trl, 1) > 1
    kickOut(1, 1) = false;
    for i = 2:size(cfg.trl, 1)
        if cfg.trl(i, 1) <= cfg.trl(i-1, 1) + 8
            kickOut(i, 1) = true;
        else
            kickOut(i, 1) = false;
        end
    end
    cfg.trl(logical(kickOut), :) = [];
end


loaded_data = ft_preprocessing(cfg);


cfg.foi = 24.8:0.4:40;
cfg.toi = -2:0.05:15.5;

cfg.output       = 'pow';
cfg.method = 'wavelet';
cfg.width = 30;
% cfg.taper        = 'hanning';
freq = ft_freqanalysis(cfg, loaded_data);


% Make the TFR Plot
cfg.baseline     = [-1 0];
cfg.baselinetype = 'relchange';
% cfg.zlim         = [0 15];
cfg.channel = freq.label(1:64);
cfg.showlabels   = 'yes';
cfg.layout       = 'biosemi64.lay';
figure;
ft_multiplotTFR(cfg, freq);
% ft_singleplotTFR(cfg, freq);