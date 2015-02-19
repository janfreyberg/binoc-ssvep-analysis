cfg = [];
cfg.continuous = 'yes';
cfg.demean    = 'yes';
cfg.detrend = 'yes';
% cfg.lpfreq = 49;
% cfg.hpfreq = 20;
cfg.reref = '1:64';
cfg.refchannel = 'all';

%% Trial definition
fileID = 'EG-CTR-0008-BinSSVEP.bdf';
cfg.dataset = fullfile('E:\Documents\Recorded Data\EEG Feb 2015', fileID);
cfg.channel = 1:64;

cfg.trialdef.eventtype = 'STATUS';
cfg.trialfun = 'ft_trialfun_general';
cfg.trialdef.prestim = 2;
cfg.trialdef.poststim = 15.5;

trial_data = cell(1, 16);
for eventvalue = 1:16;
    
    cfg.trialdef.eventvalue = 200 + eventvalue;
    
    cfg.trl = [];
    cfg = rmfield(cfg, 'trl');
    cfg = ft_definetrial(cfg);
    
    cfg.trl = remove_overlaps(cfg.trl); % this scripts removes overlapping trials
    
    trial_data{eventvalue} = ft_preprocessing(cfg);
    
end

%% TFR Parameters
cfg = [];
cfg.toi = -2:0.05:15.5;
% cfg.foi = 24.8:0.4:40;
cfg.foi = [28.8, 36];

cfg.output = 'pow';
cfg.method = 'wavelet'; % other option would be multi-taper-convolution ('mtmconvol')
cfg.width = 20; % width of the wavelet window in terms of wave cycles

freq = cell(1, 16);

for eventvalue = 1:16
    
    freq{eventvalue} = ft_freqanalysis(cfg, trial_data{eventvalue});
    
end


%% Scaling
relfreq = cell(1, 16);
zscfreq = cell(1, 16);
maxfreq = cell(1, 16);
for eventvalue = 1:16
    
    index28 = find( freq{eventvalue}.freq == 28.8 );
    index36 = find( freq{eventvalue}.freq == 36 );
    
    % baseline correction
    cfg = [];
    cfg.baseline = [-2 0];
    cfg.baselinetype = 'absolute';
    
    relfreq{eventvalue} = ft_freqbaseline(cfg, freq{eventvalue});
    
    % zscore transformation
    zscfreq{eventvalue} = relfreq{eventvalue};
    zscfreq{eventvalue}.powspctrm = zscore_transform(zscfreq{eventvalue}.powspctrm);
    
    % scaling by maximum
    maxfreq{eventvalue} = relfreq{eventvalue};
    maxfreq{eventvalue}.powspctrm = maxtransform(maxfreq{eventvalue}.powspctrm);
    
end



%% Compare zscore & max correction
% heatmaps
figure;
cfg = [];
cfg.baseline = 'no';

cfg.channel = freq{1}.label(62:64);

cfg.showlabels = 'yes';
cfg.layout = 'biosemi64.lay';

subplot(3, 1, 1);
ft_singleplotTFR(cfg, relfreq{1});

subplot(3, 1, 2);
ft_singleplotTFR(cfg, zscfreq{1});

subplot(3, 1, 3);
ft_singleplotTFR(cfg, maxfreq{1});

% line plot
figure;

x = freq{1}.time;
channels = 62:64;

subplot(3, 1, 1);
y28 = squeeze(mean(relfreq{1}.powspctrm(channels, index28, :), 1));
y36 = squeeze(mean(relfreq{1}.powspctrm(channels, index36, :), 1));
plot(x, y28, 'b', x, y36, 'g', x, y28-y36, 'r');

subplot(3, 1, 2);
y28 = squeeze(mean(zscfreq{1}.powspctrm(channels, index28, :), 1));
y36 = squeeze(mean(zscfreq{1}.powspctrm(channels, index36, :), 1));
plot(x, y28, 'b', x, y36, 'g', x, y28-y36, 'r');

subplot(3, 1, 3);
y28 = squeeze(mean(maxfreq{1}.powspctrm(channels, index28, :), 1));
y36 = squeeze(mean(maxfreq{1}.powspctrm(channels, index36, :), 1));
plot(x, y28, 'b', x, y36, 'g', x, y28-y36, 'r');


hhh;



%% Plots
cfg = [];
cfg.baseline     = [-2 -0.5];
cfg.baselinetype = 'relchange';
% cfg.zlim         = [0 15];
cfg.channel = freq.label(1:64);
cfg.showlabels   = 'yes';
cfg.layout       = 'biosemi64.lay';
% figure;
% ft_multiplotTFR(cfg, freq);

cfg.channel = freq.label(62:64);
figure;
ft_singleplotTFR(cfg, freq);


relFreq = ft_freqbaseline(cfg, freq);

pow36 = squeeze(relFreq.powspctrm(29, 29, :));
pow28 = squeeze(relFreq.powspctrm(29, 11, :));

% figure;
% plot(freq.time, pow36)
% hold on
% plot(freq.time, pow28, 'r')

% subplot(diffFig{iTrials});
% hold on
% plot(freq.time, pow36 - pow28);


%% ICA
cfg = [];
cfg.resamplerefs = 150;
cfg.detrend = 'no';

downsample_data = ft_resampledata(cfg, loaded_data);

cfg = [];
cfg.method = 'runica';
comp = ft_componentanalysis(cfg, downsample_data);

cfg = [];
cfg.component = 1:20;
cfg.layout = 'biosemi64.lay';
cfg.comment = 'no';
figure;
ft_topoplotIC(cfg, comp);


%% Include button data

% continButton = zeros(size(freq.time));
% 
% for i = 1:size(buttonCfg.trl, 1)
%     binVector = dec2bin(buttonCfg.trl(i, 4)-1, 3);
%     continButton(freq.time > buttonCfg.time(i)) = binVector(3) - binVector(1);
% end
% continButton(freq.time > 12) = 0;

% if angleOrder(1, trialOfInterest - 200) ~= flickerOrder(1, trialOfInterest - 200)
%     summedButton = summedButton * (-1);
% end
% 




%% Analyse all trials together

% cfg.trialdef.eventvalue = [12, 22];
% 
% cfg = ft_definetrial(cfg);
% loaded_data = ft_preprocessing(cfg);
% 
% cfg.toi = -2:0.05:14;
% cfg.foi = 26:0.4:40;
% 
% cfg.output       = 'pow';
% cfg.method       = 'mtmconvol';
% cfg.taper        = 'hanning';
% 
% cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;
% freq = ft_freqanalysis(cfg, loaded_data);

%% Multiplot
% cfg.baseline     = [-1 0]; 
% cfg.baselinetype = 'relative'; 
% cfg.zlim         = [0 15];	        
% cfg.channel = freq.label(1:64);
% cfg.showlabels   = 'yes';
% cfg.layout       = 'biosemi64.lay';
% figure;
% ft_multiplotTFR(cfg, freq);