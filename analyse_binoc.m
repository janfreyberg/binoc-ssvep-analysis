
if exist('E:\Documents\Recorded Data\EEG Feb 2015', 'dir') % location on desktop
    file_directory = 'E:\Documents\Recorded Data\EEG Feb 2015';
elseif exist('D:\Recorded Data', 'dir') % location on laptop
    file_directory = 'D:\Recorded Data';
else
    error('please provide directory where file is stored');
end

%% Pre-processing parameters
cfg = [];
cfg.continuous = 'yes';
cfg.demean    = 'yes';
cfg.detrend = 'yes';
% cfg.lpfreq = 49;
% cfg.hpfreq = 20;
cfg.reref = '1:64';
cfg.refchannel = 'all';


%% Trial definition & preprocessing
participantID = 8;

fileID = sprintf('EG-CTR-0008-BinSSVEP.bdf');
cfg.dataset = fullfile(file_directory, fileID);
cfg.channel = 1:64;

cfg.trialdef.eventtype = 'STATUS';
cfg.trialfun = 'ft_trialfun_general';
cfg.trialdef.prestim = 2;
cfg.trialdef.poststim = 15.5; %actual length of trial 12 s

trial_data = cell(1, 16);
for eventvalue = 1:16;
    
    cfg.trialdef.eventvalue = 200 + eventvalue;
    
    cfg.trl = [];
    cfg = rmfield(cfg, 'trl');
    cfg = ft_definetrial(cfg);
    
    cfg.trl = remove_overlaps(cfg.trl); % this scripts removes overlapping trials
    
    trial_data{eventvalue} = ft_preprocessing(cfg);
    
end


%% Button Press from trigger data

button_data = cell(1, 16);
button_time = cell(1, 16);
button_sum = cell(1, 16);
load('buttons_freqs.mat');
for trialNo = 1:16;
    cfg = [];
    fileID = 'EG-CTR-0008-BinSSVEP.bdf';
    cfg.dataset = fullfile(file_directory, fileID);
    cfg.trialdef.eventvalue = 1:8;
    cfg.trialdef.eventtype = 'STATUS';
    cfg.trialfun = 'ft_trialfun_general';
    
    cfg = ft_definetrial(cfg);
    
    % only take the buttonpresses within current trial
    
    trl_start = trial_data{trialNo}.sampleinfo(1) + trial_data{trialNo}.fsample* 2;
    trl_end = trial_data{trialNo}.sampleinfo(1) + trial_data{trialNo}.fsample* (2+12);
    
    cfg.trl(cfg.trl(:, 1) < trl_start | cfg.trl(:, 1) > trl_end, :) = [];
    
    button_data{trialNo} = cfg.trl;
    
    % make time (X) for plotting button presses
    button_time{trialNo} = zeros(1, 2*size(button_data{trialNo}, 1));
    button_time{trialNo}(1, 1:2:end-1) = button_data{trialNo}(1:end, 1);
    button_time{trialNo}(1, 2:2:end-2) = button_data{trialNo}(2:end, 1);
    button_time{trialNo}(1, end) = button_data{trialNo}(1, 1) + 12*1024;
    % start at time 0, go to 12
    button_time{trialNo} = (button_time{trialNo} - button_time{trialNo}(1, 1)) /1024;
    
    
    % duplicate buttonpresses for plotting
    
    bin_buttons = dec2bin(button_data{trialNo}(:, 4)-1);
    button_sum{trialNo} = zeros(1, 2*size(button_data{trialNo}, 1));
    
    if isequal(buttons_freqs(participantID).flickerOrder(:, trialNo), buttons_freqs(participantID).angleOrder(:, trialNo))
        button_sum{trialNo}(2:2:end) = bin_buttons(:, 1) - bin_buttons(:, 3);
        button_sum{trialNo}(1:2:end-1) = bin_buttons(:, 1) - bin_buttons(:, 3);
    elseif ~isequal(buttons_freqs(participantID).flickerOrder(:, trialNo), buttons_freqs(participantID).angleOrder(:, trialNo))
        button_sum{trialNo}(2:2:end) = bin_buttons(:, 3) - bin_buttons(:, 1);
        button_sum{trialNo}(1:2:end-1) = bin_buttons(:, 3) - bin_buttons(:, 1);
    end
    
    
    % now make sure the right numbers are subtracted from each other
    % angle order and freq order are stored for each participant here:
    
    
    
end


%% TFR Parameters
cfg = [];
cfg.toi = -2:0.05:15.5;
% cfg.foi = 24.8:0.4:40;
cfg.foi = [28.8, 36];

cfg.output = 'pow';
cfg.method = 'wavelet'; % other option would be multi-taper-convolution ('mtmconvol')
cfg.width = cfg.foi/2; % width of the wavelet window in terms of wave cycles

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
    zscfreq{eventvalue}.powspctrm = zscore_transform(zscfreq{eventvalue}.powspctrm( :, :, zscfreq{eventvalue}.time > 0 & zscfreq{eventvalue}.time < 12 ));
    
    % scaling by maximum
    maxfreq{eventvalue} = relfreq{eventvalue};
    maxfreq{eventvalue}.powspctrm = maxtransform(maxfreq{eventvalue}.powspctrm);
    
end


%% Compare zscore & max correction
figure;
for plttrl = 1:16

% heatmaps
% figure;
% cfg = [];
% cfg.baseline = 'no';
% 
% cfg.channel = freq{plttrl}.label(62:64);
% 
% cfg.showlabels = 'yes';
% cfg.layout = 'biosemi64.lay';
% 
% subplot(3, 1, 1);
% ft_singleplotTFR(cfg, relfreq{plttrl});
% title('No transform');
% 
% subplot(3, 1, 2);
% ft_singleplotTFR(cfg, zscfreq{plttrl});
% title('Z-Score Transformation');
% 
% subplot(3, 1, 3);
% ft_singleplotTFR(cfg, maxfreq{plttrl});
% title('Divide by maximum');


% line plot
subplot(8, 2, plttrl);

x = freq{1}.time( freq{1}.time > 0 & freq{1}.time < 12 );
channels = 62:64;

% subplot(3, 1, 1);
% y28 = squeeze(mean(relfreq{plttrl}.powspctrm(channels, index28, :), 1));
% y36 = squeeze(mean(relfreq{plttrl}.powspctrm(channels, index36, :), 1));
% plot(x, y28, 'b--', x, y36, 'g--', x, y28-y36, 'r', button_time{plttrl}, button_sum{plttrl}, 'm');
% title('No transform');

% subplot(3, 1, 2);
y28 = squeeze(mean(zscfreq{plttrl}.powspctrm(channels, index28, :), 1));
y36 = squeeze(mean(zscfreq{plttrl}.powspctrm(channels, index36, :), 1));
plot(x, y28, 'b--', x, y36, 'g--', x, y28-y36, 'r', button_time{plttrl}, button_sum{plttrl}, 'm');
% title('Z-Score Transformation');

% subplot(3, 1, 3);
% y28 = squeeze(mean(maxfreq{plttrl}.powspctrm(channels, index28, :), 1));
% y36 = squeeze(mean(maxfreq{plttrl}.powspctrm(channels, index36, :), 1));
% plot(x, y28, 'b--', x, y36, 'g--', x, y28-y36, 'r', button_time{plttrl}, button_sum{plttrl}, 'm');
% title('Divide by maximum');

end
hhh;


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


