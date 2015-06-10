%% Setup!
%#ok<*AGROW>
clearvars;  %#ok<*SAGROW> 
clc;
load('buttons_freqs.mat');
if exist('E:\Documents\Recorded Data\EEG Feb 2015', 'dir') % location on desktop
    file_directory = 'E:\Documents\Recorded Data\EEG Feb 2015';
elseif exist('D:\Recorded Data', 'dir') % location on laptop
    file_directory = 'D:\Recorded Data';
else
    error('please provide directory where file is stored');
end


%% File Names
filenames{1} = dir([file_directory, '\*CTR*BinSSVEP.bdf']);
filenames{2} = dir([file_directory, '\*ASC*BinSSVEP.bdf']);
n{1} = size(filenames{1}, 1);
n{2} = size(filenames{2}, 1);

%% Global Variables
trial_dur = 12;
plotting = 0;
remove_subject{1} = false(n{1}, 1);
remove_subject{2} = false(n{2}, 1);

for group = 1:2
for subject = 1:n{group}
    %% Definition of Trials
    fileID = fullfile(file_directory, filenames{group}(subject).name);
    official_ID = str2double(fileID(end-16:end-13));
    cfg = [];
    cfg.dataset = fileID;
    cfg.channel = 1:64;
    cfg.trialdef.eventtype = 'STATUS';
    cfg.trialfun = 'ft_trialfun_general';
    cfg.trialdef.prestim = 3;
    cfg.trialdef.poststim = 17; %actual length of trial 12 s
    
        cfg.trialdef.eventvalue = 201:216;

        try
            cfg = ft_definetrial(cfg);
        catch define_trial_error
            cfg.trialdef.eventvalue
            cfg.dataset
            rethrow(define_trial_error);
        end
        
        cfg.trl = remove_overlaps(cfg); % this script removes overlapping trials
        samplefreq = abs(cfg.trl(1, 3)) / cfg.trialdef.prestim;
        offset = abs(cfg.trl(1, 3));
        
    for eventvalue = 201:216
        
        cfg_trial = cfg;
        cfg_trial.trl = cfg.trl( cfg.trl(:,4)==eventvalue, : );
        trial_starts(eventvalue-200) = cfg_trial.trl(1, 1) - cfg_trial.trl(1, 3);
%         trial_data{eventvalue-200} = ft_preprocessing(cfg_trial);
        
    end
    
    %% Behavioural Analysis (from triggers)
    button_data = cell(1, 16);
    load('buttons_freqs.mat');
    
    % Extract the relevant sequence
    cfg_buttons = [];
    cfg_buttons.dataset = fileID;
    cfg_buttons.trialdef.eventvalue = 1:8;
    cfg_buttons.trialdef.eventtype = 'STATUS';
    cfg_buttons.trialfun = 'ft_trialfun_general';
    try
        cfg_buttons = ft_definetrial(cfg_buttons);
    catch button_define_error
        warning(button_define_error.message);
        remove_subject{group}(subject) = 1;
        continue
    end
    
    % Now analyse trial by trial
    for trialNo = 1:16
        trl_start = trial_starts(trialNo);
        trl_end = trl_start + samplefreq*trial_dur;
        
        button_data{trialNo} = cfg_buttons.trl( cfg_buttons.trl(:, 1) > trl_start & cfg_buttons.trl(:, 1) < trl_end, : );
        if size(button_data{trialNo}, 1) < 1
            break % if no buttons were pressed, this trial needs to be skipped
        end
        
        % convert the second column to real time values (in secs, relative
        % to trial onset)
        button_data{trialNo}(:, 2) = (button_data{trialNo}(:, 2) - trl_start)/samplefreq;
        
        percepts(trialNo).start = button_data{trialNo}(:, 2);
        percepts(trialNo).buttons = dec2bin(button_data{trialNo}(:, 4)-1, 3)-'0';
        
        % clean the buttonpresses (remove doubles)
        [percepts(trialNo).buttons, percepts(trialNo).start] = remove_doublebtns(percepts(trialNo).buttons, percepts(trialNo).start);
        
        % parse the percepts
        [percepts(trialNo).buttons, percepts(trialNo).duration] = parse_percepts(percepts(trialNo).start, percepts(trialNo).buttons, trial_dur);
        
        % clean the percepts
        [percepts(trialNo).buttons, percepts(trialNo).duration, percepts(trialNo).start] = ...
            remove_tooshort(percepts(trialNo).buttons, percepts(trialNo).duration, percepts(trialNo).start, 0.5);
        [percepts(trialNo).buttons, percepts(trialNo).duration, percepts(trialNo).start] = ...
            remove_last(percepts(trialNo).buttons, percepts(trialNo).duration, percepts(trialNo).start, trial_dur);
        
        % determine index for percepts
        [percepts(trialNo).ccw_index, percepts(trialNo).cw_index, percepts(trialNo).mix_index] = find_percept_index(percepts(trialNo).buttons);
    end
    
    % determine median and mean percepts
    all_mix = [];
    all_cw = [];
    all_ccw = [];
    all_dom = [];
    mix_prop = [];
    for trialNo = 1:16
        all_mix = [all_mix; percepts(trialNo).duration(percepts(trialNo).mix_index)];
        all_cw = [all_cw; percepts(trialNo).duration(percepts(trialNo).cw_index)];
        all_ccw = [all_ccw; percepts(trialNo).duration(percepts(trialNo).ccw_index)];
        all_dom = [all_dom; percepts(trialNo).duration(percepts(trialNo).cw_index | percepts(trialNo).ccw_index)];
        mix_prop = [mix_prop; sum(percepts(trialNo).duration(percepts(trialNo).mix_index))/sum(percepts(trialNo).duration(percepts(trialNo).cw_index | percepts(trialNo).ccw_index))];
    end
    
    dom_num{group}(subject, 1) = numel(all_dom);
    if dom_num{group}(subject, 1) < 20
        remove_subject{group}(subject) = 1;
    end
    
    if ~isempty(all_cw) && ~isempty(all_cw)
        bias{group}(subject, 1) = ttest2(all_cw, all_ccw);
    else
        bias{group}(subject, 1) = NaN;
        remove_subject{group}(subject) = 1;
    end
    
    
    mean_mixprop{group}(subject, 1) = mean(mix_prop);
    mean_mix{group}(subject, 1) = mean(all_mix);
    median_mix{group}(subject, 1) = median(all_mix);
    mean_dom{group}(subject, 1) = mean(all_dom);
    median_dom{group}(subject, 1) = median(all_dom);
    
    %% SSVEP Analysis
    
    % preprocessing
    cfg_preproc = cfg;
    cfg_preproc.continuous = 'yes';
    cfg_preproc.demean    = 'yes';
    cfg_preproc.detrend = 'yes';
    cfg_preproc.reref = 'yes';
    cfg_preproc.refchannel = 1:64;
    for eventvalue = 201:216
        cfg_preproc.trl = cfg.trl( cfg.trl(:,4)==eventvalue, : );
        trial_data{eventvalue-200} = ft_preprocessing(cfg_preproc);
    end
    
    % TFR
    
%     freq_res = 0.4;
    
    cfg_tfr = [];
    cfg_tfr.method = 'wavelet';
    cfg_tfr.width = 16;
    cfg_tfr.foi = [28.8, 36];
    cfg_tfr.toi = -3:0.05:12;
    
    for eventvalue = 1:16
    
        freqs{eventvalue} = ft_freqanalysis(cfg_tfr, trial_data{eventvalue});
        cfg_tfr.toi = freqs{eventvalue}.time;
    end
    
    %% Plotting
    if plotting
    figure;
    cfg_plot = [];
    cfg_plot.baseline = [-3, 0];
    cfg_plot.baselinetype = 'relative';
    cfg_plot.parameter = 'powspctrm';
    
    cfg_plot.channel = freqs{1}.label([27, 29, 64]);
    cfg_plot.channel = {'Oz'; 'POz'; 'O2'; 'O1'; 'PO3'; 'PO4'};

    cfg_plot.showlabels = 'yes';
    cfg_plot.layout = 'biosemi64.lay';
    for plttrl = 1:16
        subplot(8, 2, plttrl);
        ft_singleplotTFR(cfg_plot, freqs{plttrl});
    end
    end
    
    %% Analyse based on Buttons
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % timecourse around 1s of dominance
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    amps_dom{group, subject} = [];
    amps_sup{group, subject} = [];
    for trialNo = 1:16
        if ~isequal(buttons_freqs(official_ID).flickerOrder(:, trialNo), buttons_freqs(official_ID).angleOrder(:, trialNo))
            ssvep_index_36 = percepts(trialNo).start > 1 & percepts(trialNo).duration > 1 & percepts(trialNo).cw_index;
            ssvep_index_28 = percepts(trialNo).start > 1 & percepts(trialNo).duration > 1 & percepts(trialNo).ccw_index;
        else
            ssvep_index_28 = percepts(trialNo).start > 1 & percepts(trialNo).duration > 1 & percepts(trialNo).cw_index;
            ssvep_index_36 = percepts(trialNo).start > 1 & percepts(trialNo).duration > 1 & percepts(trialNo).ccw_index;
        end
        
        mean_36 = mean( mean(freqs{trialNo}.powspctrm(29:30, 2, (freqs{trialNo}.time > 1 & freqs{trialNo}.time < 12)), 1 ), 3);
        mean_28 = mean( mean(freqs{trialNo}.powspctrm(29:30, 1, (freqs{trialNo}.time > 1 & freqs{trialNo}.time < 12)), 1 ), 3);
        
        for timeSlice = find(ssvep_index_36)'
            
            t_0 = percepts(trialNo).start(timeSlice);
            
            amp = mean(freqs{trialNo}.powspctrm(29:30, 2, (freqs{trialNo}.time > t_0-1 & freqs{trialNo}.time < t_0+1.5)), 1 );
            amp = permute(amp, [1, 3, 2]) / mean_36;
            amp = padarray(amp, [0, 50-size(amp, 2)], NaN, 'post');
            amps_dom{group, subject} = [amps_dom{group, subject}; amp];
            
            amp = mean(freqs{trialNo}.powspctrm(29:30, 1, (freqs{trialNo}.time > t_0-1 & freqs{trialNo}.time < t_0+1.5)), 1 );
            amp = permute(amp, [1, 3, 2]) / mean_28;
            amp = padarray(amp, [0, 50-size(amp, 2)], NaN, 'post');
            amps_sup{group, subject} = [amps_sup{group, subject}; amp];
        end
        
        for timeSlice = find(ssvep_index_28)'
            
            t_0 = percepts(trialNo).start(timeSlice);
            
            amp = mean(freqs{trialNo}.powspctrm(29:30, 1, (freqs{trialNo}.time > t_0-1 & freqs{trialNo}.time < t_0+1.5)), 1 );
            amp = permute(amp, [1, 3, 2]) / mean_28;
            amp = padarray(amp, [0, 50-size(amp, 2)], NaN, 'post');
            amps_dom{group, subject} = [amps_dom{group, subject}; amp];
            
            amp = mean(freqs{trialNo}.powspctrm(29:30, 2, (freqs{trialNo}.time > t_0-1 & freqs{trialNo}.time < t_0+1.5)), 1 );
            amp = permute(amp, [1, 3, 2]) / mean_36;
            amp = padarray(amp, [0, 50-size(amp, 2)], NaN, 'post');
            amps_sup{group, subject} = [amps_sup{group, subject}; amp];
        end
        
    end
    
    if ~isempty( amps_dom{group, subject} ) && ~isempty( amps_sup{group, subject} )
        avg_amp_dom{group}(subject, :) = nanmedian( amps_dom{group, subject}, 1 );
        avg_amp_sup{group}(subject, :) = nanmedian( amps_sup{group, subject}, 1 );
    else
        avg_amp_dom{group}(subject, :) = NaN;
        avg_amp_sup{group}(subject, :) = NaN;
    end
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % average across all timebins
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i = 1:10
        postWhole_dom_binned{i} = [];
        postWhole_sup_binned{i} = [];
    end
    
    pre150_dom = [];
    postWhole_dom = [];
    pre150_sup = [];
    postWhole_sup = [];
    postWhole_mix36 = [];
    postWhole_mix28 = [];
    
    for trialNo = 1:16
        if ~isequal(buttons_freqs(official_ID).flickerOrder(:, trialNo), buttons_freqs(official_ID).angleOrder(:, trialNo))
            ssvep_index_36 = percepts(trialNo).start > 1 & percepts(trialNo).cw_index;
            ssvep_index_28 = percepts(trialNo).start > 1 & percepts(trialNo).ccw_index;
        else
            ssvep_index_28 = percepts(trialNo).start > 1 & percepts(trialNo).cw_index;
            ssvep_index_36 = percepts(trialNo).start > 1 & percepts(trialNo).ccw_index;
        end
        ssvep_index_mix = percepts(trialNo).start > 1 & percepts(trialNo).mix_index;
        
        mean_36 = mean( mean(freqs{trialNo}.powspctrm(29:30, 2, (freqs{trialNo}.time > 1 & freqs{trialNo}.time < 12)), 1 ), 3);
        mean_28 = mean( mean(freqs{trialNo}.powspctrm(29:30, 1, (freqs{trialNo}.time > 1 & freqs{trialNo}.time < 12)), 1 ), 3);
        
        
        for timeSlice = find(ssvep_index_36)'
            
            t_0 = percepts(trialNo).start(timeSlice);
            t_end = percepts(trialNo).start(timeSlice) + percepts(trialNo).duration(timeSlice);
            for i = 1:10
                t_end_bin(i) = percepts(trialNo).start(timeSlice) + i * 0.1 * percepts(trialNo).duration(timeSlice);
            end
            
            % %%%%%%%%%%%%%%%%%%
            % Dominant Frequency
            % %%%%%%%%%%%%%%%%%%
            % cut out relevant sample
            samples_postWhole = ...
                mean(freqs{trialNo}.powspctrm(29:30, 2, (freqs{trialNo}.time > t_0 & freqs{trialNo}.time < t_end)), 1);
            samples_pre150 = ...
                mean(freqs{trialNo}.powspctrm(29:30, 2, (freqs{trialNo}.time > t_0-150 & freqs{trialNo}.time < t_0)), 1);
            % scale to trial mean
            samples_postWhole = samples_postWhole / mean_36;
            samples_pre150 = samples_pre150 / mean_36;
            % reshape
            samples_postWhole = permute(samples_postWhole, [1, 3, 2]);
            samples_pre150 = permute(samples_pre150, [1, 3, 2]);
            
            postWhole_sup = [postWhole_dom, samples_postWhole];
            pre150_sup = [pre150_dom, samples_pre150];
            
            % do the same for the BINNED frequencies
            for i = 1:10
            samples_postWhole_binned{i} = ...
                mean(freqs{trialNo}.powspctrm(29:30, 2, (freqs{trialNo}.time > t_0 & freqs{trialNo}.time < t_end_bin(i))), 1);
            
            samples_postWhole_binned{i} = samples_postWhole_binned{i} / mean_36;
            
            samples_postWhole_binned{i} = permute(samples_postWhole_binned{i}, [1 3 2]);
            
            postWhole_sup_binned{i} = [postWhole_sup_binned{i}, samples_postWhole_binned{i}];
            end
            
            
            
            % %%%%%%%%%%%%%%%%%%%%
            % Suppressed Frequency
            % %%%%%%%%%%%%%%%%%%%%
            % cut out relevant sample
            samples_postWhole = ...
                mean(freqs{trialNo}.powspctrm(29:30, 1, (freqs{trialNo}.time > t_0 & freqs{trialNo}.time < t_end)), 1);
            samples_pre150 = ...
                mean(freqs{trialNo}.powspctrm(29:30, 1, (freqs{trialNo}.time > t_0-150 & freqs{trialNo}.time < t_0)), 1);
            % scale to trial mean
            samples_postWhole = samples_postWhole / mean_28;
            samples_pre150 = samples_pre150 / mean_28;
            % reshape
            samples_postWhole = permute(samples_postWhole, [1, 3, 2]);
            samples_pre150 = permute(samples_pre150, [1, 3, 2]);
            
            postWhole_sup = [postWhole_sup, samples_postWhole];
            pre150_sup = [pre150_sup, samples_pre150];
            
            % do the same for the BINNED frequencies
            for i = 1:10
            samples_postWhole_binned{i} = ...
                mean(freqs{trialNo}.powspctrm(29:30, 2, (freqs{trialNo}.time > t_0 & freqs{trialNo}.time < t_end_bin(i))), 1);
            
            samples_postWhole_binned{i} = samples_postWhole_binned{i} / mean_28;
            
            samples_postWhole_binned{i} = permute(samples_postWhole_binned{i}, [1 3 2]);
            
            postWhole_sup_binned{i} = [postWhole_sup_binned{i}, samples_postWhole_binned{i}];
            end
        end
        
        for timeSlice = find(ssvep_index_28)'
            j = j+1;
            
            t_0 = percepts(trialNo).start(timeSlice);
            t_end = percepts(trialNo).start(timeSlice) + percepts(trialNo).duration(timeSlice)-0.15;
            for i = 1:10
                t_end_bin(i) = percepts(trialNo).start(timeSlice) + i * 0.1 * percepts(trialNo).duration(timeSlice);
            end
            
            % %%%%%%%%%%%%%%%%%%%%
            % Suppressed Frequency
            % %%%%%%%%%%%%%%%%%%%%
            % cut out relevant sample
            samples_postWhole = ...
                mean(freqs{trialNo}.powspctrm(29:30, 2, (freqs{trialNo}.time > t_0 & freqs{trialNo}.time < t_end)), 1);
            samples_pre150 = ...
                mean(freqs{trialNo}.powspctrm(29:30, 2, (freqs{trialNo}.time > t_0-150 & freqs{trialNo}.time < t_0)), 1);
            % scale to trial mean
            samples_postWhole = samples_postWhole / mean_36;
            samples_pre150 = samples_pre150 / mean_36;
            % reshape
            samples_postWhole = permute(samples_postWhole, [1, 3, 2]);
            samples_pre150 = permute(samples_pre150, [1, 3, 2]);
            
            postWhole_sup = [postWhole_sup, samples_postWhole];
            pre150_sup = [pre150_sup, samples_pre150];
            
            % do the same for the BINNED frequencies
            for i = 1:10
            samples_postWhole_binned{i} = ...
                mean(freqs{trialNo}.powspctrm(29:30, 2, (freqs{trialNo}.time > t_0 & freqs{trialNo}.time < t_end_bin(i))), 1);
            
            samples_postWhole_binned{i} = samples_postWhole_binned{i} / mean_36;
            
            samples_postWhole_binned{i} = permute(samples_postWhole_binned{i}, [1 3 2]);
            
            postWhole_sup_binned{i} = [postWhole_sup_binned{i}, samples_postWhole_binned{i}];
            end
            
            % %%%%%%%%%%%%%%%%%%
            % Dominant Frequency
            % %%%%%%%%%%%%%%%%%%
            % cut out relevant sample
            samples_postWhole = ...
                mean(freqs{trialNo}.powspctrm(29:30, 1, (freqs{trialNo}.time > t_0 & freqs{trialNo}.time < t_end)), 1);
            samples_pre150 = ...
                mean(freqs{trialNo}.powspctrm(29:30, 1, (freqs{trialNo}.time > t_0-150 & freqs{trialNo}.time < t_0)), 1);
            % scale to trial mean
            samples_postWhole = samples_postWhole / mean_28;
            samples_pre150 = samples_pre150 / mean_28;
            % reshape
            samples_postWhole = permute(samples_postWhole, [1, 3, 2]);
            samples_pre150 = permute(samples_pre150, [1, 3, 2]);
            
            postWhole_dom = [postWhole_dom, samples_postWhole];
            pre150_dom = [pre150_dom, samples_pre150];
            
            % do the same for the BINNED frequencies
            for i = 1:10
            samples_postWhole_binned{i} = ...
                mean(freqs{trialNo}.powspctrm(29:30, 2, (freqs{trialNo}.time > t_0 & freqs{trialNo}.time < t_end_bin(i))), 1);
            
            samples_postWhole_binned{i} = samples_postWhole_binned{i} / mean_28;
            
            samples_postWhole_binned{i} = permute(samples_postWhole_binned{i}, [1 3 2]);
            
            postWhole_dom_binned{i} = [postWhole_dom_binned{i}, samples_postWhole_binned{i}];
            end
        end
        
        for timeSlice = find(ssvep_index_mix)'
            j = j+1;
            
            t_0 = percepts(trialNo).start(timeSlice);
            t_end = percepts(trialNo).start(timeSlice) + percepts(trialNo).duration(timeSlice)-0.15;
            
            % %%%%%%%%%%%%%%%
            % 36 Hz Frequency
            % %%%%%%%%%%%%%%%
            % cut out relevant sample
            samples_postWhole = ...
                mean(freqs{trialNo}.powspctrm(29:30, 2, (freqs{trialNo}.time > t_0 & freqs{trialNo}.time < t_end)), 1);
            % scale to trial mean
            samples_postWhole = samples_postWhole / mean_36;
            % reshape
            samples_postWhole = permute(samples_postWhole, [1, 3, 2]);
            
            postWhole_mix36 = [postWhole_mix36, samples_postWhole];
            
            % %%%%%%%%%%%%%%%
            % 28 Hz Frequency
            % %%%%%%%%%%%%%%%
            % cut out relevant sample
            samples_postWhole = ...
                mean(freqs{trialNo}.powspctrm(29:30, 1, (freqs{trialNo}.time > t_0 & freqs{trialNo}.time < t_end)), 1);
            % scale to trial mean
            samples_postWhole = samples_postWhole / mean_28;
            % reshape
            samples_postWhole = permute(samples_postWhole, [1, 3, 2]);
            
            postWhole_mix28 = [postWhole_mix28, samples_postWhole];
            
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Was the previous percept a dominant one?
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        end
        
    end
    
    % store for group
%     pre150_dom_av{group}(subject, :) = nanmean( pre150_dom );
%     pre150_sup_av{group}(subject, :) = nanmean( pre150_sup );
    
    postWhole_dom_av{group}(subject, :) = nanmean( postWhole_dom );
    postWhole_sup_av{group}(subject, :) = nanmean( postWhole_sup );
%     postWhole_mix28_av{group}(subject, :) = nanmean( postWhole_mix28 );
%     postWhole_mix36_av{group}(subject, :) = nanmean( postWhole_mix36 );
    
%     for i = 1:10
%         postWhole_dom_binned_av{group}(subject, i) = nanmean( postWhole_dom_binned{i} );
%         postWhole_sup_binned_av{group}(subject, i) = nanmean( postWhole_sup_binned{i} );
%     end
    
    %% Analyse Correlation between frequencies
    % concatenate each trial sequence, excluding the first second of a
    % trial
    all_samples_36 = [];
    all_samples_28 = [];
    for trialNo = 1:16
        % scale each frequency to the trial mean
        mean_36 = mean( mean(freqs{trialNo}.powspctrm(29:30, 2, (freqs{trialNo}.time > 1 & freqs{trialNo}.time < 12)), 1 ), 3);
        mean_28 = mean( mean(freqs{trialNo}.powspctrm(29:30, 1, (freqs{trialNo}.time > 1 & freqs{trialNo}.time < 12)), 1 ), 3);
        
%         normalised_36 = permute( mean(freqs{trialNo}.powspctrm(29:30, 2, (freqs{trialNo}.time > 1 & freqs{trialNo}.time < 12)), 1 ) / mean_36, [3, 1, 2]);
%         normalised_28 = permute( mean(freqs{trialNo}.powspctrm(29:30, 1, (freqs{trialNo}.time > 1 & freqs{trialNo}.time < 12)), 1 ) / mean_28, [3, 1, 2]);
        not_normalised_36 = permute( mean(freqs{trialNo}.powspctrm(29:30, 2, (freqs{trialNo}.time > 1 & freqs{trialNo}.time < 12)), 1 ), [3, 1, 2]);
        not_normalised_28 = permute( mean(freqs{trialNo}.powspctrm(29:30, 1, (freqs{trialNo}.time > 1 & freqs{trialNo}.time < 12)), 1 ), [3, 1, 2]);
        
        all_samples_36 = [all_samples_36; not_normalised_36];
        all_samples_28 = [all_samples_28; not_normalised_28];
    end
    
    pearson_r{group}(subject, 1) = corr(all_samples_36, all_samples_28);
    
end
end



    %% Plot based on buttons
    % whole percept
    figure;
    hold on;
    linespec_dom{1} = '--ko'; linespec_dom{2} = ':kv';
    linespec_sup{1} = '--ro'; linespec_sup{2} = ':rv';
    for group = 1:2
    errorbar(0.8, nanmean(postWhole_dom_av{group}(~remove_subject{group})), nansem(postWhole_dom_av{group}(~remove_subject{group})), linespec_dom{group});
    errorbar(1.2, nanmean(postWhole_sup_av{group}(~remove_subject{group})), nansem(postWhole_sup_av{group}(~remove_subject{group})), linespec_sup{group});
%     errorbar(1.05, nanmean(postWhole_mix28_av{group}(~remove_subject{group})), nansem(postWhole_mix28_av{group}(~remove_subject{group})), 'go');
%     errorbar(0.95, nanmean(postWhole_mix36_av{group}(~remove_subject{group})), nansem(postWhole_mix36_av{group}(~remove_subject{group})), 'go');
    % 150ms before percept
%     errorbar(-1.1, nanmean(pre150_dom_av{group}(~remove_subject{group})), nansem(pre150_dom_av{group}(~remove_subject{group})), 'ko');
%     errorbar(-0.9, nanmean(pre150_sup_av{group}(~remove_subject{group})), nansem(pre150_sup_av{group}(~remove_subject{group})), 'ro');
    end
    
    hhh;
    
    % suppression index?
    figure;
    hold on;
    for group = 1:2
    suppression_index{group} = postWhole_dom_av{group} ./ postWhole_sup_av{group};
%     bar(group, nanmedian(suppression_index(~remove_subject{group})), 0.4);
    errorbar(group, nanmean(suppression_index{group}(~remove_subject{group})), nansem(suppression_index{group}(~remove_subject{group})), 'o');
    xlim([0, 3]);
    end
    
    % time course around switch
    figure;
    hold on;
    plot( nanmedian(avg_amp_sup{group}( ~remove_subject{group}, : ), 1), 'r' );
    plot( nanmedian(avg_amp_sup{group}( ~remove_subject{group}, : ), 1) - nanstd(avg_amp_sup{1}, 1)/sqrt(size(avg_amp_sup{1}, 1)), 'r:');
    plot( nanmedian(avg_amp_sup{group}( ~remove_subject{group}, : ), 1) + nanstd(avg_amp_sup{1}, 1)/sqrt(size(avg_amp_sup{1}, 1)), 'r:');
    plot( nanmedian(avg_amp_dom{group}( ~remove_subject{group}, : ), 1), 'g' );
    plot( nanmedian(avg_amp_dom{group}( ~remove_subject{group}, : ), 1) - nanstd(avg_amp_dom{1}, 1)/sqrt(size(avg_amp_dom{1}, 1)), 'g:');
    plot( nanmedian(avg_amp_dom{group}( ~remove_subject{group}, : ), 1) + nanstd(avg_amp_dom{1}, 1)/sqrt(size(avg_amp_dom{1}, 1)), 'g:');
    
    % Old style first second
    figure;
    hold on;
    for group = 1:2
        
        first_second_dom = nanmean(nanmean(avg_amp_dom{group}( ~remove_subject{group}, 20:40 ), 2), 1);
        first_second_sup = nanmean(nanmean(avg_amp_sup{group}( ~remove_subject{group}, 20:40 ), 2), 1);
        first_second_dom_sem = nansem(nanmean(avg_amp_dom{group}( ~remove_subject{group}, 20:40 ), 2), 1);
        first_second_sup_sem = nansem(nanmean(avg_amp_sup{group}( ~remove_subject{group}, 20:40 ), 2), 1);
        
        errorbar([1, 2], [first_second_dom, first_second_sup], [first_second_dom_sem, first_second_sup_sem], linespec_dom{group});
        
    end
    
    % time course across a dominant percept
    figure;
    hold on;
    linespec_dom{1} = '--ko'; linespec_dom{2} = ':kv';
    linespec_sup{1} = '--ro'; linespec_sup{2} = ':rv';
    for group = 1:2
        errorbar(1:10, nanmean(postWhole_dom_binned_av{group}, 1), nansem(postWhole_dom_binned_av{group}, 1), linespec_dom{group} );
        errorbar(1:10, nanmean(postWhole_sup_binned_av{group}, 1), nansem(postWhole_sup_binned_av{group}, 1), linespec_sup{group} );
    end
    figure;
    hold on;
    for group = 1:2
        errorbar( nanmean(postWhole_dom_binned_av{group}-postWhole_sup_binned_av{group}, 1), nansem(postWhole_dom_binned_av{group}-postWhole_sup_binned_av{group}, 1), linespec_dom{group} );
    end
    

%% Button Press from trigger data

button_data = cell(1, 16);
button_time = cell(1, 16);
button_sum = cell(1, 16);
load('buttons_freqs.mat');
for trialNo = 1:16;
    cfg = [];
    cfg.dataset = fullfile(file_directory, fileID);
    cfg.trialdef.eventvalue = 1:8;
    cfg.trialdef.eventtype = 'STATUS';
    cfg.trialfun = 'ft_trialfun_general';
    
    try
        cfg = ft_definetrial(cfg);
    catch err
        warning(err.message);
        continue
    end
    % only take the buttonpresses within current trial
    
    trl_start = trial_data{trialNo}.sampleinfo(1) + trial_data{trialNo}.fsample* 2;
    trl_end = trial_data{trialNo}.sampleinfo(1) + trial_data{trialNo}.fsample* (2+12);
    
    cfg.trl(cfg.trl(:, 1) < trl_start | cfg.trl(:, 1) > trl_end, :) = [];
    
    button_data{trialNo} = cfg.trl;
    if size(button_data{trialNo}, 1) < 1
        break
    end
    button_data{trialNo}(:, 2) = (button_data{trialNo}(:, 2) - trl_start)/trial_data{trialNo}.fsample;
    
    % analyse the behavioural response
    [percepts(trialNo).buttons, percepts(trialNo).duration] = parse_percepts(button_data{trialNo}(:, 2), dec2bin(button_data{trialNo}(:, 4)-1, 3), 12);
    
    % make time (X) for plotting button presses
    button_time{trialNo} = zeros(1, 2*size(button_data{trialNo}, 1));
    button_time{trialNo}(1, 1:2:end-1) = button_data{trialNo}(1:end, 1);
    button_time{trialNo}(1, 2:2:end-2) = button_data{trialNo}(2:end, 1);
    button_time{trialNo}(1, end) = button_data{trialNo}(1, 1) + 12*1024;
    % start at time 0, go to 12
    button_time{trialNo} = (button_time{trialNo} - button_time{trialNo}(1, 1)) /1024;
    
    
    % duplicate buttonpresses for plotting
    
    bin_buttons = dec2bin(button_data{trialNo}(:, 4)-1, 3);
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




%% Scaling
% In this part we scale either by maximum, or do a Z-Score transform!
% 
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


%% tfr plot
% This makes heatmaps of the frequencies
figure;
cfg = [];
cfg.baseline = [-5, 0-1/freq_res];
cfg.baselinetype = 'relative';
cfg.parameter = 'powspctrm';

cfg.channel = freq{1}.label([27, 29, 64]);
cfg.channel = {'Oz'; 'POz'; 'O2'; 'O1'; 'PO3'; 'PO4'};

cfg.showlabels = 'yes';
cfg.layout = 'biosemi64.lay';
for plttrl = 1:16
    subplot(8, 2, plttrl);
    ft_singleplotTFR(cfg, freq{plttrl});
    
end

hhh;


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
y28 = squeeze(nanmean(maxfreq{plttrl}.powspctrm(channels, index28, :), 1));
y36 = squeeze(nanmean(maxfreq{plttrl}.powspctrm(channels, index36, :), 1));
% plot(x, y28, 'b--', x, y36, 'g--', x, y28-y36, 'r', button_time{plttrl}, button_sum{plttrl}, 'm');
% title('Divide by maximum');



end
hhh;


%% Analyse correlation between two frequencies
for trl = 1:16

% Identify stimulated timecourse
x = freq{1}.time( freq{1}.time > 0 & freq{1}.time < 12 );

pow28 = squeeze(nanmean(maxfreq{trl}.powspctrm(channels, index28, ( freq{1}.time > 0 & freq{1}.time < 12 )), 1));
pow36 = squeeze(nanmean(maxfreq{trl}.powspctrm(channels, index36, ( freq{1}.time > 0 & freq{1}.time < 12 )), 1));

[r, p] = corrcoef(pow28, pow36);
disp(['Trial ' num2str(trl)])
disp(['r=' num2str(r(2)) ' p=' num2str(p(2))])

end