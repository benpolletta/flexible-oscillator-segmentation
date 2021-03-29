function runFig6simulations

[model_vary, model_names, include_deepFS] = getModelParams;

no_bands = 8;
cochlearBands = getCochlearBands(no_bands)';
% cochlearBands = cochlearBands(:, 1:2:end); % cochlearBands = cochlearBands(:, 4:8:end); %8:16:end); % 

simulation_vary = {%'deepRS', 'mechanism_list', '+iSpeechInput';...
    'deepRS', 'SentenceIndex', .025:.025:1;... .001:.001:1;... .2:.2:1;... [.5 1];... 
    'deepRS', 'gSpeech', 2.2:0.2:3;... [0 .5 1 2];... 0:2;... 0:.1:2;... 0.1:0.1:0.5;... 
    'deepRS', '(SpeechLowFreq, SpeechHighFreq)', cochlearBands;...
    'deepRS', 'SpeechNorm', 0;...
    };

names = cell(length(model_names), 1);

Today = datestr(datenum(date),'yy-mm-dd');

Now = clock;

for m = 1:length(model_vary)
    
    model_name = model_names{m}
    
    vary = [model_vary{m}; simulation_vary]
    
    sim_struct = init_sim_struct('vary', vary, 'include_deepRS', 1, 'NdeepRS', 128/no_bands,...
        'include_deepFS', include_deepFS(m), 'NdeepFS', 128/no_bands, 'tspan', [0 6000],...
        'cluster_flag', 1, 'parallel_flag', 0, 'save_data_flag', 0, 'sim_mode', 3,...
        'analysis_functions', {@boundary_analysis, @phase_metrics}, 'analysis_options',...
        {{'threshold', .333, 'sum_window', 25, 'synchrony_window', 50, 'vp_norm', 3},...
        {'i_field', 'deepRS_iSpeechInput_input', 'i_transform', 'wavelet_reduce',...
        'no_periods', -2, 'p_field', 'deepRS_iSpeechInput_syllables'}},...
        'note', ['TIMIT segmentation simulations for model ', model_name],...
        'parent', model_names{m});
    
    [~, names{m}] = simulateModel(sim_struct);

end

sim_name = sprintf('Fig6sims_%s_%d+%d', Today, Now([4 5]));

save(['Sims/', Today, '/', sim_name, '.mat'], 'sim_name', 'names')