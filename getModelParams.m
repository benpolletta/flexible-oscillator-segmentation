function [model_vary, model_names, include_deepFS] = getModelParams

model_names = {'M', 'MI', 'I', 'IS', 'MIS', 'MS'};

include_deepFS = [0 1 1 1 1 0];

model_vary{1} = {'deepRS', 'I_app', -7.100000;...
    'deepRS', 'gKCa', 0.000000;...
    'deepRS', 'gl', 0.310000;...
    'deepRS', 'Inoise', 0.250000};

model_vary{2} = {'deepRS', 'I_app', -6.500000;...
    'deepRS->deepFS', 'g_SYN', 0.200000;...
    'deepFS->deepRS', 'g_SYN', 0.075000;...
    'deepFS->deepRS', 'tauRx', 0.250000;...
    'deepFS->deepRS', 'tauDx', 2.500000;...
    'deepRS->deepFS', 'tauDx', 50.000000;...
    'deepRS->deepFS', 'tauRx', 0.250000;...
    'deepFS', 'stim', 0.950000;...
    'deepRS', 'gl', 0.270000;...
    'deepRS', 'gKCa', 0.000000;...
    'deepRS', 'Inoise', 0.250000};

model_vary{3} = {'deepRS', 'I_app', -7.600000;...
    'deepRS->deepFS', 'g_SYN', 0.200000;...
    'deepFS->deepRS', 'g_SYN', 0.075000;...
    'deepFS->deepRS', 'tauRx', 0.250000;...
    'deepFS->deepRS', 'tauDx', 2.500000;...
    'deepRS->deepFS', 'tauDx', 50.000000;...
    'deepRS->deepFS', 'tauRx', 0.250000;...
    'deepFS', 'stim', 0.950000;...
    'deepRS', 'gl', 0.780000;...
    'deepRS', 'gKs', 0.000000;...
    'deepRS', 'gKCa', 0.000000;...
    'deepRS', 'Inoise', 0.250000};

model_vary{4} = {'deepRS', 'I_app', -10.500000;...
    'deepRS->deepFS', 'g_SYN', 0.200000;...
    'deepFS->deepRS', 'g_SYN', 0.075000;...
    'deepFS->deepRS', 'tauRx', 0.250000;...
    'deepFS->deepRS', 'tauDx', 2.500000;...
    'deepRS->deepFS', 'tauDx', 50.000000;...
    'deepRS->deepFS', 'tauRx', 0.250000;...
    'deepFS', 'stim', 0.950000;...
    'deepRS', 'gl', 0.780000;...
    'deepRS', 'gKs', 0.000000;...
    'deepRS', 'Inoise', 0.250000};

model_vary{5} = {'deepRS', 'I_app', -9.800000;...
    'deepRS->deepFS', 'g_SYN', 0.200000;...
    'deepFS->deepRS', 'g_SYN', 0.075000;...
    'deepFS->deepRS', 'tauRx', 0.250000;...
    'deepFS->deepRS', 'tauDx', 2.500000;...
    'deepRS->deepFS', 'tauDx', 50.000000;...
    'deepRS->deepFS', 'tauRx', 0.250000;...
    'deepFS', 'stim', 0.950000;...
    'deepRS', 'gl', 0.270000;...
    'deepRS', 'Inoise', 0.250000};

model_vary{6} = {'deepRS', 'I_app', -9.200000;...
    'deepRS', 'Inoise', 0.250000};

end

