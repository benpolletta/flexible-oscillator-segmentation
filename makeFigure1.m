function sim_struct = makeFigure1

[model_vary, model_names, include_deepFS] = getModelParams;

for m = 1:length(model_vary)
    
    Iapp=0:.1:20;

    simulation_vary = {'deepRS', 'Iapp', Iapp};
    
    vary = [model_vary{m}; simulation_vary];
    
    sim_struct{m} = initSimStruct('vary', vary, 'name', model_names{m}, 'include_deepFS', include_deepFS(m),...
        'analysis_functions', @spike_metrics);
    
    simulateModel(sim_struct{m})
    
end