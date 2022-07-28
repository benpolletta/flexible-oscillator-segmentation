function [data, name, sim_spec] = simulateModel(sim_struct)

if nargin < 1; sim_struct = []; end
if isempty(sim_struct); sim_struct = init_sim_struct; end

%% Setting up files to save simulations.
Today = datestr(datenum(date),'yy-mm-dd');

start_dir = pwd;

savepath = fullfile(pwd, 'Sims', Today);
mkdir(savepath);

Now = clock;
name = sprintf('FOS_%g_%g_%.4g', Now(4), Now(5), Now(6));

%% Setting up simulation specifications.

sim_spec = getSimSpec(sim_struct.include_deepRS, sim_struct.include_deepFS);

switch isfield(sim_struct, 'vary')
    
    case 0
        
        if sim_struct.cluster_flag
            
            sim_struct.one_solve_file_flag = 0;
            
        end
    
        save(fullfile(savepath, [name, '_sim_spec.mat']), 'sim_spec', 'sim_struct', 'name');
        
    case 1
        
        if sim_struct.cluster_flag
            
            sim_struct.one_solve_file_flag = 1;
            
        end
        
        vary = sim_struct.vary;
    
        save(fullfile(savepath, [name, '_sim_spec.mat']), 'sim_spec', 'sim_struct', 'vary', 'name');
        
end

dsSim_struct = rmfield(sim_struct, {'include_deepRS', 'include_deepFS'});
sim_cell = [fieldnames(dsSim_struct), struct2cell(dsSim_struct)]';
    
if sim_struct.cluster_flag
    
    data = dsSimulate(sim_spec, sim_cell{:},'study_dir',fullfile(savepath, name));
    
    cd (start_dir)
    
    return
    
else
    
    tic;
    
    data = dsSimulate(sim_spec, sim_cell{:},'study_dir',fullfile(savepath, name));
    
    toc;
    
end

close('all')

dsPlot(data)

figHandles = findobj('Type', 'Figure');

for f = 1:length(figHandles)

    save_as_pdf(figHandles(f), fullfile(savepath, [name, '_', num2str(f)]))
    
end

cd (start_dir)

end