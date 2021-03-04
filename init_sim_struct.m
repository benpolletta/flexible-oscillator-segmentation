function sim_struct = init_sim_struct(varargin)

keys = {'tspan', 'sim_mode', 'pulse_mode',...
    'save_data_flag', 'parallel_flag', 'num_cores', 'compile_flag', 'cluster_flag', 'qsub_mode', 'one_solve_file_flag',...
    'include_IB', 'include_deepNG', 'include_deepFS', 'include_deepRS',...
    'include_NG', 'include_RS', 'include_FS', 'include_LTS',...
    'debug_flag', 'analysis_functions', 'analysis_options',...
    };

values = {[0 6000], 1, 0,...
    0, 1, 16, 0, 0, 'array', 1,...
    0, 0, 0, 0,...
    0, 0, 0, 0,...
    0,{},{},...
    };

sim_struct = init_struct(keys, values, varargin{:});