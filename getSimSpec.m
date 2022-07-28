function sim_spec = getSimSpec(NdeepRS, NdeepFS)

i = 0;

% % % % % % % % % % % % % % Deep RS cells % % % % % % % % % % % % % % % % %
    
i = i+1;
spec.populations(i).name = 'deepRS';
spec.populations(i).size = NdeepRS; % 1;
spec.populations(i).equations = {['V''=(I_const+@current)/Cm; V(0)=' num2str(-65) '; monitor functions; monitor V.spikes(0);']};
spec.populations(i).mechanism_list = {'iNaP','iKs','iKsConst','iKDRG','iNaG','iLeak',...
    'CaDynT','iCaT','iKCaT','iPeriodicPulsesBen','itonic', 'iSpeechInput', 'iVariedPulses'}; 
spec.populations(i).parameters = {...
    'Cm', 2.7, 'gKs', 1.4472, 'gNaP', 0.4307, 'gKCa', 0.1512, 'bKCa', .002,...
    'gCa', 0.54, 'CAF', 2.2222,...
    'gl', 0.27, 'gNaG', 135, 'gKDR', 54, 'tau_h', 5, 'tau_m', 5,...
    'I_const', 0, 'ton', 500, 'toff', Inf, 'I_app', -9.2, 'Inoise', 0.25,...    
    'PPstim', 0, 'PPcenter', 0, 'PPnorm', 0, 'PPonset', 10,... 
    'FMPstim', 0, 'STPstim', 0, 'VPstim', 0, 'gSpike', 0, 'gCar', 0, 'gSpeech', 0,... 
    };

if NdeepFS > 0

% % % % % % % % % % % % % % Deep FS cells % % % % % % % % % % % % % % % % %
    
    i=i+1;
    spec.populations(i).name = 'deepFS';
    spec.populations(i).size = NdeepFS; % 1;
    spec.populations(i).equations = {['V''=(@current)/Cm; V(0)=' num2str(-65) '; monitor functions;']};
    spec.populations(i).mechanism_list = {'IBitonic','IBnoise','FSiNaF','FSiKDR','IBleak','iSpeechInput'};
    spec.populations(i).parameters = {...
        'V_IC', -65, 'IC_noise', 0, 'Cm', 0.9, 'E_l', -70, 'g_l', 0.1,...
        'gNaF', 100, 'E_NaF', 50,...
        'gKDR', 80, 'E_KDR', -95,...
        'stim', 0.95, 'onset', 0, 'offset', Inf,...
        'V_noise', 0,...
        };

% % % % % % % % % % % %  Connections  % % % % % % % % % % % % % %

    i=0;

    % % deepRS->deepFS synaptic connection
    i=i+1;
    spec.connections(i).direction = 'deepRS->deepFS';
    spec.connections(i).mechanism_list = {'IBaIBdbiSYNseed'};
    spec.connections(i).parameters = {'g_SYN', 0.075, 'E_SYN', 0,'tauDx', 2.5,...
        'tauRx', 0.25, 'fanout', 0, 'IC_noise', 0, 'g_SYN_hetero', 0,...
        };

    % % deepFS->deepRS Synaptic connections
    i=i+1;
    spec.connections(i).direction = 'deepFS->deepRS';                  
    spec.connections(i).mechanism_list = {'IBaIBdbiSYNseed'};
    spec.connections(i).parameters = {'g_SYN', 0.15, 'E_SYN', -95, 'tauDx', 50,...
        'tauRx', 0.25, 'fanout', 0, 'IC_noise', 0, 'g_SYN_hetero', 0,...
        };
    
end

sim_spec = spec;

end