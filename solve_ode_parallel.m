function [T,deepFS_V,deepFS_FSiNaF_hNaF,deepFS_FSiKDR_mKDR,deepRS_V,deepRS_iNaP_m,deepRS_iKs_n,deepRS_iKDRG_mKDR,deepRS_iNaG_h,deepRS_CaDynT_cai,deepRS_iCaT_s,deepRS_iKCaT_q,deepFS_deepRS_IBaIBdbiSYNseed_s,deepRS_deepFS_IBaIBdbiSYNseed_s,deepFS_FSiKDR_IKDR,deepFS_FSiKDR_aM,deepFS_FSiKDR_bM,deepFS_FSiKDR_minf,deepFS_FSiKDR_mtau,deepFS_FSiNaF_INaF,deepFS_FSiNaF_aH,deepFS_FSiNaF_bH,deepFS_FSiNaF_hinf,deepFS_FSiNaF_htau,deepFS_FSiNaF_m0,deepFS_IBitonic_Itonic,deepFS_IBleak_Ileak,deepFS_IBnoise_Inoise,deepFS_deepRS_IBaIBdbiSYNseed_ISYN,deepRS_deepFS_IBaIBdbiSYNseed_ISYN,deepRS_iCaT_ICa,deepRS_iCaT_as,deepRS_iCaT_bs,deepRS_iCaT_sinf,deepRS_iCaT_stau,deepRS_iFMPulses_Iext,deepRS_iFMPulses_input,deepRS_iKCaT_IKCa,deepRS_iKCaT_aKCa,deepRS_iKCaT_infKCa,deepRS_iKCaT_tauKCa,deepRS_iKDRG_IKDR,deepRS_iKDRG_aM,deepRS_iKDRG_bM,deepRS_iKs_IKs,deepRS_iKs_ninf,deepRS_iKs_tauKs,deepRS_iLeak_Il,deepRS_iNaG_INa,deepRS_iNaG_ah,deepRS_iNaG_am,deepRS_iNaG_bh,deepRS_iNaG_bm,deepRS_iNaG_minf,deepRS_iNaP_INaP,deepRS_iNaP_minf,deepRS_iPeriodicPulsesBen_Iext,deepRS_iPeriodicPulsesBen_input,deepRS_itonicBen_IBen,deepRS_iKs_tau_div,deepRS_iKDRG_tau_m,deepRS_iNaG_tau_h,deepRS_iPeriodicPulsesBen_PPwidth,deepRS_iPeriodicPulsesBen_Tend,deepRS_iPeriodicPulsesBen_s2,deepRS_itonicBen_toff,deepRS_iFMPulses_FMPwindowlength,deepRS_iFMPulses_Tend,deepRS_iFMPulses_sampling_freq,deepRS_iFMPulses_s2,deepFS_deepRS_IBaIBdbiSYNseed_UB,deepFS_deepRS_IBaIBdbiSYNseed_Xpre,deepFS_deepRS_IBaIBdbiSYNseed_Xpost,deepFS_deepRS_IBaIBdbiSYNseed_mask,deepFS_deepRS_IBaIBdbiSYNseed_gsyn,deepRS_deepFS_IBaIBdbiSYNseed_UB,deepRS_deepFS_IBaIBdbiSYNseed_Xpre,deepRS_deepFS_IBaIBdbiSYNseed_Xpost,deepRS_deepFS_IBaIBdbiSYNseed_mask,deepRS_deepFS_IBaIBdbiSYNseed_gsyn]=solve_ode
% ------------------------------------------------------------
% Parameters:
% ------------------------------------------------------------
params = load('params.mat','p');
p = params.p;
downsample_factor=p.downsample_factor;
dt=p.dt;
T=(p.tspan(1):dt:p.tspan(2))';
ntime=length(T);
nsamp=length(1:downsample_factor:ntime);
% ------------------------------------------------------------
% Fixed variables:
% ------------------------------------------------------------
deepRS_iKs_tau_div =  3^((p.deepRS_iKs_temp-22)/10);
deepRS_iKDRG_tau_m = +5.00000000;
deepRS_iNaG_tau_h = +5.00000000;
deepRS_iPeriodicPulsesBen_PPwidth =  (1000/p.deepRS_iPeriodicPulsesBen_PPfreq)*p.deepRS_iPeriodicPulsesBen_PPduty;
deepRS_iPeriodicPulsesBen_Tend = T(end);
deepRS_iPeriodicPulsesBen_s2 =  getPeriodicPulseFastBensVersion(p.deepRS_iPeriodicPulsesBen_PPfreq,deepRS_iPeriodicPulsesBen_PPwidth,p.deepRS_iPeriodicPulsesBen_PPshift,deepRS_iPeriodicPulsesBen_Tend,p.deepRS_iPeriodicPulsesBen_dt,p.deepRS_iPeriodicPulsesBen_PPonset,p.deepRS_iPeriodicPulsesBen_PPoffset,p.deepRS_Npop,p.deepRS_iPeriodicPulsesBen_kernel_type,p.deepRS_iPeriodicPulsesBen_width2_rise,p.deepRS_iPeriodicPulsesBen_PPcenter,p.deepRS_iPeriodicPulsesBen_PPnorm,0);
deepRS_itonicBen_toff = +20000.00000000;
deepRS_iFMPulses_FMPwindowlength =  (1000/mean([p.deepRS_iFMPulses_FMPlowfreq p.deepRS_iFMPulses_FMPhighfreq]));
deepRS_iFMPulses_Tend = T(end);
deepRS_iFMPulses_sampling_freq =  1000/p.deepRS_iFMPulses_dt;
deepRS_iFMPulses_s2 =  phase_modulated(p.deepRS_iFMPulses_dt, deepRS_iFMPulses_Tend, deepRS_iFMPulses_sampling_freq, [p.deepRS_iFMPulses_FMPlowfreq p.deepRS_iFMPulses_FMPhighfreq], deepRS_iFMPulses_FMPwindowlength, 1);
deepFS_deepRS_IBaIBdbiSYNseed_UB =  max(p.deepRS_Npop,p.deepFS_Npop);
deepFS_deepRS_IBaIBdbiSYNseed_Xpre =  linspace(1,deepFS_deepRS_IBaIBdbiSYNseed_UB,p.deepRS_Npop)'*ones(1,p.deepFS_Npop);
deepFS_deepRS_IBaIBdbiSYNseed_Xpost =  (linspace(1,deepFS_deepRS_IBaIBdbiSYNseed_UB,p.deepFS_Npop)'*ones(1,p.deepRS_Npop))';
deepFS_deepRS_IBaIBdbiSYNseed_mask =  (abs(deepFS_deepRS_IBaIBdbiSYNseed_Xpre-deepFS_deepRS_IBaIBdbiSYNseed_Xpost)<=p.deepFS_deepRS_IBaIBdbiSYNseed_fanout);
deepFS_deepRS_IBaIBdbiSYNseed_gsyn = unifrnd(p.deepFS_deepRS_IBaIBdbiSYNseed_g_SYN-p.deepFS_deepRS_IBaIBdbiSYNseed_gSYNhetero,p.deepFS_deepRS_IBaIBdbiSYNseed_g_SYN+p.deepFS_deepRS_IBaIBdbiSYNseed_gSYNhetero,[p.deepFS_Npop 1])';
deepRS_deepFS_IBaIBdbiSYNseed_UB =  max(p.deepFS_Npop,p.deepRS_Npop);
deepRS_deepFS_IBaIBdbiSYNseed_Xpre =  linspace(1,deepRS_deepFS_IBaIBdbiSYNseed_UB,p.deepFS_Npop)'*ones(1,p.deepRS_Npop);
deepRS_deepFS_IBaIBdbiSYNseed_Xpost =  (linspace(1,deepRS_deepFS_IBaIBdbiSYNseed_UB,p.deepRS_Npop)'*ones(1,p.deepFS_Npop))';
deepRS_deepFS_IBaIBdbiSYNseed_mask =  (abs(deepRS_deepFS_IBaIBdbiSYNseed_Xpre-deepRS_deepFS_IBaIBdbiSYNseed_Xpost)<=p.deepRS_deepFS_IBaIBdbiSYNseed_fanout);
deepRS_deepFS_IBaIBdbiSYNseed_gsyn = unifrnd(p.deepRS_deepFS_IBaIBdbiSYNseed_g_SYN-p.deepRS_deepFS_IBaIBdbiSYNseed_gSYNhetero,p.deepRS_deepFS_IBaIBdbiSYNseed_g_SYN+p.deepRS_deepFS_IBaIBdbiSYNseed_gSYNhetero,[p.deepRS_Npop 1])';
% ------------------------------------------------------------
% Initial conditions:
% ------------------------------------------------------------
% seed the random number generator
rng(p.random_seed);
t=0; k=1;

% STATE_VARIABLES:
deepFS_V_last = -65*ones(1,p.deepFS_Npop);
deepFS_V = zeros(nsamp,p.deepFS_Npop);
  deepFS_V(1,:) = deepFS_V_last;
deepFS_FSiNaF_hNaF_last =  p.deepFS_FSiNaF_IC+p.deepFS_FSiNaF_IC_noise.*rand(1,p.deepFS_Npop);
deepFS_FSiNaF_hNaF = zeros(nsamp,p.deepFS_Npop);
  deepFS_FSiNaF_hNaF(1,:) = deepFS_FSiNaF_hNaF_last;
deepFS_FSiKDR_mKDR_last =  p.deepFS_FSiKDR_IC+p.deepFS_FSiKDR_IC_noise.*rand(1,p.deepFS_Npop);
deepFS_FSiKDR_mKDR = zeros(nsamp,p.deepFS_Npop);
  deepFS_FSiKDR_mKDR(1,:) = deepFS_FSiKDR_mKDR_last;
deepRS_V_last = -65*ones(1,p.deepRS_Npop);
deepRS_V = zeros(nsamp,p.deepRS_Npop);
  deepRS_V(1,:) = deepRS_V_last;
deepRS_iNaP_m_last = p.deepRS_iNaP_m_IC+p.deepRS_iNaP_IC_noise*rand(1,p.deepRS_Npop);
deepRS_iNaP_m = zeros(nsamp,p.deepRS_Npop);
  deepRS_iNaP_m(1,:) = deepRS_iNaP_m_last;
deepRS_iKs_n_last = p.deepRS_iKs_IC+p.deepRS_iKs_IC_noise.*rand(1,p.deepRS_Npop);
deepRS_iKs_n = zeros(nsamp,p.deepRS_Npop);
  deepRS_iKs_n(1,:) = deepRS_iKs_n_last;
deepRS_iKDRG_mKDR_last =  p.deepRS_iKDRG_IC+p.deepRS_iKDRG_IC_noise.*rand(1,p.deepRS_Npop);
deepRS_iKDRG_mKDR = zeros(nsamp,p.deepRS_Npop);
  deepRS_iKDRG_mKDR(1,:) = deepRS_iKDRG_mKDR_last;
deepRS_iNaG_h_last =  p.deepRS_iNaG_IC+p.deepRS_iNaG_IC_noise.*rand(1,p.deepRS_Npop);
deepRS_iNaG_h = zeros(nsamp,p.deepRS_Npop);
  deepRS_iNaG_h(1,:) = deepRS_iNaG_h_last;
deepRS_CaDynT_cai_last =  p.deepRS_CaDynT_cainf + p.deepRS_CaDynT_IC_noise.*rand(1,p.deepRS_Npop);
deepRS_CaDynT_cai = zeros(nsamp,p.deepRS_Npop);
  deepRS_CaDynT_cai(1,:) = deepRS_CaDynT_cai_last;
deepRS_iCaT_s_last =  (( 1.6./(1+exp(-.072*(-65-p.deepRS_iCaT_Caoffset-65)))))/((( 1.6./(1+exp(-.072*(-65-p.deepRS_iCaT_Caoffset-65)))))+(( 0.02*(-65-p.deepRS_iCaT_Caoffset-51.1)./(exp((-65-p.deepRS_iCaT_Caoffset-51.1)/5)-1))))+p.deepRS_iCaT_IC_noise.*rand(1,p.deepRS_Npop);
deepRS_iCaT_s = zeros(nsamp,p.deepRS_Npop);
  deepRS_iCaT_s(1,:) = deepRS_iCaT_s_last;
deepRS_iKCaT_q_last =  (( (( min(p.deepRS_iKCaT_aKCa_factor.*0,1)))./((( min(p.deepRS_iKCaT_aKCa_factor.*0,1)))+p.deepRS_iKCaT_bKCa)))+p.deepRS_iKCaT_IC_noise.*rand(1,p.deepRS_Npop);
deepRS_iKCaT_q = zeros(nsamp,p.deepRS_Npop);
  deepRS_iKCaT_q(1,:) = deepRS_iKCaT_q_last;
deepFS_deepRS_IBaIBdbiSYNseed_s_last =  p.deepFS_deepRS_IBaIBdbiSYNseed_IC*ones(1,p.deepRS_Npop)+p.deepFS_deepRS_IBaIBdbiSYNseed_IC_noise.*rand(1,p.deepRS_Npop);
deepFS_deepRS_IBaIBdbiSYNseed_s = zeros(nsamp,p.deepRS_Npop);
  deepFS_deepRS_IBaIBdbiSYNseed_s(1,:) = deepFS_deepRS_IBaIBdbiSYNseed_s_last;
deepRS_deepFS_IBaIBdbiSYNseed_s_last =  p.deepRS_deepFS_IBaIBdbiSYNseed_IC*ones(1,p.deepFS_Npop)+p.deepRS_deepFS_IBaIBdbiSYNseed_IC_noise.*rand(1,p.deepFS_Npop);
deepRS_deepFS_IBaIBdbiSYNseed_s = zeros(nsamp,p.deepFS_Npop);
  deepRS_deepFS_IBaIBdbiSYNseed_s(1,:) = deepRS_deepFS_IBaIBdbiSYNseed_s_last;

% MONITORS:
  deepFS_FSiKDR_IKDR_last=p.deepFS_FSiKDR_gKDR.*deepFS_FSiKDR_mKDR_last.^4.*(deepFS_V_last-p.deepFS_FSiKDR_E_KDR);
deepFS_FSiKDR_IKDR = zeros(nsamp,p.deepFS_Npop);
  deepFS_FSiKDR_IKDR(1,:)=p.deepFS_FSiKDR_gKDR.*deepFS_FSiKDR_mKDR(1,:).^4.*(deepFS_V(1,:)-p.deepFS_FSiKDR_E_KDR);
  deepFS_FSiKDR_aM_last=(( 1./(1+exp((-deepFS_V_last-p.deepFS_FSiKDR_KDR_V1+p.deepFS_FSiKDR_KDR_offset)/p.deepFS_FSiKDR_KDR_d1)))) ./ (( .25+4.35*exp(-abs(deepFS_V_last+p.deepFS_FSiKDR_KDR_V2-p.deepFS_FSiKDR_KDR_offset)/p.deepFS_FSiKDR_KDR_d2)));
deepFS_FSiKDR_aM = zeros(nsamp,p.deepFS_Npop);
  deepFS_FSiKDR_aM(1,:)=(( 1./(1+exp((-deepFS_V(1,:)-p.deepFS_FSiKDR_KDR_V1+p.deepFS_FSiKDR_KDR_offset)/p.deepFS_FSiKDR_KDR_d1)))) ./ (( .25+4.35*exp(-abs(deepFS_V(1,:)+p.deepFS_FSiKDR_KDR_V2-p.deepFS_FSiKDR_KDR_offset)/p.deepFS_FSiKDR_KDR_d2)));
  deepFS_FSiKDR_bM_last=(1-(( 1./(1+exp((-deepFS_V_last-p.deepFS_FSiKDR_KDR_V1+p.deepFS_FSiKDR_KDR_offset)/p.deepFS_FSiKDR_KDR_d1)))))./(( .25+4.35*exp(-abs(deepFS_V_last+p.deepFS_FSiKDR_KDR_V2-p.deepFS_FSiKDR_KDR_offset)/p.deepFS_FSiKDR_KDR_d2)));
deepFS_FSiKDR_bM = zeros(nsamp,p.deepFS_Npop);
  deepFS_FSiKDR_bM(1,:)=(1-(( 1./(1+exp((-deepFS_V(1,:)-p.deepFS_FSiKDR_KDR_V1+p.deepFS_FSiKDR_KDR_offset)/p.deepFS_FSiKDR_KDR_d1)))))./(( .25+4.35*exp(-abs(deepFS_V(1,:)+p.deepFS_FSiKDR_KDR_V2-p.deepFS_FSiKDR_KDR_offset)/p.deepFS_FSiKDR_KDR_d2)));
  deepFS_FSiKDR_minf_last=1./(1+exp((-deepFS_V_last-p.deepFS_FSiKDR_KDR_V1+p.deepFS_FSiKDR_KDR_offset)/p.deepFS_FSiKDR_KDR_d1));
deepFS_FSiKDR_minf = zeros(nsamp,p.deepFS_Npop);
  deepFS_FSiKDR_minf(1,:)=1./(1+exp((-deepFS_V(1,:)-p.deepFS_FSiKDR_KDR_V1+p.deepFS_FSiKDR_KDR_offset)/p.deepFS_FSiKDR_KDR_d1));
  deepFS_FSiKDR_mtau_last=.25+4.35*exp(-abs(deepFS_V_last+p.deepFS_FSiKDR_KDR_V2-p.deepFS_FSiKDR_KDR_offset)/p.deepFS_FSiKDR_KDR_d2);
deepFS_FSiKDR_mtau = zeros(nsamp,p.deepFS_Npop);
  deepFS_FSiKDR_mtau(1,:)=.25+4.35*exp(-abs(deepFS_V(1,:)+p.deepFS_FSiKDR_KDR_V2-p.deepFS_FSiKDR_KDR_offset)/p.deepFS_FSiKDR_KDR_d2);
  deepFS_FSiNaF_INaF_last=p.deepFS_FSiNaF_gNaF.*(( 1./(1+exp((-deepFS_V_last-p.deepFS_FSiNaF_NaF_V0+p.deepFS_FSiNaF_NaF_offset)/10)))).^3.*deepFS_FSiNaF_hNaF_last.*(deepFS_V_last-p.deepFS_FSiNaF_E_NaF);
deepFS_FSiNaF_INaF = zeros(nsamp,p.deepFS_Npop);
  deepFS_FSiNaF_INaF(1,:)=p.deepFS_FSiNaF_gNaF.*(( 1./(1+exp((-deepFS_V(1,:)-p.deepFS_FSiNaF_NaF_V0+p.deepFS_FSiNaF_NaF_offset)/10)))).^3.*deepFS_FSiNaF_hNaF(1,:).*(deepFS_V(1,:)-p.deepFS_FSiNaF_E_NaF);
  deepFS_FSiNaF_aH_last=(( 1./(1+exp((deepFS_V_last+p.deepFS_FSiNaF_NaF_V1-p.deepFS_FSiNaF_NaF_offset)/p.deepFS_FSiNaF_NaF_d1)))) ./ (( p.deepFS_FSiNaF_NaF_c0 + p.deepFS_FSiNaF_NaF_c1./(1+exp((deepFS_V_last+p.deepFS_FSiNaF_NaF_V2-p.deepFS_FSiNaF_NaF_offset)/p.deepFS_FSiNaF_NaF_d2))));
deepFS_FSiNaF_aH = zeros(nsamp,p.deepFS_Npop);
  deepFS_FSiNaF_aH(1,:)=(( 1./(1+exp((deepFS_V(1,:)+p.deepFS_FSiNaF_NaF_V1-p.deepFS_FSiNaF_NaF_offset)/p.deepFS_FSiNaF_NaF_d1)))) ./ (( p.deepFS_FSiNaF_NaF_c0 + p.deepFS_FSiNaF_NaF_c1./(1+exp((deepFS_V(1,:)+p.deepFS_FSiNaF_NaF_V2-p.deepFS_FSiNaF_NaF_offset)/p.deepFS_FSiNaF_NaF_d2))));
  deepFS_FSiNaF_bH_last=(1-(( 1./(1+exp((deepFS_V_last+p.deepFS_FSiNaF_NaF_V1-p.deepFS_FSiNaF_NaF_offset)/p.deepFS_FSiNaF_NaF_d1)))))./(( p.deepFS_FSiNaF_NaF_c0 + p.deepFS_FSiNaF_NaF_c1./(1+exp((deepFS_V_last+p.deepFS_FSiNaF_NaF_V2-p.deepFS_FSiNaF_NaF_offset)/p.deepFS_FSiNaF_NaF_d2))));
deepFS_FSiNaF_bH = zeros(nsamp,p.deepFS_Npop);
  deepFS_FSiNaF_bH(1,:)=(1-(( 1./(1+exp((deepFS_V(1,:)+p.deepFS_FSiNaF_NaF_V1-p.deepFS_FSiNaF_NaF_offset)/p.deepFS_FSiNaF_NaF_d1)))))./(( p.deepFS_FSiNaF_NaF_c0 + p.deepFS_FSiNaF_NaF_c1./(1+exp((deepFS_V(1,:)+p.deepFS_FSiNaF_NaF_V2-p.deepFS_FSiNaF_NaF_offset)/p.deepFS_FSiNaF_NaF_d2))));
  deepFS_FSiNaF_hinf_last=1./(1+exp((deepFS_V_last+p.deepFS_FSiNaF_NaF_V1-p.deepFS_FSiNaF_NaF_offset)/p.deepFS_FSiNaF_NaF_d1));
deepFS_FSiNaF_hinf = zeros(nsamp,p.deepFS_Npop);
  deepFS_FSiNaF_hinf(1,:)=1./(1+exp((deepFS_V(1,:)+p.deepFS_FSiNaF_NaF_V1-p.deepFS_FSiNaF_NaF_offset)/p.deepFS_FSiNaF_NaF_d1));
  deepFS_FSiNaF_htau_last=p.deepFS_FSiNaF_NaF_c0 + p.deepFS_FSiNaF_NaF_c1./(1+exp((deepFS_V_last+p.deepFS_FSiNaF_NaF_V2-p.deepFS_FSiNaF_NaF_offset)/p.deepFS_FSiNaF_NaF_d2));
deepFS_FSiNaF_htau = zeros(nsamp,p.deepFS_Npop);
  deepFS_FSiNaF_htau(1,:)=p.deepFS_FSiNaF_NaF_c0 + p.deepFS_FSiNaF_NaF_c1./(1+exp((deepFS_V(1,:)+p.deepFS_FSiNaF_NaF_V2-p.deepFS_FSiNaF_NaF_offset)/p.deepFS_FSiNaF_NaF_d2));
  deepFS_FSiNaF_m0_last=1./(1+exp((-deepFS_V_last-p.deepFS_FSiNaF_NaF_V0+p.deepFS_FSiNaF_NaF_offset)/10));
deepFS_FSiNaF_m0 = zeros(nsamp,p.deepFS_Npop);
  deepFS_FSiNaF_m0(1,:)=1./(1+exp((-deepFS_V(1,:)-p.deepFS_FSiNaF_NaF_V0+p.deepFS_FSiNaF_NaF_offset)/10));
  deepFS_IBitonic_Itonic_last=p.deepFS_IBitonic_stim*(t>p.deepFS_IBitonic_onset & t<p.deepFS_IBitonic_offset);
deepFS_IBitonic_Itonic = zeros(nsamp,p.deepFS_Npop);
  deepFS_IBitonic_Itonic(1,:)=p.deepFS_IBitonic_stim*(t>p.deepFS_IBitonic_onset & t<p.deepFS_IBitonic_offset);
  deepFS_IBleak_Ileak_last=p.deepFS_IBleak_g_l.*(deepFS_V_last-p.deepFS_IBleak_E_l);
deepFS_IBleak_Ileak = zeros(nsamp,p.deepFS_Npop);
  deepFS_IBleak_Ileak(1,:)=p.deepFS_IBleak_g_l.*(deepFS_V(1,:)-p.deepFS_IBleak_E_l);
  deepFS_IBnoise_Inoise_last=p.deepFS_IBnoise_V_noise.*randn(1,p.deepFS_Npop);
deepFS_IBnoise_Inoise = zeros(nsamp,p.deepFS_Npop);
  deepFS_IBnoise_Inoise(1,:)=p.deepFS_IBnoise_V_noise.*randn(1,p.deepFS_Npop);
  deepFS_deepRS_IBaIBdbiSYNseed_ISYN_last=(deepFS_deepRS_IBaIBdbiSYNseed_gsyn.*(deepFS_deepRS_IBaIBdbiSYNseed_s_last*deepFS_deepRS_IBaIBdbiSYNseed_mask).*(deepFS_V_last-p.deepFS_deepRS_IBaIBdbiSYNseed_E_SYN));
deepFS_deepRS_IBaIBdbiSYNseed_ISYN = zeros(nsamp,p.deepFS_Npop);
  deepFS_deepRS_IBaIBdbiSYNseed_ISYN(1,:)=(deepFS_deepRS_IBaIBdbiSYNseed_gsyn.*(deepFS_deepRS_IBaIBdbiSYNseed_s(1,:)*deepFS_deepRS_IBaIBdbiSYNseed_mask).*(deepFS_V(1,:)-p.deepFS_deepRS_IBaIBdbiSYNseed_E_SYN));
  deepRS_deepFS_IBaIBdbiSYNseed_ISYN_last=(deepRS_deepFS_IBaIBdbiSYNseed_gsyn.*(deepRS_deepFS_IBaIBdbiSYNseed_s_last*deepRS_deepFS_IBaIBdbiSYNseed_mask).*(deepRS_V_last-p.deepRS_deepFS_IBaIBdbiSYNseed_E_SYN));
deepRS_deepFS_IBaIBdbiSYNseed_ISYN = zeros(nsamp,p.deepRS_Npop);
  deepRS_deepFS_IBaIBdbiSYNseed_ISYN(1,:)=(deepRS_deepFS_IBaIBdbiSYNseed_gsyn.*(deepRS_deepFS_IBaIBdbiSYNseed_s(1,:)*deepRS_deepFS_IBaIBdbiSYNseed_mask).*(deepRS_V(1,:)-p.deepRS_deepFS_IBaIBdbiSYNseed_E_SYN));
  deepRS_iCaT_ICa_last=p.deepRS_iCaT_gCa.*(deepRS_iCaT_s_last.^2).*(deepRS_V_last-p.deepRS_iCaT_ECa);
deepRS_iCaT_ICa = zeros(nsamp,p.deepRS_Npop);
  deepRS_iCaT_ICa(1,:)=p.deepRS_iCaT_gCa.*(deepRS_iCaT_s(1,:).^2).*(deepRS_V(1,:)-p.deepRS_iCaT_ECa);
  deepRS_iCaT_as_last=1.6./(1+exp(-.072*(deepRS_V_last-p.deepRS_iCaT_Caoffset-65)));
deepRS_iCaT_as = zeros(nsamp,p.deepRS_Npop);
  deepRS_iCaT_as(1,:)=1.6./(1+exp(-.072*(deepRS_V(1,:)-p.deepRS_iCaT_Caoffset-65)));
  deepRS_iCaT_bs_last=0.02*(deepRS_V_last-p.deepRS_iCaT_Caoffset-51.1)./(exp((deepRS_V_last-p.deepRS_iCaT_Caoffset-51.1)/5)-1);
deepRS_iCaT_bs = zeros(nsamp,p.deepRS_Npop);
  deepRS_iCaT_bs(1,:)=0.02*(deepRS_V(1,:)-p.deepRS_iCaT_Caoffset-51.1)./(exp((deepRS_V(1,:)-p.deepRS_iCaT_Caoffset-51.1)/5)-1);
  deepRS_iCaT_sinf_last=(( 1.6./(1+exp(-.072*(deepRS_V_last-p.deepRS_iCaT_Caoffset-65)))))./((( 1.6./(1+exp(-.072*(deepRS_V_last-p.deepRS_iCaT_Caoffset-65)))))+(( 0.02*(deepRS_V_last-p.deepRS_iCaT_Caoffset-51.1)./(exp((deepRS_V_last-p.deepRS_iCaT_Caoffset-51.1)/5)-1))));
deepRS_iCaT_sinf = zeros(nsamp,p.deepRS_Npop);
  deepRS_iCaT_sinf(1,:)=(( 1.6./(1+exp(-.072*(deepRS_V(1,:)-p.deepRS_iCaT_Caoffset-65)))))./((( 1.6./(1+exp(-.072*(deepRS_V(1,:)-p.deepRS_iCaT_Caoffset-65)))))+(( 0.02*(deepRS_V(1,:)-p.deepRS_iCaT_Caoffset-51.1)./(exp((deepRS_V(1,:)-p.deepRS_iCaT_Caoffset-51.1)/5)-1))));
  deepRS_iCaT_stau_last=1./((( 1.6./(1+exp(-.072*(deepRS_V_last-p.deepRS_iCaT_Caoffset-65)))))+(( 0.02*(deepRS_V_last-p.deepRS_iCaT_Caoffset-51.1)./(exp((deepRS_V_last-p.deepRS_iCaT_Caoffset-51.1)/5)-1))));
deepRS_iCaT_stau = zeros(nsamp,p.deepRS_Npop);
  deepRS_iCaT_stau(1,:)=1./((( 1.6./(1+exp(-.072*(deepRS_V(1,:)-p.deepRS_iCaT_Caoffset-65)))))+(( 0.02*(deepRS_V(1,:)-p.deepRS_iCaT_Caoffset-51.1)./(exp((deepRS_V(1,:)-p.deepRS_iCaT_Caoffset-51.1)/5)-1))));
  deepRS_iFMPulses_Iext_last=p.deepRS_iFMPulses_FMPstim*(( deepRS_iFMPulses_s2(k,:)));
deepRS_iFMPulses_Iext = zeros(nsamp,p.deepRS_Npop);
  deepRS_iFMPulses_Iext(1,:)=p.deepRS_iFMPulses_FMPstim*(( deepRS_iFMPulses_s2(k,:)));
  deepRS_iFMPulses_input_last=deepRS_iFMPulses_s2(k,:);
deepRS_iFMPulses_input = zeros(nsamp,p.deepRS_Npop);
  deepRS_iFMPulses_input(1,:)=deepRS_iFMPulses_s2(k,:);
  deepRS_iKCaT_IKCa_last=p.deepRS_iKCaT_gKCa*deepRS_iKCaT_q_last*(deepRS_V_last-p.deepRS_iKCaT_EK);
deepRS_iKCaT_IKCa = zeros(nsamp,p.deepRS_Npop);
  deepRS_iKCaT_IKCa(1,:)=p.deepRS_iKCaT_gKCa*deepRS_iKCaT_q(1,:)*(deepRS_V(1,:)-p.deepRS_iKCaT_EK);
  deepRS_iKCaT_aKCa_last=min(p.deepRS_iKCaT_aKCa_factor.*deepRS_CaDynT_cai_last,1);
deepRS_iKCaT_aKCa = zeros(nsamp,p.deepRS_Npop);
  deepRS_iKCaT_aKCa(1,:)=min(p.deepRS_iKCaT_aKCa_factor.*deepRS_CaDynT_cai(1,:),1);
  deepRS_iKCaT_infKCa_last=(( min(p.deepRS_iKCaT_aKCa_factor.*deepRS_CaDynT_cai_last,1)))./((( min(p.deepRS_iKCaT_aKCa_factor.*deepRS_CaDynT_cai_last,1)))+p.deepRS_iKCaT_bKCa);
deepRS_iKCaT_infKCa = zeros(nsamp,p.deepRS_Npop);
  deepRS_iKCaT_infKCa(1,:)=(( min(p.deepRS_iKCaT_aKCa_factor.*deepRS_CaDynT_cai(1,:),1)))./((( min(p.deepRS_iKCaT_aKCa_factor.*deepRS_CaDynT_cai(1,:),1)))+p.deepRS_iKCaT_bKCa);
  deepRS_iKCaT_tauKCa_last=1./((( min(p.deepRS_iKCaT_aKCa_factor.*deepRS_CaDynT_cai_last,1)))+p.deepRS_iKCaT_bKCa);
deepRS_iKCaT_tauKCa = zeros(nsamp,p.deepRS_Npop);
  deepRS_iKCaT_tauKCa(1,:)=1./((( min(p.deepRS_iKCaT_aKCa_factor.*deepRS_CaDynT_cai(1,:),1)))+p.deepRS_iKCaT_bKCa);
  deepRS_iKDRG_IKDR_last=(p.deepRS_iKDRG_gKDR/p.deepRS_iKDRG_gdenomK).*deepRS_iKDRG_mKDR_last.^4.*(deepRS_V_last-p.deepRS_iKDRG_E_KDR);
deepRS_iKDRG_IKDR = zeros(nsamp,p.deepRS_Npop);
  deepRS_iKDRG_IKDR(1,:)=(p.deepRS_iKDRG_gKDR/p.deepRS_iKDRG_gdenomK).*deepRS_iKDRG_mKDR(1,:).^4.*(deepRS_V(1,:)-p.deepRS_iKDRG_E_KDR);
  deepRS_iKDRG_aM_last=-.01*(deepRS_V_last+20-p.deepRS_iKDRG_Koffset)./(exp(-(deepRS_V_last+20-p.deepRS_iKDRG_Koffset)/10)-1);
deepRS_iKDRG_aM = zeros(nsamp,p.deepRS_Npop);
  deepRS_iKDRG_aM(1,:)=-.01*(deepRS_V(1,:)+20-p.deepRS_iKDRG_Koffset)./(exp(-(deepRS_V(1,:)+20-p.deepRS_iKDRG_Koffset)/10)-1);
  deepRS_iKDRG_bM_last=.125*exp(-(deepRS_V_last+30-p.deepRS_iKDRG_Koffset)/80);
deepRS_iKDRG_bM = zeros(nsamp,p.deepRS_Npop);
  deepRS_iKDRG_bM(1,:)=.125*exp(-(deepRS_V(1,:)+30-p.deepRS_iKDRG_Koffset)/80);
  deepRS_iKs_IKs_last=p.deepRS_iKs_gKs.*deepRS_iKs_n_last.*(deepRS_V_last-p.deepRS_iKs_E_Ks);
deepRS_iKs_IKs = zeros(nsamp,p.deepRS_Npop);
  deepRS_iKs_IKs(1,:)=p.deepRS_iKs_gKs.*deepRS_iKs_n(1,:).*(deepRS_V(1,:)-p.deepRS_iKs_E_Ks);
  deepRS_iKs_ninf_last=1./(1+exp(-(deepRS_V_last-(p.deepRS_iKs_halfKs+p.deepRS_iKs_Ks_offset))/10));
deepRS_iKs_ninf = zeros(nsamp,p.deepRS_Npop);
  deepRS_iKs_ninf(1,:)=1./(1+exp(-(deepRS_V(1,:)-(p.deepRS_iKs_halfKs+p.deepRS_iKs_Ks_offset))/10));
  deepRS_iKs_tauKs_last=p.deepRS_iKs_tau_mult*1000./(3.3*deepRS_iKs_tau_div*(exp((deepRS_V_last-(p.deepRS_iKs_halfKs+p.deepRS_iKs_Ks_offset))/40)+exp(-(deepRS_V_last-(p.deepRS_iKs_halfKs+p.deepRS_iKs_Ks_offset))/20)));
deepRS_iKs_tauKs = zeros(nsamp,p.deepRS_Npop);
  deepRS_iKs_tauKs(1,:)=p.deepRS_iKs_tau_mult*1000./(3.3*deepRS_iKs_tau_div*(exp((deepRS_V(1,:)-(p.deepRS_iKs_halfKs+p.deepRS_iKs_Ks_offset))/40)+exp(-(deepRS_V(1,:)-(p.deepRS_iKs_halfKs+p.deepRS_iKs_Ks_offset))/20)));
  deepRS_iLeak_Il_last=(p.deepRS_iLeak_gl/p.deepRS_iLeak_gdenoml).*(deepRS_V_last-p.deepRS_iLeak_El);
deepRS_iLeak_Il = zeros(nsamp,p.deepRS_Npop);
  deepRS_iLeak_Il(1,:)=(p.deepRS_iLeak_gl/p.deepRS_iLeak_gdenoml).*(deepRS_V(1,:)-p.deepRS_iLeak_El);
  deepRS_iNaG_INa_last=(p.deepRS_iNaG_gNa/p.deepRS_iNaG_gdenomN)*((( (( -(deepRS_V_last+16-p.deepRS_iNaG_Noffset)./(10*(exp(-(deepRS_V_last+16-p.deepRS_iNaG_Noffset)/10)-1))))./((( -(deepRS_V_last+16-p.deepRS_iNaG_Noffset)./(10*(exp(-(deepRS_V_last+16-p.deepRS_iNaG_Noffset)/10)-1))))+(( 4*exp(-(deepRS_V_last+41-p.deepRS_iNaG_Noffset)/18)))))).^3).*deepRS_iNaG_h_last.*(deepRS_V_last-p.deepRS_iNaG_ENa);
deepRS_iNaG_INa = zeros(nsamp,p.deepRS_Npop);
  deepRS_iNaG_INa(1,:)=(p.deepRS_iNaG_gNa/p.deepRS_iNaG_gdenomN)*((( (( -(deepRS_V(1,:)+16-p.deepRS_iNaG_Noffset)./(10*(exp(-(deepRS_V(1,:)+16-p.deepRS_iNaG_Noffset)/10)-1))))./((( -(deepRS_V(1,:)+16-p.deepRS_iNaG_Noffset)./(10*(exp(-(deepRS_V(1,:)+16-p.deepRS_iNaG_Noffset)/10)-1))))+(( 4*exp(-(deepRS_V(1,:)+41-p.deepRS_iNaG_Noffset)/18)))))).^3).*deepRS_iNaG_h(1,:).*(deepRS_V(1,:)-p.deepRS_iNaG_ENa);
  deepRS_iNaG_ah_last=.07*exp(-(deepRS_V_last+30-p.deepRS_iNaG_Noffset)/20);
deepRS_iNaG_ah = zeros(nsamp,p.deepRS_Npop);
  deepRS_iNaG_ah(1,:)=.07*exp(-(deepRS_V(1,:)+30-p.deepRS_iNaG_Noffset)/20);
  deepRS_iNaG_am_last=-(deepRS_V_last+16-p.deepRS_iNaG_Noffset)./(10*(exp(-(deepRS_V_last+16-p.deepRS_iNaG_Noffset)/10)-1));
deepRS_iNaG_am = zeros(nsamp,p.deepRS_Npop);
  deepRS_iNaG_am(1,:)=-(deepRS_V(1,:)+16-p.deepRS_iNaG_Noffset)./(10*(exp(-(deepRS_V(1,:)+16-p.deepRS_iNaG_Noffset)/10)-1));
  deepRS_iNaG_bh_last=1./(exp(-(deepRS_V_last-p.deepRS_iNaG_Noffset)/10)+1);
deepRS_iNaG_bh = zeros(nsamp,p.deepRS_Npop);
  deepRS_iNaG_bh(1,:)=1./(exp(-(deepRS_V(1,:)-p.deepRS_iNaG_Noffset)/10)+1);
  deepRS_iNaG_bm_last=4*exp(-(deepRS_V_last+41-p.deepRS_iNaG_Noffset)/18);
deepRS_iNaG_bm = zeros(nsamp,p.deepRS_Npop);
  deepRS_iNaG_bm(1,:)=4*exp(-(deepRS_V(1,:)+41-p.deepRS_iNaG_Noffset)/18);
  deepRS_iNaG_minf_last=(( -(deepRS_V_last+16-p.deepRS_iNaG_Noffset)./(10*(exp(-(deepRS_V_last+16-p.deepRS_iNaG_Noffset)/10)-1))))./((( -(deepRS_V_last+16-p.deepRS_iNaG_Noffset)./(10*(exp(-(deepRS_V_last+16-p.deepRS_iNaG_Noffset)/10)-1))))+(( 4*exp(-(deepRS_V_last+41-p.deepRS_iNaG_Noffset)/18))));
deepRS_iNaG_minf = zeros(nsamp,p.deepRS_Npop);
  deepRS_iNaG_minf(1,:)=(( -(deepRS_V(1,:)+16-p.deepRS_iNaG_Noffset)./(10*(exp(-(deepRS_V(1,:)+16-p.deepRS_iNaG_Noffset)/10)-1))))./((( -(deepRS_V(1,:)+16-p.deepRS_iNaG_Noffset)./(10*(exp(-(deepRS_V(1,:)+16-p.deepRS_iNaG_Noffset)/10)-1))))+(( 4*exp(-(deepRS_V(1,:)+41-p.deepRS_iNaG_Noffset)/18))));
  deepRS_iNaP_INaP_last=p.deepRS_iNaP_gNaP.*deepRS_iNaP_m_last.*(deepRS_V_last-p.deepRS_iNaP_ENa);
deepRS_iNaP_INaP = zeros(nsamp,p.deepRS_Npop);
  deepRS_iNaP_INaP(1,:)=p.deepRS_iNaP_gNaP.*deepRS_iNaP_m(1,:).*(deepRS_V(1,:)-p.deepRS_iNaP_ENa);
  deepRS_iNaP_minf_last=1./(1+exp(-(deepRS_V_last-(p.deepRS_iNaP_halfNaP+p.deepRS_iNaP_NaP_offset))/5));
deepRS_iNaP_minf = zeros(nsamp,p.deepRS_Npop);
  deepRS_iNaP_minf(1,:)=1./(1+exp(-(deepRS_V(1,:)-(p.deepRS_iNaP_halfNaP+p.deepRS_iNaP_NaP_offset))/5));
  deepRS_iPeriodicPulsesBen_Iext_last=p.deepRS_iPeriodicPulsesBen_PPstim*(( deepRS_iPeriodicPulsesBen_s2(k,:)));
deepRS_iPeriodicPulsesBen_Iext = zeros(nsamp,p.deepRS_Npop);
  deepRS_iPeriodicPulsesBen_Iext(1,:)=p.deepRS_iPeriodicPulsesBen_PPstim*(( deepRS_iPeriodicPulsesBen_s2(k,:)));
  deepRS_iPeriodicPulsesBen_input_last=deepRS_iPeriodicPulsesBen_s2(k,:);
deepRS_iPeriodicPulsesBen_input = zeros(nsamp,p.deepRS_Npop);
  deepRS_iPeriodicPulsesBen_input(1,:)=deepRS_iPeriodicPulsesBen_s2(k,:);
  deepRS_itonicBen_IBen_last=p.deepRS_itonicBen_I_app*((t/p.deepRS_itonicBen_ton)*(t<=p.deepRS_itonicBen_ton)+(p.deepRS_itonicBen_ton<t&t<deepRS_itonicBen_toff)+rand(1,p.deepRS_Npop)*p.deepRS_itonicBen_Inoise);
deepRS_itonicBen_IBen = zeros(nsamp,p.deepRS_Npop);
  deepRS_itonicBen_IBen(1,:)=p.deepRS_itonicBen_I_app*((t/p.deepRS_itonicBen_ton)*(t<=p.deepRS_itonicBen_ton)+(p.deepRS_itonicBen_ton<t&t<deepRS_itonicBen_toff)+rand(1,p.deepRS_Npop)*p.deepRS_itonicBen_Inoise);
% ###########################################################
% Numerical integration:
% ###########################################################
n=2;
for k=2:ntime
  t=T(k-1);
  deepFS_V_k1=(((-(( p.deepFS_IBitonic_stim*(t>p.deepFS_IBitonic_onset & t<p.deepFS_IBitonic_offset))))+(((( p.deepFS_IBnoise_V_noise.*randn(1,p.deepFS_Npop))))+((-(( p.deepFS_FSiNaF_gNaF.*(( 1./(1+exp((-deepFS_V_last-p.deepFS_FSiNaF_NaF_V0+p.deepFS_FSiNaF_NaF_offset)/10)))).^3.*deepFS_FSiNaF_hNaF_last.*(deepFS_V_last-p.deepFS_FSiNaF_E_NaF))))+((-(( p.deepFS_FSiKDR_gKDR.*deepFS_FSiKDR_mKDR_last.^4.*(deepFS_V_last-p.deepFS_FSiKDR_E_KDR))))+((-(( p.deepFS_IBleak_g_l.*(deepFS_V_last-p.deepFS_IBleak_E_l))))+((-(( (deepFS_deepRS_IBaIBdbiSYNseed_gsyn.*(deepFS_deepRS_IBaIBdbiSYNseed_s_last*deepFS_deepRS_IBaIBdbiSYNseed_mask).*(deepFS_V_last-p.deepFS_deepRS_IBaIBdbiSYNseed_E_SYN))))))))))))/p.deepFS_Cm;
  deepFS_FSiNaF_hNaF_k1= (( (( 1./(1+exp((deepFS_V_last+p.deepFS_FSiNaF_NaF_V1-p.deepFS_FSiNaF_NaF_offset)/p.deepFS_FSiNaF_NaF_d1)))) ./ (( p.deepFS_FSiNaF_NaF_c0 + p.deepFS_FSiNaF_NaF_c1./(1+exp((deepFS_V_last+p.deepFS_FSiNaF_NaF_V2-p.deepFS_FSiNaF_NaF_offset)/p.deepFS_FSiNaF_NaF_d2)))))).*(1-deepFS_FSiNaF_hNaF_last)-(( (1-(( 1./(1+exp((deepFS_V_last+p.deepFS_FSiNaF_NaF_V1-p.deepFS_FSiNaF_NaF_offset)/p.deepFS_FSiNaF_NaF_d1)))))./(( p.deepFS_FSiNaF_NaF_c0 + p.deepFS_FSiNaF_NaF_c1./(1+exp((deepFS_V_last+p.deepFS_FSiNaF_NaF_V2-p.deepFS_FSiNaF_NaF_offset)/p.deepFS_FSiNaF_NaF_d2)))))).*deepFS_FSiNaF_hNaF_last;
  deepFS_FSiKDR_mKDR_k1= (( (( 1./(1+exp((-deepFS_V_last-p.deepFS_FSiKDR_KDR_V1+p.deepFS_FSiKDR_KDR_offset)/p.deepFS_FSiKDR_KDR_d1)))) ./ (( .25+4.35*exp(-abs(deepFS_V_last+p.deepFS_FSiKDR_KDR_V2-p.deepFS_FSiKDR_KDR_offset)/p.deepFS_FSiKDR_KDR_d2))))).*(1-deepFS_FSiKDR_mKDR_last)-(( (1-(( 1./(1+exp((-deepFS_V_last-p.deepFS_FSiKDR_KDR_V1+p.deepFS_FSiKDR_KDR_offset)/p.deepFS_FSiKDR_KDR_d1)))))./(( .25+4.35*exp(-abs(deepFS_V_last+p.deepFS_FSiKDR_KDR_V2-p.deepFS_FSiKDR_KDR_offset)/p.deepFS_FSiKDR_KDR_d2))))).*deepFS_FSiKDR_mKDR_last;
  deepRS_V_k1=(p.deepRS_I_const+((-((p.deepRS_iNaP_gNaP.*deepRS_iNaP_m_last.*(deepRS_V_last-p.deepRS_iNaP_ENa))))+((-((p.deepRS_iKs_gKs.*deepRS_iKs_n_last.*(deepRS_V_last-p.deepRS_iKs_E_Ks))))+((-(( (p.deepRS_iKDRG_gKDR/p.deepRS_iKDRG_gdenomK).*deepRS_iKDRG_mKDR_last.^4.*(deepRS_V_last-p.deepRS_iKDRG_E_KDR))))+((-(( (p.deepRS_iNaG_gNa/p.deepRS_iNaG_gdenomN)*((( (( -(deepRS_V_last+16-p.deepRS_iNaG_Noffset)./(10*(exp(-(deepRS_V_last+16-p.deepRS_iNaG_Noffset)/10)-1))))./((( -(deepRS_V_last+16-p.deepRS_iNaG_Noffset)./(10*(exp(-(deepRS_V_last+16-p.deepRS_iNaG_Noffset)/10)-1))))+(( 4*exp(-(deepRS_V_last+41-p.deepRS_iNaG_Noffset)/18)))))).^3).*deepRS_iNaG_h_last.*(deepRS_V_last-p.deepRS_iNaG_ENa))))+((-(( (p.deepRS_iLeak_gl/p.deepRS_iLeak_gdenoml).*(deepRS_V_last-p.deepRS_iLeak_El))))+((-(( p.deepRS_iCaT_gCa.*(deepRS_iCaT_s_last.^2).*(deepRS_V_last-p.deepRS_iCaT_ECa))))+((-(( p.deepRS_iKCaT_gKCa*deepRS_iKCaT_q_last*(deepRS_V_last-p.deepRS_iKCaT_EK))))+((-(( p.deepRS_iPeriodicPulsesBen_PPstim*(( deepRS_iPeriodicPulsesBen_s2(k,:))))))+((-((p.deepRS_itonicBen_I_app*((t/p.deepRS_itonicBen_ton)*(t<=p.deepRS_itonicBen_ton)+(p.deepRS_itonicBen_ton<t&t<deepRS_itonicBen_toff)+rand(1,p.deepRS_Npop)*p.deepRS_itonicBen_Inoise))))+((-(( p.deepRS_iFMPulses_FMPstim*(( deepRS_iFMPulses_s2(k,:))))))+((-(( (deepRS_deepFS_IBaIBdbiSYNseed_gsyn.*(deepRS_deepFS_IBaIBdbiSYNseed_s_last*deepRS_deepFS_IBaIBdbiSYNseed_mask).*(deepRS_V_last-p.deepRS_deepFS_IBaIBdbiSYNseed_E_SYN)))))))))))))))))/p.deepRS_Cm;
  deepRS_iNaP_m_k1=(((1./(1+exp(-(deepRS_V_last-(p.deepRS_iNaP_halfNaP+p.deepRS_iNaP_NaP_offset))/5))))-deepRS_iNaP_m_last)/p.deepRS_iNaP_tauNaP;
  deepRS_iKs_n_k1=(((1./(1+exp(-(deepRS_V_last-(p.deepRS_iKs_halfKs+p.deepRS_iKs_Ks_offset))/10))))-deepRS_iKs_n_last)/((p.deepRS_iKs_tau_mult*1000./(3.3*deepRS_iKs_tau_div*(exp((deepRS_V_last-(p.deepRS_iKs_halfKs+p.deepRS_iKs_Ks_offset))/40)+exp(-(deepRS_V_last-(p.deepRS_iKs_halfKs+p.deepRS_iKs_Ks_offset))/20)))));
  deepRS_iKDRG_mKDR_k1= deepRS_iKDRG_tau_m*((( -.01*(deepRS_V_last+20-p.deepRS_iKDRG_Koffset)./(exp(-(deepRS_V_last+20-p.deepRS_iKDRG_Koffset)/10)-1))).*(1-deepRS_iKDRG_mKDR_last)-(( .125*exp(-(deepRS_V_last+30-p.deepRS_iKDRG_Koffset)/80))).*deepRS_iKDRG_mKDR_last);
  deepRS_iNaG_h_k1= deepRS_iNaG_tau_h*((( .07*exp(-(deepRS_V_last+30-p.deepRS_iNaG_Noffset)/20))).*(1-deepRS_iNaG_h_last)-(( 1./(exp(-(deepRS_V_last-p.deepRS_iNaG_Noffset)/10)+1))).*deepRS_iNaG_h_last);
  deepRS_CaDynT_cai_k1= - p.deepRS_CaDynT_CAF.*(((( p.deepRS_iCaT_gCa.*(deepRS_iCaT_s_last.^2).*(deepRS_V_last-p.deepRS_iCaT_ECa))))) - deepRS_CaDynT_cai_last/p.deepRS_CaDynT_tauCa;
  deepRS_iCaT_s_k1= ((( (( 1.6./(1+exp(-.072*(deepRS_V_last-p.deepRS_iCaT_Caoffset-65)))))./((( 1.6./(1+exp(-.072*(deepRS_V_last-p.deepRS_iCaT_Caoffset-65)))))+(( 0.02*(deepRS_V_last-p.deepRS_iCaT_Caoffset-51.1)./(exp((deepRS_V_last-p.deepRS_iCaT_Caoffset-51.1)/5)-1))))))-deepRS_iCaT_s_last)./(( 1./((( 1.6./(1+exp(-.072*(deepRS_V_last-p.deepRS_iCaT_Caoffset-65)))))+(( 0.02*(deepRS_V_last-p.deepRS_iCaT_Caoffset-51.1)./(exp((deepRS_V_last-p.deepRS_iCaT_Caoffset-51.1)/5)-1))))));
  deepRS_iKCaT_q_k1= ((( (( min(p.deepRS_iKCaT_aKCa_factor.*deepRS_CaDynT_cai_last,1)))./((( min(p.deepRS_iKCaT_aKCa_factor.*deepRS_CaDynT_cai_last,1)))+p.deepRS_iKCaT_bKCa)))-deepRS_iKCaT_q_last)./(( 1./((( min(p.deepRS_iKCaT_aKCa_factor.*deepRS_CaDynT_cai_last,1)))+p.deepRS_iKCaT_bKCa)));
  deepFS_deepRS_IBaIBdbiSYNseed_s_k1= -deepFS_deepRS_IBaIBdbiSYNseed_s_last./p.deepFS_deepRS_IBaIBdbiSYNseed_tauDx + ((1-deepFS_deepRS_IBaIBdbiSYNseed_s_last)/p.deepFS_deepRS_IBaIBdbiSYNseed_tauRx).*(1+tanh(deepRS_V_last/10));
  deepRS_deepFS_IBaIBdbiSYNseed_s_k1= -deepRS_deepFS_IBaIBdbiSYNseed_s_last./p.deepRS_deepFS_IBaIBdbiSYNseed_tauDx + ((1-deepRS_deepFS_IBaIBdbiSYNseed_s_last)/p.deepRS_deepFS_IBaIBdbiSYNseed_tauRx).*(1+tanh(deepFS_V_last/10));
  % ------------------------------------------------------------
  % Update state variables:
  % ------------------------------------------------------------
  deepFS_V_last=deepFS_V_last+dt*deepFS_V_k1;
  deepFS_FSiNaF_hNaF_last=deepFS_FSiNaF_hNaF_last+dt*deepFS_FSiNaF_hNaF_k1;
  deepFS_FSiKDR_mKDR_last=deepFS_FSiKDR_mKDR_last+dt*deepFS_FSiKDR_mKDR_k1;
  deepRS_V_last=deepRS_V_last+dt*deepRS_V_k1;
  deepRS_iNaP_m_last=deepRS_iNaP_m_last+dt*deepRS_iNaP_m_k1;
  deepRS_iKs_n_last=deepRS_iKs_n_last+dt*deepRS_iKs_n_k1;
  deepRS_iKDRG_mKDR_last=deepRS_iKDRG_mKDR_last+dt*deepRS_iKDRG_mKDR_k1;
  deepRS_iNaG_h_last=deepRS_iNaG_h_last+dt*deepRS_iNaG_h_k1;
  deepRS_CaDynT_cai_last=deepRS_CaDynT_cai_last+dt*deepRS_CaDynT_cai_k1;
  deepRS_iCaT_s_last=deepRS_iCaT_s_last+dt*deepRS_iCaT_s_k1;
  deepRS_iKCaT_q_last=deepRS_iKCaT_q_last+dt*deepRS_iKCaT_q_k1;
  deepFS_deepRS_IBaIBdbiSYNseed_s_last=deepFS_deepRS_IBaIBdbiSYNseed_s_last+dt*deepFS_deepRS_IBaIBdbiSYNseed_s_k1;
  deepRS_deepFS_IBaIBdbiSYNseed_s_last=deepRS_deepFS_IBaIBdbiSYNseed_s_last+dt*deepRS_deepFS_IBaIBdbiSYNseed_s_k1;
if mod(k,downsample_factor)==0 % store this time point
  % ------------------------------------------------------------
  % Store state variables:
  % ------------------------------------------------------------
  deepFS_V(n)=deepFS_V_last;
  deepFS_FSiNaF_hNaF(n)=deepFS_FSiNaF_hNaF_last;
  deepFS_FSiKDR_mKDR(n)=deepFS_FSiKDR_mKDR_last;
  deepRS_V(n)=deepRS_V_last;
  deepRS_iNaP_m(n)=deepRS_iNaP_m_last;
  deepRS_iKs_n(n)=deepRS_iKs_n_last;
  deepRS_iKDRG_mKDR(n)=deepRS_iKDRG_mKDR_last;
  deepRS_iNaG_h(n)=deepRS_iNaG_h_last;
  deepRS_CaDynT_cai(n)=deepRS_CaDynT_cai_last;
  deepRS_iCaT_s(n)=deepRS_iCaT_s_last;
  deepRS_iKCaT_q(n)=deepRS_iKCaT_q_last;
  deepFS_deepRS_IBaIBdbiSYNseed_s(n)=deepFS_deepRS_IBaIBdbiSYNseed_s_last;
  deepRS_deepFS_IBaIBdbiSYNseed_s(n)=deepRS_deepFS_IBaIBdbiSYNseed_s_last;
  % ------------------------------------------------------------
  % Update monitors:
  % ------------------------------------------------------------
  deepFS_FSiKDR_IKDR(n,:)=p.deepFS_FSiKDR_gKDR.*deepFS_FSiKDR_mKDR(n,:).^4.*(deepFS_V(n,:)-p.deepFS_FSiKDR_E_KDR);
  deepFS_FSiKDR_aM(n,:)=(( 1./(1+exp((-deepFS_V(n,:)-p.deepFS_FSiKDR_KDR_V1+p.deepFS_FSiKDR_KDR_offset)/p.deepFS_FSiKDR_KDR_d1)))) ./ (( .25+4.35*exp(-abs(deepFS_V(n,:)+p.deepFS_FSiKDR_KDR_V2-p.deepFS_FSiKDR_KDR_offset)/p.deepFS_FSiKDR_KDR_d2)));
  deepFS_FSiKDR_bM(n,:)=(1-(( 1./(1+exp((-deepFS_V(n,:)-p.deepFS_FSiKDR_KDR_V1+p.deepFS_FSiKDR_KDR_offset)/p.deepFS_FSiKDR_KDR_d1)))))./(( .25+4.35*exp(-abs(deepFS_V(n,:)+p.deepFS_FSiKDR_KDR_V2-p.deepFS_FSiKDR_KDR_offset)/p.deepFS_FSiKDR_KDR_d2)));
  deepFS_FSiKDR_minf(n,:)=1./(1+exp((-deepFS_V(n,:)-p.deepFS_FSiKDR_KDR_V1+p.deepFS_FSiKDR_KDR_offset)/p.deepFS_FSiKDR_KDR_d1));
  deepFS_FSiKDR_mtau(n,:)=.25+4.35*exp(-abs(deepFS_V(n,:)+p.deepFS_FSiKDR_KDR_V2-p.deepFS_FSiKDR_KDR_offset)/p.deepFS_FSiKDR_KDR_d2);
  deepFS_FSiNaF_INaF(n,:)=p.deepFS_FSiNaF_gNaF.*(( 1./(1+exp((-deepFS_V(n,:)-p.deepFS_FSiNaF_NaF_V0+p.deepFS_FSiNaF_NaF_offset)/10)))).^3.*deepFS_FSiNaF_hNaF(n,:).*(deepFS_V(n,:)-p.deepFS_FSiNaF_E_NaF);
  deepFS_FSiNaF_aH(n,:)=(( 1./(1+exp((deepFS_V(n,:)+p.deepFS_FSiNaF_NaF_V1-p.deepFS_FSiNaF_NaF_offset)/p.deepFS_FSiNaF_NaF_d1)))) ./ (( p.deepFS_FSiNaF_NaF_c0 + p.deepFS_FSiNaF_NaF_c1./(1+exp((deepFS_V(n,:)+p.deepFS_FSiNaF_NaF_V2-p.deepFS_FSiNaF_NaF_offset)/p.deepFS_FSiNaF_NaF_d2))));
  deepFS_FSiNaF_bH(n,:)=(1-(( 1./(1+exp((deepFS_V(n,:)+p.deepFS_FSiNaF_NaF_V1-p.deepFS_FSiNaF_NaF_offset)/p.deepFS_FSiNaF_NaF_d1)))))./(( p.deepFS_FSiNaF_NaF_c0 + p.deepFS_FSiNaF_NaF_c1./(1+exp((deepFS_V(n,:)+p.deepFS_FSiNaF_NaF_V2-p.deepFS_FSiNaF_NaF_offset)/p.deepFS_FSiNaF_NaF_d2))));
  deepFS_FSiNaF_hinf(n,:)=1./(1+exp((deepFS_V(n,:)+p.deepFS_FSiNaF_NaF_V1-p.deepFS_FSiNaF_NaF_offset)/p.deepFS_FSiNaF_NaF_d1));
  deepFS_FSiNaF_htau(n,:)=p.deepFS_FSiNaF_NaF_c0 + p.deepFS_FSiNaF_NaF_c1./(1+exp((deepFS_V(n,:)+p.deepFS_FSiNaF_NaF_V2-p.deepFS_FSiNaF_NaF_offset)/p.deepFS_FSiNaF_NaF_d2));
  deepFS_FSiNaF_m0(n,:)=1./(1+exp((-deepFS_V(n,:)-p.deepFS_FSiNaF_NaF_V0+p.deepFS_FSiNaF_NaF_offset)/10));
  deepFS_IBitonic_Itonic(n,:)=p.deepFS_IBitonic_stim*(t>p.deepFS_IBitonic_onset & t<p.deepFS_IBitonic_offset);
  deepFS_IBleak_Ileak(n,:)=p.deepFS_IBleak_g_l.*(deepFS_V(n,:)-p.deepFS_IBleak_E_l);
  deepFS_IBnoise_Inoise(n,:)=p.deepFS_IBnoise_V_noise.*randn(1,p.deepFS_Npop);
  deepFS_deepRS_IBaIBdbiSYNseed_ISYN(n,:)=(deepFS_deepRS_IBaIBdbiSYNseed_gsyn.*(deepFS_deepRS_IBaIBdbiSYNseed_s(n,:)*deepFS_deepRS_IBaIBdbiSYNseed_mask).*(deepFS_V(n,:)-p.deepFS_deepRS_IBaIBdbiSYNseed_E_SYN));
  deepRS_deepFS_IBaIBdbiSYNseed_ISYN(n,:)=(deepRS_deepFS_IBaIBdbiSYNseed_gsyn.*(deepRS_deepFS_IBaIBdbiSYNseed_s(n,:)*deepRS_deepFS_IBaIBdbiSYNseed_mask).*(deepRS_V(n,:)-p.deepRS_deepFS_IBaIBdbiSYNseed_E_SYN));
  deepRS_iCaT_ICa(n,:)=p.deepRS_iCaT_gCa.*(deepRS_iCaT_s(n,:).^2).*(deepRS_V(n,:)-p.deepRS_iCaT_ECa);
  deepRS_iCaT_as(n,:)=1.6./(1+exp(-.072*(deepRS_V(n,:)-p.deepRS_iCaT_Caoffset-65)));
  deepRS_iCaT_bs(n,:)=0.02*(deepRS_V(n,:)-p.deepRS_iCaT_Caoffset-51.1)./(exp((deepRS_V(n,:)-p.deepRS_iCaT_Caoffset-51.1)/5)-1);
  deepRS_iCaT_sinf(n,:)=(( 1.6./(1+exp(-.072*(deepRS_V(n,:)-p.deepRS_iCaT_Caoffset-65)))))./((( 1.6./(1+exp(-.072*(deepRS_V(n,:)-p.deepRS_iCaT_Caoffset-65)))))+(( 0.02*(deepRS_V(n,:)-p.deepRS_iCaT_Caoffset-51.1)./(exp((deepRS_V(n,:)-p.deepRS_iCaT_Caoffset-51.1)/5)-1))));
  deepRS_iCaT_stau(n,:)=1./((( 1.6./(1+exp(-.072*(deepRS_V(n,:)-p.deepRS_iCaT_Caoffset-65)))))+(( 0.02*(deepRS_V(n,:)-p.deepRS_iCaT_Caoffset-51.1)./(exp((deepRS_V(n,:)-p.deepRS_iCaT_Caoffset-51.1)/5)-1))));
  deepRS_iFMPulses_Iext(n,:)=p.deepRS_iFMPulses_FMPstim*(( deepRS_iFMPulses_s2(k,:)));
  deepRS_iFMPulses_input(n,:)=deepRS_iFMPulses_s2(k,:);
  deepRS_iKCaT_IKCa(n,:)=p.deepRS_iKCaT_gKCa*deepRS_iKCaT_q(n,:)*(deepRS_V(n,:)-p.deepRS_iKCaT_EK);
  deepRS_iKCaT_aKCa(n,:)=min(p.deepRS_iKCaT_aKCa_factor.*deepRS_CaDynT_cai(n,:),1);
  deepRS_iKCaT_infKCa(n,:)=(( min(p.deepRS_iKCaT_aKCa_factor.*deepRS_CaDynT_cai(n,:),1)))./((( min(p.deepRS_iKCaT_aKCa_factor.*deepRS_CaDynT_cai(n,:),1)))+p.deepRS_iKCaT_bKCa);
  deepRS_iKCaT_tauKCa(n,:)=1./((( min(p.deepRS_iKCaT_aKCa_factor.*deepRS_CaDynT_cai(n,:),1)))+p.deepRS_iKCaT_bKCa);
  deepRS_iKDRG_IKDR(n,:)=(p.deepRS_iKDRG_gKDR/p.deepRS_iKDRG_gdenomK).*deepRS_iKDRG_mKDR(n,:).^4.*(deepRS_V(n,:)-p.deepRS_iKDRG_E_KDR);
  deepRS_iKDRG_aM(n,:)=-.01*(deepRS_V(n,:)+20-p.deepRS_iKDRG_Koffset)./(exp(-(deepRS_V(n,:)+20-p.deepRS_iKDRG_Koffset)/10)-1);
  deepRS_iKDRG_bM(n,:)=.125*exp(-(deepRS_V(n,:)+30-p.deepRS_iKDRG_Koffset)/80);
  deepRS_iKs_IKs(n,:)=p.deepRS_iKs_gKs.*deepRS_iKs_n(n,:).*(deepRS_V(n,:)-p.deepRS_iKs_E_Ks);
  deepRS_iKs_ninf(n,:)=1./(1+exp(-(deepRS_V(n,:)-(p.deepRS_iKs_halfKs+p.deepRS_iKs_Ks_offset))/10));
  deepRS_iKs_tauKs(n,:)=p.deepRS_iKs_tau_mult*1000./(3.3*deepRS_iKs_tau_div*(exp((deepRS_V(n,:)-(p.deepRS_iKs_halfKs+p.deepRS_iKs_Ks_offset))/40)+exp(-(deepRS_V(n,:)-(p.deepRS_iKs_halfKs+p.deepRS_iKs_Ks_offset))/20)));
  deepRS_iLeak_Il(n,:)=(p.deepRS_iLeak_gl/p.deepRS_iLeak_gdenoml).*(deepRS_V(n,:)-p.deepRS_iLeak_El);
  deepRS_iNaG_INa(n,:)=(p.deepRS_iNaG_gNa/p.deepRS_iNaG_gdenomN)*((( (( -(deepRS_V(n,:)+16-p.deepRS_iNaG_Noffset)./(10*(exp(-(deepRS_V(n,:)+16-p.deepRS_iNaG_Noffset)/10)-1))))./((( -(deepRS_V(n,:)+16-p.deepRS_iNaG_Noffset)./(10*(exp(-(deepRS_V(n,:)+16-p.deepRS_iNaG_Noffset)/10)-1))))+(( 4*exp(-(deepRS_V(n,:)+41-p.deepRS_iNaG_Noffset)/18)))))).^3).*deepRS_iNaG_h(n,:).*(deepRS_V(n,:)-p.deepRS_iNaG_ENa);
  deepRS_iNaG_ah(n,:)=.07*exp(-(deepRS_V(n,:)+30-p.deepRS_iNaG_Noffset)/20);
  deepRS_iNaG_am(n,:)=-(deepRS_V(n,:)+16-p.deepRS_iNaG_Noffset)./(10*(exp(-(deepRS_V(n,:)+16-p.deepRS_iNaG_Noffset)/10)-1));
  deepRS_iNaG_bh(n,:)=1./(exp(-(deepRS_V(n,:)-p.deepRS_iNaG_Noffset)/10)+1);
  deepRS_iNaG_bm(n,:)=4*exp(-(deepRS_V(n,:)+41-p.deepRS_iNaG_Noffset)/18);
  deepRS_iNaG_minf(n,:)=(( -(deepRS_V(n,:)+16-p.deepRS_iNaG_Noffset)./(10*(exp(-(deepRS_V(n,:)+16-p.deepRS_iNaG_Noffset)/10)-1))))./((( -(deepRS_V(n,:)+16-p.deepRS_iNaG_Noffset)./(10*(exp(-(deepRS_V(n,:)+16-p.deepRS_iNaG_Noffset)/10)-1))))+(( 4*exp(-(deepRS_V(n,:)+41-p.deepRS_iNaG_Noffset)/18))));
  deepRS_iNaP_INaP(n,:)=p.deepRS_iNaP_gNaP.*deepRS_iNaP_m(n,:).*(deepRS_V(n,:)-p.deepRS_iNaP_ENa);
  deepRS_iNaP_minf(n,:)=1./(1+exp(-(deepRS_V(n,:)-(p.deepRS_iNaP_halfNaP+p.deepRS_iNaP_NaP_offset))/5));
  deepRS_iPeriodicPulsesBen_Iext(n,:)=p.deepRS_iPeriodicPulsesBen_PPstim*(( deepRS_iPeriodicPulsesBen_s2(k,:)));
  deepRS_iPeriodicPulsesBen_input(n,:)=deepRS_iPeriodicPulsesBen_s2(k,:);
  deepRS_itonicBen_IBen(n,:)=p.deepRS_itonicBen_I_app*((t/p.deepRS_itonicBen_ton)*(t<=p.deepRS_itonicBen_ton)+(p.deepRS_itonicBen_ton<t&t<deepRS_itonicBen_toff)+rand(1,p.deepRS_Npop)*p.deepRS_itonicBen_Inoise);
  n=n+1;
end
end
T=T(1:downsample_factor:ntime);

end

