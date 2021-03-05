function [T,pop1_v,pop1_iNa_m,pop1_iNa_h,pop1_iK_n]=solve_ode
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
% Initial conditions:
% ------------------------------------------------------------
% seed the random number generator
rng(p.random_seed);
t=0; k=1;

% STATE_VARIABLES:
pop1_v = zeros(nsamp,p.pop1_Npop);
  pop1_v(1,:) = zeros(1,p.pop1_Npop);
pop1_iNa_m = zeros(nsamp,p.pop1_Npop);
  pop1_iNa_m(1,:) = p.pop1_iNa_m_IC+p.pop1_iNa_IC_noise*rand(1,p.pop1_Npop);
pop1_iNa_h = zeros(nsamp,p.pop1_Npop);
  pop1_iNa_h(1,:) = p.pop1_iNa_h_IC+p.pop1_iNa_IC_noise*rand(1,p.pop1_Npop);
pop1_iK_n = zeros(nsamp,p.pop1_Npop);
  pop1_iK_n(1,:) = p.pop1_iK_n_IC+p.pop1_iK_IC_noise*rand(1,p.pop1_Npop);


% ###########################################################
% Numerical integration:
% ###########################################################
% seed the random number generator
rng(p.random_seed);
n=2;
for k=2:ntime
  t=T(k-1);
  pop1_v_k1 =(((( -p.pop1_iNa_gNa.*pop1_iNa_m(n-1).^3.*pop1_iNa_h(n-1).*(pop1_v(n-1)-p.pop1_iNa_ENa))))+(((( -p.pop1_iK_gK.*pop1_iK_n(n-1).^4.*(pop1_v(n-1)-p.pop1_iK_EK))))))+p.pop1_I;
  pop1_iNa_m_k1 = (( (2.5-.1*(pop1_v(n-1)+65))./(exp(2.5-.1*(pop1_v(n-1)+65))-1))).*(1-pop1_iNa_m(n-1))-(( 4*exp(-(pop1_v(n-1)+65)/18))).*pop1_iNa_m(n-1);
  pop1_iNa_h_k1 = (( .07*exp(-(pop1_v(n-1)+65)/20))).*(1-pop1_iNa_h(n-1))-(( 1./(exp(3-.1*(pop1_v(n-1)+65))+1))).*pop1_iNa_h(n-1);
  pop1_iK_n_k1 = (( (.1-.01*(pop1_v(n-1)+65))./(exp(1-.1*(pop1_v(n-1)+65))-1))).*(1-pop1_iK_n(n-1))-(( .125*exp(-(pop1_v(n-1)+65)/80))).*pop1_iK_n(n-1);

  t = t + .5*dt;
  pop1_v_k2 =(((( -p.pop1_iNa_gNa.*((pop1_iNa_m(n-1) + .5*dt*pop1_iNa_m_k1)).^3.*((pop1_iNa_h(n-1) + .5*dt*pop1_iNa_h_k1)).*(((pop1_v(n-1) + .5*dt*pop1_v_k1))-p.pop1_iNa_ENa))))+(((( -p.pop1_iK_gK.*((pop1_iK_n(n-1) + .5*dt*pop1_iK_n_k1)).^4.*(((pop1_v(n-1) + .5*dt*pop1_v_k1))-p.pop1_iK_EK))))))+p.pop1_I;
  pop1_iNa_m_k2 = (( (2.5-.1*(((pop1_v(n-1) + .5*dt*pop1_v_k1))+65))./(exp(2.5-.1*(((pop1_v(n-1) + .5*dt*pop1_v_k1))+65))-1))).*(1-((pop1_iNa_m(n-1) + .5*dt*pop1_iNa_m_k1)))-(( 4*exp(-(((pop1_v(n-1) + .5*dt*pop1_v_k1))+65)/18))).*((pop1_iNa_m(n-1) + .5*dt*pop1_iNa_m_k1));
  pop1_iNa_h_k2 = (( .07*exp(-(((pop1_v(n-1) + .5*dt*pop1_v_k1))+65)/20))).*(1-((pop1_iNa_h(n-1) + .5*dt*pop1_iNa_h_k1)))-(( 1./(exp(3-.1*(((pop1_v(n-1) + .5*dt*pop1_v_k1))+65))+1))).*((pop1_iNa_h(n-1) + .5*dt*pop1_iNa_h_k1));
  pop1_iK_n_k2 = (( (.1-.01*(((pop1_v(n-1) + .5*dt*pop1_v_k1))+65))./(exp(1-.1*(((pop1_v(n-1) + .5*dt*pop1_v_k1))+65))-1))).*(1-((pop1_iK_n(n-1) + .5*dt*pop1_iK_n_k1)))-(( .125*exp(-(((pop1_v(n-1) + .5*dt*pop1_v_k1))+65)/80))).*((pop1_iK_n(n-1) + .5*dt*pop1_iK_n_k1));

  pop1_v_k3 =(((( -p.pop1_iNa_gNa.*((pop1_iNa_m(n-1) + .5*dt*pop1_iNa_m_k2)).^3.*((pop1_iNa_h(n-1) + .5*dt*pop1_iNa_h_k2)).*(((pop1_v(n-1) + .5*dt*pop1_v_k2))-p.pop1_iNa_ENa))))+(((( -p.pop1_iK_gK.*((pop1_iK_n(n-1) + .5*dt*pop1_iK_n_k2)).^4.*(((pop1_v(n-1) + .5*dt*pop1_v_k2))-p.pop1_iK_EK))))))+p.pop1_I;
  pop1_iNa_m_k3 = (( (2.5-.1*(((pop1_v(n-1) + .5*dt*pop1_v_k2))+65))./(exp(2.5-.1*(((pop1_v(n-1) + .5*dt*pop1_v_k2))+65))-1))).*(1-((pop1_iNa_m(n-1) + .5*dt*pop1_iNa_m_k2)))-(( 4*exp(-(((pop1_v(n-1) + .5*dt*pop1_v_k2))+65)/18))).*((pop1_iNa_m(n-1) + .5*dt*pop1_iNa_m_k2));
  pop1_iNa_h_k3 = (( .07*exp(-(((pop1_v(n-1) + .5*dt*pop1_v_k2))+65)/20))).*(1-((pop1_iNa_h(n-1) + .5*dt*pop1_iNa_h_k2)))-(( 1./(exp(3-.1*(((pop1_v(n-1) + .5*dt*pop1_v_k2))+65))+1))).*((pop1_iNa_h(n-1) + .5*dt*pop1_iNa_h_k2));
  pop1_iK_n_k3 = (( (.1-.01*(((pop1_v(n-1) + .5*dt*pop1_v_k2))+65))./(exp(1-.1*(((pop1_v(n-1) + .5*dt*pop1_v_k2))+65))-1))).*(1-((pop1_iK_n(n-1) + .5*dt*pop1_iK_n_k2)))-(( .125*exp(-(((pop1_v(n-1) + .5*dt*pop1_v_k2))+65)/80))).*((pop1_iK_n(n-1) + .5*dt*pop1_iK_n_k2));

  t = t + .5*dt;
  pop1_v_k4 =(((( -p.pop1_iNa_gNa.*((pop1_iNa_m(n-1) + dt*pop1_iNa_m_k3)).^3.*((pop1_iNa_h(n-1) + dt*pop1_iNa_h_k3)).*(((pop1_v(n-1) + dt*pop1_v_k3))-p.pop1_iNa_ENa))))+(((( -p.pop1_iK_gK.*((pop1_iK_n(n-1) + dt*pop1_iK_n_k3)).^4.*(((pop1_v(n-1) + dt*pop1_v_k3))-p.pop1_iK_EK))))))+p.pop1_I;
  pop1_iNa_m_k4 = (( (2.5-.1*(((pop1_v(n-1) + dt*pop1_v_k3))+65))./(exp(2.5-.1*(((pop1_v(n-1) + dt*pop1_v_k3))+65))-1))).*(1-((pop1_iNa_m(n-1) + dt*pop1_iNa_m_k3)))-(( 4*exp(-(((pop1_v(n-1) + dt*pop1_v_k3))+65)/18))).*((pop1_iNa_m(n-1) + dt*pop1_iNa_m_k3));
  pop1_iNa_h_k4 = (( .07*exp(-(((pop1_v(n-1) + dt*pop1_v_k3))+65)/20))).*(1-((pop1_iNa_h(n-1) + dt*pop1_iNa_h_k3)))-(( 1./(exp(3-.1*(((pop1_v(n-1) + dt*pop1_v_k3))+65))+1))).*((pop1_iNa_h(n-1) + dt*pop1_iNa_h_k3));
  pop1_iK_n_k4 = (( (.1-.01*(((pop1_v(n-1) + dt*pop1_v_k3))+65))./(exp(1-.1*(((pop1_v(n-1) + dt*pop1_v_k3))+65))-1))).*(1-((pop1_iK_n(n-1) + dt*pop1_iK_n_k3)))-(( .125*exp(-(((pop1_v(n-1) + dt*pop1_v_k3))+65)/80))).*((pop1_iK_n(n-1) + dt*pop1_iK_n_k3));

  % ------------------------------------------------------------
  % Update state variables:
  % ------------------------------------------------------------
  pop1_v(n) = pop1_v(n-1)+(dt/6)*(pop1_v_k1 + 2*(pop1_v_k2 + pop1_v_k3) + pop1_v_k4);
  pop1_iNa_m(n) = pop1_iNa_m(n-1)+(dt/6)*(pop1_iNa_m_k1 + 2*(pop1_iNa_m_k2 + pop1_iNa_m_k3) + pop1_iNa_m_k4);
  pop1_iNa_h(n) = pop1_iNa_h(n-1)+(dt/6)*(pop1_iNa_h_k1 + 2*(pop1_iNa_h_k2 + pop1_iNa_h_k3) + pop1_iNa_h_k4);
  pop1_iK_n(n) = pop1_iK_n(n-1)+(dt/6)*(pop1_iK_n_k1 + 2*(pop1_iK_n_k2 + pop1_iK_n_k3) + pop1_iK_n_k4);
  n=n+1;
end

T=T(1:downsample_factor:ntime);

end

