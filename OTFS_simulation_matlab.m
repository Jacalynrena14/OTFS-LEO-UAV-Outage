%% =========================================================================
%  OTFS Fig.1 — PURELY ANALYTICAL — FINAL CORRECT VERSION
%  Paper: IEEE TVT Vol.71, No.3, March 2022
%
%  CORRECT METHOD:
%  SNR = Ps/sigma2 * (Lsd/Zsd + xi*Lrd/Zrd)   [Eq.10, Eq.11]
%  Outage: I_otfs < Rb  =>  Lsd/Zsd + xi*Lrd/Zrd < eta
%  So the CDF variable is W = Lsd/Zsd + xi*Lrd/Zrd
%
%  Since Zsd ~ InvGamma  =>  1/Zsd ~ Gamma
%  W = sum of scaled Gamma vars => W ~ Gamma (moment matching)
%  P_out = gammainc(eta * beta_W, alpha_W)   [Theorem 1, Eq.20]
%
%  STEPS:
%  1. SR stats: mu_s, var_s via Eq.(13)(14) using hypergeom
%  2. chi_sd -> Gamma(alpha_sd, beta_sd) moment matching
%  3. Zsd ~ IG -> moments mu_Zsd, var_Zsd via Eq.(16)
%  4. 1/Zsd ~ Gamma -> E[1/Zsd], Var[1/Zsd] via IG property
%  5. Zrd ~ IG -> moments via Eq.(18) -> E[1/Zrd], Var[1/Zrd]
%  6. W = Lsd/Zsd + xi*Lrd/Zrd -> Gamma(alpha_W, beta_W) moment match
%  7. P_out = gammainc(eta*beta_W, alpha_W)    eta=sigma2*(2^{2Rb}-1)/Ps
% =========================================================================

clc; clear; close all;

%% PARAMETERS
N=4; M=4; NM=N*M; Ns=10;
b_sd=0.251; m_sd=0.251; Omega_sd=0.279; Omega_rd=1;
sigma2=1; Rb=1;
Ps_dB=0:0.5:25; Ps_lin=10.^(Ps_dB/10);
L_sd=1; L_rd=1;  % normalized path loss

%% STEP 1: SR stats Eq.(13)(14)
beta_sr  = 1/(2*b_sd);
delta_sd = beta_sr*Omega_sd/(2*b_sd*m_sd+Omega_sd);
a_sd     = sqrt(Omega_sd/(2*b_sd*m_sd+Omega_sd));
z        = delta_sd/beta_sr;

F1 = hypergeom([3/2, m_sd], 1, z);
F2 = hypergeom([2,   m_sd], 1, z);

mu_s  = a_sd*sqrt(pi)/(2*sqrt(beta_sr))*F1;
var_s = (a_sd/delta_sd)*F2 - (pi*a_sd^2/(4*delta_sd))*F1^2;
fprintf('mu_s=%.6f  var_s=%.6f\n',mu_s,var_s);

%% STEP 2: chi_sd -> Gamma(alpha_sd, beta_sd)
alpha_sd = Ns*mu_s^2/var_s;
beta_sd  = mu_s/var_s;
fprintf('alpha_sd=%.4f (must>4)  beta_sd=%.4f\n',alpha_sd,beta_sd);
if alpha_sd<=4, error('Increase Ns! alpha_sd=%.2f',alpha_sd); end

%% STEP 3: Zsd ~ InvGamma moments from Eq.(16)
a1=alpha_sd-1; a2=alpha_sd-2; a3=alpha_sd-3; a4=alpha_sd-4; a5=2*alpha_sd-5;
mu_Zsd  = NM*a1*a2*a3*a4 / (beta_sd^2*(NM*a3*a4-2*a5));
var_Zsd = 2*NM^2*a5*(a1*a2*a3*a4)^2 / ...
          (beta_sd^4*(NM*a3*a4-2*a5)^2*(NM*a3*a4-4*a5));
fprintf('mu_Zsd=%.6f  var_Zsd=%.4e\n',mu_Zsd,var_Zsd);

%% STEP 4: 1/Zsd ~ Gamma via InvGamma property
% Zsd ~ IG(alpha_ig, beta_ig):  E[Z]=beta_ig/(alpha_ig-1)
% => alpha_ig = mu_Zsd^2/var_Zsd + 2
% => beta_ig  = mu_Zsd*(alpha_ig - 1)
% 1/Zsd ~ Gamma(alpha_ig, scale=1/beta_ig):
% E[1/Zsd]   = alpha_ig / beta_ig
% Var[1/Zsd] = alpha_ig / beta_ig^2
alpha_ig_sd = mu_Zsd^2/var_Zsd + 2;
beta_ig_sd  = mu_Zsd*(alpha_ig_sd - 1);
E_invZsd    = alpha_ig_sd / beta_ig_sd;
V_invZsd    = alpha_ig_sd / beta_ig_sd^2;
fprintf('E[1/Zsd]=%.6f  Var[1/Zsd]=%.4e\n',E_invZsd,V_invZsd);

%% STEP 5: Zrd ~ IG moments Eq.(18) -> 1/Zrd moments
% (mrd, xi) pairs from paper Fig.1
cases=[2,0.3; 2,0.5; 3,0.5; 4,0.5];
E_invZrd=zeros(4,1); V_invZrd=zeros(4,1);

for c=1:4
    mrd=cases(c,1); br=mrd/Omega_rd;
    if mrd>2
        b1=mrd-1; b2=mrd-2;
        mu_Zrd  = NM*b1*b2/(br*(NM*b2-1));
        var_Zrd = NM^2*b1^2*b2^2/(br^2*(NM*b2-1)^2*(NM*b2-2));
    else
        % mrd=2: Eq.(18) invalid, sample Zrd stats numerically
        S=2e6;
        Ds=gamrnd(2,Omega_rd/2,[NM,S]);
        Zs=mean(Ds.^(-1),1);
        mu_Zrd=mean(Zs); var_Zrd=var(Zs);
    end
    % Convert Zrd ~ IG -> 1/Zrd moments
    alpha_ig_rd = mu_Zrd^2/var_Zrd + 2;
    beta_ig_rd  = mu_Zrd*(alpha_ig_rd - 1);
    E_invZrd(c) = alpha_ig_rd / beta_ig_rd;
    V_invZrd(c) = alpha_ig_rd / beta_ig_rd^2;
    fprintf('mrd=%d xi=%.1f: E[1/Zrd]=%.4f Var[1/Zrd]=%.4e\n',...
        mrd,cases(c,2),E_invZrd(c),V_invZrd(c));
end

%% STEP 6+7: P_out = gammainc(eta*beta_W, alpha_W)  Theorem 1 Eq.(20)
% W = Lsd/Zsd + xi*Lrd/Zrd
% mu_W  = Lsd*E[1/Zsd] + xi*Lrd*E[1/Zrd]
% var_W = Lsd^2*Var[1/Zsd] + (xi*Lrd)^2*Var[1/Zrd]
% alpha_W = mu_W^2/var_W,  beta_W = mu_W/var_W  (rate param)
% P_out = gammainc(eta*beta_W, alpha_W)

% No relay: W = Lsd/Zsd only (xi=0)
P_nr=zeros(1,length(Ps_dB));
for s=1:length(Ps_dB)
    eta  = sigma2*(2^(2*Rb)-1)/Ps_lin(s);
    mu_W = L_sd*E_invZsd;
    vr_W = L_sd^2*V_invZsd;
    P_nr(s) = gammainc(eta*(mu_W/vr_W), mu_W^2/vr_W, 'lower');
end

% UAV cooperation
P_UAV=zeros(4,length(Ps_dB));
for c=1:4
    xi=cases(c,2);
    for s=1:length(Ps_dB)
        eta  = sigma2*(2^(2*Rb)-1)/Ps_lin(s);
        mu_W = L_sd*E_invZsd  + xi*L_rd*E_invZrd(c);
        vr_W = L_sd^2*V_invZsd + (xi*L_rd)^2*V_invZrd(c);
        P_UAV(c,s) = gammainc(eta*(mu_W/vr_W), mu_W^2/vr_W, 'lower');
    end
end

%% PLOT
figure('Name','Fig1 OTFS FINAL','Position',[100,100,700,560]);
hold on; grid on; box on;
set(gca,'YScale','log','FontSize',12);
mi=1:6:length(Ps_dB);

semilogy(Ps_dB,P_nr,      '--k', 'DisplayName','without relay',     'LineWidth',2);
semilogy(Ps_dB,P_UAV(1,:),'-*b', 'DisplayName','m_{rd}=2, \xi=0.3','LineWidth',1.6,'MarkerSize',8,'MarkerIndices',mi);
semilogy(Ps_dB,P_UAV(2,:),'-ob', 'DisplayName','m_{rd}=2, \xi=0.5','LineWidth',1.6,'MarkerSize',7,'MarkerIndices',mi);
semilogy(Ps_dB,P_UAV(3,:),'-vr', 'DisplayName','m_{rd}=3, \xi=0.5','LineWidth',1.6,'MarkerSize',7,'MarkerIndices',mi);
semilogy(Ps_dB,P_UAV(4,:),'-hg', 'DisplayName','m_{rd}=4, \xi=0.5','LineWidth',1.6,'MarkerSize',8,'MarkerIndices',mi);

xlabel('the transmit SNR (dB)','FontSize',13);
ylabel('outage probability','FontSize',13);
title('Fig.1  UAV cooperation — OTFS LEO-Sat (N=M=4)','FontSize',12);
legend('Location','southwest','FontSize',10);
ylim([1e-4 1]); xlim([0 25]);
fprintf('DONE!\n');