clear all
clc
close all

%% Set precision and define parameters
digits(2); % Limit variable display precision to 2 digits
n = 5;
downsample_factor = 100; % Downsampling factor for saving output data
p = parse_inputs(); % Load model parameters

%% Load maximum values for boundary conditions
load('WT_max.mat')
load('WSS_max.mat');
WSS_max = WSS(:,2); % Extract max wall shear stress values

%% Build simplified cerebrovascular network for hemodynamic simulation
% The model includes 4 penetrating arterioles (PAs) to reduce computational load
% One of the PAs is selected to simulate cellular-level hemodynamics
[H, BC] = Graph_dynamic(); % Generate vascular graph and boundary conditions

% Compute edge center coordinates
X_Edges = (H.Nodes.X(H.Edges.EndNodes(:,1)) + H.Nodes.X(H.Edges.EndNodes(:,2))) / 2;
Z_Edges = (H.Nodes.Z(H.Edges.EndNodes(:,1)) + H.Nodes.Z(H.Edges.EndNodes(:,2))) / 2;

% Identify penetrating arterioles (PAs) in the specified region
PA = find(X_Edges > 150 & X_Edges < 250 & H.Edges.Type == 2);
PA = PA(2:end); % Remove the first one (assumed to be excluded)

NOC = length(PA); % Number of selected PAs
all = 1:length(H.Edges.D);
No_PA = find(all(~ismember(all, PA))); % Indices of non-PA edges

% Identify capillary segments aligned with the selected PAs
Cap1 = find(X_Edges > 150 & X_Edges < 250 & H.Edges.Type == 6);

% Segment capillaries based on depth
Cap1_L1 = Cap1(H.Edges.ZEdges(Cap1,1) > -210);                           % Superficial layer
Cap1_L2 = Cap1(H.Edges.ZEdges(Cap1,1) > -420 & H.Edges.ZEdges(Cap1,1) <= -210); % Layer 2
Cap1_L3 = Cap1(H.Edges.ZEdges(Cap1,1) > -630 & H.Edges.ZEdges(Cap1,1) <= -420); % Layer 3
Cap1_L4 = Cap1(H.Edges.ZEdges(Cap1,1) <= -630);                          % Deepest layer

%% Conduct initial hemodynamic analysis
Param = input_Boas();                    % Load hemodynamic parameters
BC(1,3) = 40;                            % Set inlet pressure (mmHg)
[q, nodpress, ~, vis] = flow_Boas_new_vis(H, Param, BC(:,1), BC(:,2), BC(:,3), BC(:,4), []);
seg = (nodpress(:, H.Edges.EndNodes(:,1)) + nodpress(:, H.Edges.EndNodes(:,2))) / 2;
% Compute average pressure across each vascular segment

%% Define vessel properties for PA segments
% Estimate relative contractility (ctl) using the pressure of the middle PA segment as reference
H.Edges.ctl = seg' / seg(671);          % Normalized contractility

RC_PA = H.Edges.ctl(PA);                % Relative contractility of selected PA segments

% Assign passive distension properties (arbitrary units) to emulate similar compliance profiles
R_PA_x = linspace(42000, 22000, length(PA)); 

% Define baseline radius and wall thickness assuming zero wall tension (WT = 0)
R_P = 0.91 * H.Edges.D(PA) * 1e-6;       % Baseline radius in meters
h = 0.2 * R_P;                           % Wall thickness as 20% of baseline radius


%% Define mechanotransduction parameters and initialize signaling pathway states
% Time constant for mechanotransdcutiuon pathways: τ_msc = -Δt / ln(1 - α)
% A small alpha leads to slower response;
% In this simulation, we begin with a large alpha (fast adaptation) to stabilize the system
% during the first 20 seconds before applying neurogenic input.

alpha = 2e-3;      % Mechanotransduction gain (used for initial stabilization)
alpha1 = 2e-3;     % Alternative gain (e.g., for other signaling branch if needed)

% After 20 seconds, alphas will be decreased to 0.03e-3 to reflect physiological myogenic dynamics.

% Initial values for signaling pathway intermediates
deact = 0.3;       % Initial deactivation level (e.g., for calcium desensitization)
deact1 = 0.3;      % Optional second pathway
DAG = 0.05;        % Initial Diacylglycerol concentration
IP3 = 0.05;        % Initial IP3 concentration

%% Initialize matrices for storing hemodynamic data over time
Tmax = 80;                    % Total simulation duration [s]
dtp = 0.1;                     % Time step [s]
tp = 0.1:dtp:Tmax+0.1;         % Time vector

% Hemodynamic outputs for each PA segment
PP = zeros(NOC, length(tp)+1);       
qq = zeros(NOC, length(tp)+1);      
WS = zeros(NOC, length(tp));         

% Layer-wise capillary flows
qq_L1 = zeros(1, length(tp)+1);
qq_L2 = zeros(1, length(tp)+1);
qq_L3 = zeros(1, length(tp)+1);
qq_L4 = zeros(1, length(tp)+1);

% Initial values at t = 0
WS(:,1) = vis(PA)' .* 32 / pi .* abs(q(PA)') ./ (H.Edges.D(PA)).^3;  % WSS using Hagen–Poiseuille formula
PP(:,1) = seg(PA);                                                  % Initial segment pressures
qq(:,1) = abs(q(PA));                                               % Initial flow magnitudes

% Average capillary flows per depth layer at t = 0
qq_L1(:,1) = mean(abs(q(Cap1_L1)));
qq_L2(:,1) = mean(abs(q(Cap1_L2)));
qq_L3(:,1) = mean(abs(q(Cap1_L3)));
qq_L4(:,1) = mean(abs(q(Cap1_L4)));

%% Define matrices for solving the segmented vessel model
dt = 1e-4;                          % Time step for ODE solver [s]
t = 0:dt:Tmax;                      % Time vector
t = floor(t * 1e5) / 1e5;           % Avoid numerical round-off error

NOV = 25;                           % Number of ODEs per vessel (state variables per segment)

% Prepare indexing and initial state matrices
id = out_indices(NOC);             % Output indexing structure
idx = indices(NOC, NOV);           % Indexing for internal variables
u = zeros(length(fieldnames(idx)) * length(PA), length(t));     % Internal state matrix
u = initial_conditions(u, idx, R_P);                             % Initialize

x = zeros(length(fieldnames(id)) * length(PA), length(t)-1);    % Output variables (e.g., WSS, Ca2+, etc.)

%% Pressure-dependent boundary conditions and corresponding muscle tone
P_BC = [40 50 60 70 80 90 100 110 120 130];   % Range of ABNP values [mmHg]
MT = [0 -0.091 -0.152 -0.196 -0.2272 -0.249 -0.268 -0.282 -0.291 -0.30]; 
% Corresponding normalized muscle tone

% Set boundary condition to ABNP = 60 mmHg
MTT = -0.152;              % Muscle tone at ABNP = 60 mmHg
BCC = 53.3;                % Adjusted inlet pressure for small cerebrovascular network

% Preallocate arrays for model outputs 
Q_Cap_L1 = zeros(1, length(BCC));
Q_Cap_L2 = zeros(1, length(BCC));
Q_Cap_L3 = zeros(1, length(BCC));
Q_Cap_L4 = zeros(1, length(BCC));

R_PAs = zeros(NOC, length(BCC));         % Radius of PAs [m]
WSS_PAs = zeros(NOC, length(BCC));       % Wall shear stress [Pa]
WSS_std_PAs = zeros(1, length(BCC));     % Standard deviation of WSS
WSS_mean_PAs = zeros(1, length(BCC));    % Mean WSS across PAs

WSS = 0;                                 % Temporary variable (can be used to accumulate)

WT_PAs = zeros(NOC, length(BCC));        % Wall tension [N/m]
P_PAs = zeros(NOC, length(BCC));         % Pressure [mmHg]
V_PAs = zeros(NOC, length(BCC));         % Volume [µm³]
Ca_PAs = zeros(NOC, length(BCC));        % Intracellular calcium concentration [a.u.]

% Preallocate dynamic variables over time (downsampled)
t_d = downsample(t, downsample_factor);  % Downsampled time vector for saving outputs
R_all = zeros(NOC, length(t_d), length(BCC));   % Radius evolution
V_all = zeros(NOC, length(t_d), length(BCC));   % Volume evolution
Ca_all = zeros(NOC, length(t_d), length(BCC));  % Calcium dynamics



formatSpec = '%.2f';



%%
Glu=0;
NOy=0;
%%
for j=1:length(BCC)
    ctr=1;
    MT=MTT(j);
    BC(1,3)=BCC(j);
    % Update diameter of non-PA vessels based on muscle tone and relative contractility
    H.Edges.D(No_PA)=H.Edges.D(No_PA)+MT*(H.Edges.D(No_PA)).*H.Edges.ctl(No_PA);
    for i=2:length(t)
        %% Extract variables from previous time step (t = i-1)
            
            V=u(idx.V, i-1);
            Cl_in=u(idx.Cl_in, i-1);
            BUF_1=u(idx.BUF_1, i-1);
            Na_in=u(idx.Na_in, i-1);
            K_in=u(idx.K_in, i-1);
            Ca_in=u(idx.Ca_in, i-1);
            Ca_NSCstr=u(idx.Ca_NSCstr, i-1);
            BUF_NSCstr1=u(idx.BUF_NSCstr1, i-1);
            NaKCl_PKC=u(idx.NaKCl_PKC, i-1);
            PP_P_NaKCl=u(idx.PP_P_NaKCl, i-1);
            Mp = u(idx.Mp, i-1);
            AMp = u(idx.AMp, i-1);
            AM = u(idx.AM, i-1);
            R= u(idx.R, i-1);
            Ca_juns1=u(idx.Ca_juns1,i-1);
            w_i=u(idx.w_i,i-1);
            Ve=u(idx.Ve,i-1);
            cGMP = u(idx.cGMP, i-1);
            GC0 = u(idx.GC0, i-1);
            GC1 = u(idx.GC1, i-1);
            GC2 = u(idx.GC2, i-1);
            GC1_act = u(idx.GC1_act, i-1);
            PDE_cGMP_P=u(idx.PDE_cGMP_P, i-1);
            CSQ_SRcen1=u(idx.CSQ_SRcen1, i-1);
                    
                    

           %% Time-dependent mechanotransduction adaptation
           % After 20 seconds, reduce alpha to represent slower, physiological myogenic response      
            if (t(i)>=20)  
                alpha = 0.025e-3;
                alpha1 = 0.025e-3;
            end

            % Define time intervals for stimulation phase
            t_k1 = 20;
            t_k2 = 20.5;
            t_k3 = 30;
            t_k4 = Tmax;

            %% EC Extracellular K⁺ dynamics
            %Note: In this simulation, extracellular K⁺ is not actively increased
            if (t(i)<=t_k1)  
                Ke_o=p.Ke_o0;
            elseif (t(i)>= t_k1 & t(i)< t_k2)
                Ke_o=p.Ke_o0+(p.Ke_o4-p.Ke_o0)*(1-exp(-p.k_alphaA1*(t(i)-t_k1)))^2;
            elseif (t(i)>= t_k2 & t(i)< t_k3)
                Ke_o=p.Ke_o4; 
            elseif (t(i)>= t_k3 & t(i)< t_k4)
                Ke_o=p.Ke_o0;           
            end
       
            %% Nitric Oxide (NO) signaling dynamics
            % Simulates spatially delayed NO production across PA segments.
            % The deeper the segment, the earlier it receives NO due to smaller delay.         
            Delay=((Z_Edges(PA)+600)/400)+0.2;
            Delay(end-8:end)=0.2;
            Amp(:)=6.8; % NO peak amplitude
            if(t(i)<p.t_1)
                NO1=p.NO0*ones(length(PA),1).*(t(i)*ones(length(PA),1)<(Delay+p.t_1));
                NO2=p.NO0*ones(length(PA),1).*(t(i)*ones(length(PA),1)<(Delay+p.t_1));
                NO3=p.NO0*ones(length(PA),1).*(t(i)*ones(length(PA),1)<(Delay+p.t_1));
                NO4=p.NO0*ones(length(PA),1).*(t(i)*ones(length(PA),1)<(Delay+p.t_1));
            else
                NO1=p.NO0*ones(length(PA),1).*(t(i)*ones(length(PA),1)<(Delay+p.t_1));
                %             elseif (t(i)>= p.t_1 & t(i)< p.t_7)
                NO2=p.NO0+(Amp-p.NO0).*(t(i)*ones(length(PA),1)>=(Delay+p.t_1)).*(t(i)*ones(length(PA),1)<(Delay+p.t_2)).*(1-exp(-p.k_NO1.*(t(i)-Delay-p.t_1))).^2;
                %             elseif (t(i)>= p.t_2 & t(i)< p.t_8)
                NO3=p.NO0+(Amp-p.NO0).*(t(i)*ones(length(PA),1)>=(Delay+p.t_2)).*(t(i)*ones(length(PA),1)<(Delay+p.t_3)).*(1-exp(-p.k_NO1*(p.t_2-p.t_1))).^2;
                %             elseif (t(i)>=p.t_3)
                NO4=p.NO0+(Amp-p.NO0).*(t(i)*ones(length(PA),1)>=(Delay+p.t_3)).*((1-exp(-p.k_NO1*(p.t_2-p.t_1)))^2).*exp((-p.k_NO2)*(t(i)-Delay-p.t_3));
            end
            NOx=NO1+NO2+NO3+NO4;
            NO=NOx+1.2; % Add basal NO level due to assumed background neuronal activity


            %% Glutamate dynamics
            % Glutamate concentration remains at baseline in this simulation            
            Ampglu(:)=0;            
            Glu=p.NO0+(Ampglu-p.NO0).*(t(i)*ones(length(PA),1)>=(Delay+p.t_1)).*(t(i)*ones(length(PA),1)<(Delay+p.t_3));
           
            %% Hemodynamic analysis        
             
             if(t(i)<tp(1))
                    P=PP(:,1);
                elseif (t(i)>=tp(1) && t(i)<=tp(end))
                    t_i=num2str(t(i),formatSpec);
                    if (ctr<=length(tp))
                        tp_ctr=num2str(tp(ctr),formatSpec);
                    end
                    if(strcmp(t_i,tp_ctr))
                        H.Edges.D(PA)=R*1e6;
                        [q,nodpress,~,vis]=flow_Boas_new_vis(H,Param,BC(:,1),BC(:,2),BC(:,3),BC(:,4),[]);
                        seg=(nodpress(:,H.Edges.EndNodes(:,1)')+nodpress(:,H.Edges.EndNodes(:,2)'))/2;
                        WS(:,ctr)=vis(PA)'.*32/pi.*abs(q(PA)')./((H.Edges.D(PA))).^3;
                        WSS=WS(:,ctr);
                        NOy=WS(:,ctr)./WSS_max;
                        P=seg(PA)';
                        ctr=ctr+1;
                        PP(:,ctr)=P;
                        qq(:,ctr)=q(PA);
                        qq_L1(:,ctr)=mean(abs(q(Cap1_L1)));
                        qq_L2(:,ctr)=mean(abs(q(Cap1_L2)));
                        qq_L3(:,ctr)=mean(abs(q(Cap1_L3)));
                        qq_L4(:,ctr)=mean(abs(q(Cap1_L4)));

                    end

             end

%% SMC fluxes
            %% BK channel
            K_act_i = (Ca_juns1).^2 ./ ((Ca_juns1).^2 + p.beta_i * exp(-(V - p.v_Ca3_i) / p.R_K_i)); 
            if (V==0)
                J_BK=1e9*p.G_K_i*p.perm_BK.*w_i.*(p.K_ext-K_in);
            else
                J_BK=1e9*p.G_K_i*p.perm_BK.*w_i*p.z_K*p.beta.*V.*((p.K_ext-K_in.*exp(p.z_K*p.beta.*V))./(exp(p.z_K*p.beta.*V)-1));
            end
            I_BK=-p.z_K*p.F*J_BK/1000;


            %% SMC kir channel  modulation with WT
            k = 0.025;  % Steepness of the curve
            x0 = WT_max'/2;  % Midpoint of the curve
            % Sigmoid function definition
            sigmoid1 = @(x) 1 ./ (1 + exp(-k * (x - x0)));
            deact_n=(1-(sigmoid1(R.*P./h)));
            deact= alpha1 * deact_n + (1 - alpha1) * deact;

            G_kirbar = p.G_kir/sqrt(p.K_ext);  % [nS/mM^0.5] inward rectifier constant
            E_K= ((1/p.beta)/p.z_K)*log(p.K_ext./K_in);     %[mV]

            I_kir = G_kirbar*(p.K_ext).^p.n_kir .*((V - E_K)./(1 + exp((V - E_K - p.delta_V_kir)/p.k_kir))).*deact;   %[pA] whole cell kir current

            %% Endothelial cells      
            k = 5.25;  % Steepness of the curve
            x0 = WSS_max/2;  % Midpoint of the curve
            % Sigmoid function definition
            sigmoid1 = @(x) 1 ./ (1 + exp(-k * (x - x0)));
            deact_n=(sigmoid1(WS(:,end-2)));
            deact1= alpha1 * deact_n + (1 - alpha1) * deact1;

            % Potassium reversal potential
            E_K = p.RT_F/p.z_K*log(Ke_o./p.Ke_i);     %[mV]
            I_kire = p.G_kirbare.*(Ke_o).^p.n_kire .*((Ve - E_K)./(1 + exp((Ve - E_K - p.delta_V_kire)./p.k_kire))).*deact1;   %[pA] whole cell kir current
                


            %% CaV channel
            CaV=p.CaV_T;
            X_CaV_a=1./(1+exp(-(V-p.V_CaV_act)./(p.zai_CaV_act)));
            Y_CaV_i=1./(1+(Ca_in/p.K_CaV_inh));
            P_CaV=X_CaV_a.*Y_CaV_i;

            if (V==0)
                J_CaV=1e6*p.CaV_T.*p.perm_CaV.*P_CaV.*(p.Ca_ex-Ca_in);
            else
                J_CaV=1e6*p.CaV_T.*p.perm_CaV.*P_CaV*p.z_Ca*p.beta.*V.*((p.Ca_ex-Ca_in.*exp(p.z_Ca.*p.beta.*V))./(exp(p.z_Ca.*p.beta.*V)-1));
            end
            I_CaV=-p.z_Ca*p.F*J_CaV/1000;

            I_CaV_ALL=(1/p.CaV_T)*(CaV.*I_CaV);
            %% Cav 3.2 channels
            CaV3=p.CaV3_T;
            X_CaV_a=1./(1+exp(-(V-p.V_CaV_act)./(p.zai_CaV_act)));
            Y_CaV_i=1./(1+(Ca_in/p.K_CaV_inh));
            P_CaV=X_CaV_a.*Y_CaV_i;

            if (V==0)
                J_CaV3=1e6*p.CaV3_T.*p.perm_CaV.*P_CaV.*(p.Ca_ex-Ca_in);
            else
                J_CaV3=1e6*p.CaV3_T.*p.perm_CaV.*P_CaV*p.z_Ca*p.beta.*V.*((p.Ca_ex-Ca_in.*exp(p.z_Ca.*p.beta.*V))./(exp(p.z_Ca.*p.beta.*V)-1));
            end
            I_CaV3=-p.z_Ca*p.F*J_CaV3/1000;

            I_CaV3=(I_CaV3);
            %% TMEM16A channel
            X_ClA=1./(1+exp(-p.zeta_ClA*p.beta*(V-p.V_ClA)));
            K_ClA_Ca_act=1./(1/(p.K_ClA_Ca_act_max)+X_ClA.*((1/(p.K_ClA_Ca_act_min))-(1/(p.K_ClA_Ca_act_max))));
            P_ClA=((Ca_NSCstr.^p.n_ClA)./((Ca_NSCstr.^p.n_ClA)+(K_ClA_Ca_act.^p.n_ClA))).*((p.K_ClA_Ca_inh^p.n_ClA)./((Ca_NSCstr.^p.n_ClA)+(p.K_ClA_Ca_inh^p.n_ClA)));

            if (V==0)
                J_ClA=1e9*p.ClA_T.*p.perm_ClA.*P_ClA.*(p.Cl_ex-Cl_in);
            else
                J_ClA=1e9*p.ClA_T.*p.perm_ClA.*P_ClA*p.z_Cl*p.beta.*V.*((p.Cl_ex-Cl_in.*exp(p.z_Cl.*p.beta.*V))./(exp(p.z_Cl.*p.beta.*V)-1));
            end

            I_ClA=-p.z_Cl*p.F*J_ClA/1000;

            %% KATP Channel
            KATP=p.KATP_T;
            if (V==0)
                J_KATP=1e9*KATP.*p.perm_KATP.*p.P_open_KATP*(p.K_ext-K_in);
            else
                J_KATP=1e9*KATP.*p.perm_KATP.*p.P_open_KATP*p.z_K*p.beta.*V.*((p.K_ext-K_in.*exp(p.z_K.*p.beta.*V))./(exp(p.z_K.*p.beta.*V)-1));
            end
            I_KATP=-p.z_K*p.F*J_KATP/1000;

            I_KATP_ALL=(I_KATP);

            %% Kv channel
            P_Kv=1./(1+exp((-(V-p.V_Kv))./(p.zeta_Kv)));
            if (V==0)
                J_Kv=1e9*p.Kv_T*p.perm_Kv.*P_Kv.*(p.K_ext-K_in);
            else
                J_Kv=1e9*p.Kv_T*p.perm_Kv.*P_Kv*p.z_K*p.beta.*V.*((p.K_ext-K_in.*exp(p.z_K*p.beta.*V))./(exp(p.z_K*p.beta.*V)-1));
            end
            I_Kv=-p.z_K*p.F*J_Kv/1000;

            %% Cl leak
            if (V==0)
                J_Cl_leak=1e9*p.Cl_leak_T.*p.perm_Cl_leak.*(p.Cl_ex-Cl_in);
            else
                J_Cl_leak=1e9*p.Cl_leak_T.*p.perm_Cl_leak.*p.z_Cl*p.beta.*V.*((p.Cl_ex-Cl_in.*exp(p.z_Cl.*p.beta.*V))./(exp(p.z_Cl.*p.beta.*V)-1));
            end

            I_Cl_leak=-p.z_Cl*p.F*J_Cl_leak/1000;

            %% K leak

            if (V==0)
                J_K_leak=1e9*p.K_leak_T.*p.perm_K_leak.*(p.K_ext-K_in);
            else
                J_K_leak=1e9*p.K_leak_T.*p.perm_K_leak.*p.z_K*p.beta.*V.*((p.K_ext-K_in.*exp(p.z_K.*p.beta.*V))./(exp(p.z_K.*p.beta.*V)-1));
            end

            I_K_leak=-p.z_K*p.F*J_K_leak/1000;

            %% Na leak
            if (V==0)
                J_Na_leak=1e9*p.Na_leak_T.*p.perm_Na_leak.*(p.Na_ex-Na_in);
            else
                J_Na_leak=1e9*p.Na_leak_T.*p.perm_Na_leak.*p.z_Na*p.beta.*V.*((p.Na_ex-Na_in.*exp(p.z_Na.*p.beta.*V))./(exp(p.z_Na.*p.beta.*V)-1));
            end

            I_Na_leak=-p.z_Na*p.F*J_Na_leak/1000;

            %% Ca leak
            if (V==0)
                J_Ca_leak=1e6*p.Ca_leak_T.*p.perm_Ca_leak.*(p.Ca_ex-Ca_in);
            else
                J_Ca_leak=1e6*p.Ca_leak_T.*p.perm_Ca_leak.*p.z_Ca*p.beta.*V.*((p.Ca_ex-Ca_in.*exp(p.z_Ca.*p.beta.*V))./(exp(p.z_Ca.*p.beta.*V)-1));
            end

            I_Ca_leak=-p.z_Ca*p.F*J_Ca_leak/1000;

            %% I_ALL_leak
            I_ALL_leak=I_Cl_leak+I_K_leak+I_Na_leak+I_Ca_leak;


            %% I_NSCne
            P_NSCne=(((p.DAG.^p.n_NSCne)./(p.DAG.^p.n_NSCne+p.K_NSCne_DAG^p.n_NSCne)).*((Ca_in.^p.n_NSCne)./(Ca_in.^p.n_NSCne+p.K_NSCne_Ca_act^p.n_NSCne))).*(1./(1+((Ca_in.^p.n_NSCne)./(p.K_NSCne_Ca_inh^p.n_NSCne))));
            if (V==0)
                J_Na_NSCne=1e9*p.NSCne_T.*p.perm_Na_NSCne*P_NSCne.*(p.Na_ex-Na_in);
            else
                J_Na_NSCne=1e9*p.NSCne_T.*p.perm_Na_NSCne*P_NSCne.*p.z_Na*p.beta.*V.*((p.Na_ex-Na_in.*exp(p.z_Na.*p.beta.*V))./(exp(p.z_Na.*p.beta.*V)-1));
            end

            I_Na_NSCne=-p.z_Na*p.F*J_Na_NSCne/1000;

            if (V==0)
                J_Ca_NSCne=1e6*p.NSCne_T.*p.perm_Ca_NSCne*P_NSCne.*(p.Ca_ex-Ca_in);
            else
                J_Ca_NSCne=1e6*p.NSCne_T.*p.perm_Ca_NSCne*P_NSCne.*p.z_Ca*p.beta.*V.*((p.Ca_ex-Ca_in.*exp(p.z_Ca.*p.beta.*V))./(exp(p.z_Ca.*p.beta.*V)-1));
            end

            I_Ca_NSCne=-p.z_Ca*p.F*J_Ca_NSCne/1000;

            if (V==0)
                J_K_NSCne=1e9*p.NSCne_T.*p.perm_K_NSCne*P_NSCne.*(p.K_ext-K_in);
            else
                J_K_NSCne=1e9*p.NSCne_T.*p.perm_K_NSCne*P_NSCne.*p.z_K*p.beta.*V.*((p.K_ext-K_in.*exp(p.z_K.*p.beta.*V))./(exp(p.z_K.*p.beta.*V)-1));
            end

             I_K_NSCne=-p.z_K*p.F*J_K_NSCne/1000;

            I_NSCne_ALL=(I_Na_NSCne+I_Ca_NSCne+I_K_NSCne);
            %%
            k = 0.01;  % Steepness of the curve
            x0 = WT_max'/2;  % Midpoint of the curve
            % Sigmoid function definition
            sigmoid = @(x) (1 ./ (1 + exp(-k * (x - x0))));  

            nDAG1=0.1*sigmoid(R.*P./h)+0.02;
            nIP31=0.1*sigmoid(R.*P./h)+0.02;

            DAG = alpha * nDAG1 + (1 - alpha) * DAG;
            IP3 = alpha * nIP31 + (1 - alpha) * IP3;

            %% I_TRPM4
            TRPM4=p.TRPM4_T;

            X_TRPM4_act=((Ca_NSCstr.^p.n_TRPM4_Ca)./((Ca_NSCstr.^p.n_TRPM4_Ca)+(p.K_TRPM4_Ca_act^p.n_TRPM4_Ca)));
            X_IP3R=(IP3.^p.n_TRPM4_IP3)./((IP3.^p.n_TRPM4_IP3)+(p.K_TRPM4_IP3.^p.n_TRPM4_IP3));
            P_TRPM4=X_TRPM4_act.*X_IP3R;


            if (V==0)
                J_Na_TRPM4=1e9*TRPM4.*p.perm_Na_TRPM4.*P_TRPM4.*(p.Na_ex-Na_in);
            else
                J_Na_TRPM4=1e9*TRPM4.*p.perm_Na_TRPM4.*P_TRPM4.*p.z_Na*p.beta.*V.*((p.Na_ex-Na_in.*exp(p.z_Na.*p.beta.*V))./(exp(p.z_Na.*p.beta.*V)-1));
            end

            I_Na_TRPM4=-p.z_Na*p.F*J_Na_TRPM4/1000;


            I_Na_TRPM4_ALL=I_Na_TRPM4;


            I_TRPM4_ALL=I_Na_TRPM4_ALL;

            %% I_NSCstr (TRPC6)
            NSCstr6=p.NSCstr_T;

            S_NSCstr_act6=((DAG.^p.n_NSCstr_DAG)./((DAG.^p.n_NSCstr_DAG)+(p.K_NSCstr_DAG^p.n_NSCstr_DAG)));
            X_NSCstr_act6=((Ca_NSCstr.^p.n_NSCstr_Ca)./((Ca_NSCstr.^p.n_NSCstr_Ca)+(p.K_NSCstr_Ca_act^p.n_NSCstr_Ca)));

            Y_NSCstr_inh6=(p.K_NSCstr_Ca_inh^p.n_NSCstr_Ca)./((p.K_NSCstr_Ca_inh^p.n_NSCstr_Ca)+(Ca_NSCstr.^p.n_NSCstr_Ca));
            P_NSCstr6=S_NSCstr_act6.*X_NSCstr_act6.*Y_NSCstr_inh6;



            % Ca I_NSCstr
            if (V==0)
                J_Ca_NSCstr6=1e6*NSCstr6.*p.perm_Ca_NSCstr.*P_NSCstr6.*(p.Ca_ex-Ca_in);
            else
                J_Ca_NSCstr6=1e6*NSCstr6.*p.perm_Ca_NSCstr.*P_NSCstr6.*p.z_Ca*p.beta.*V.*((p.Ca_ex-Ca_in.*exp(p.z_Ca.*p.beta.*V))./(exp(p.z_Ca.*p.beta.*V)-1));
            end

            I_Ca_NSCstr6=-p.z_Ca*p.F*J_Ca_NSCstr6/1000;


            J_Ca_NSCstr_ALL6=J_Ca_NSCstr6;
            I_Ca_NSCstr_ALL6=I_Ca_NSCstr6;



            I_NSCstr_ALL2=I_Ca_NSCstr_ALL6;
%% gap junctions
        I_gj_SMC=1/p.R_gj_SMC*(V-[V(2:end);V(end)]);
        I_gj_EC=1/p.R_gj_SMC*(Ve-[Ve(2:end);Ve(end)]);
        I_gj=1/p.R_gj*(V-Ve);
%% DIffusion of Ca between Compartments
            J_juns_cyt1=1e9*p.lambda_juns_cyt.*(Ca_juns1-Ca_in);
          
            J_juns_cyt_ALL=J_juns_cyt1;
           
            J_CaNSCstr_cyt=1e9*p.lambda_NSCstr_cyt*(Ca_NSCstr-Ca_in); 
%% NGVC, Glutumate receptors Ca influx
            

            P_mG=((Glu.^p.n_NSCstr_DAG)./((Glu.^p.n_NSCstr_DAG)+(p.K_NSCstr_DAG^p.n_NSCstr_DAG)));
            
            if (V==0)
                J_Ca_GlM=1e6*p.Ca_GlR_T.*p.perm_Ca_NSCne*P_mG.*(p.Ca_ex-Ca_in);
            else
                J_Ca_GlM=1e6*p.Ca_GlR_T.*p.perm_Ca_NSCne*P_mG.*p.z_Ca*p.beta.*V.*((p.Ca_ex-Ca_in.*exp(p.z_Ca.*p.beta.*V))./(exp(p.z_Ca.*p.beta.*V)-1));
            end
            
            I_Ca_GlM=-p.z_Ca*p.F*J_Ca_GlM/1000;
%% Transporters
            %% Na,K ATPase
            reldeltau_NaK=(3*p.R*p.T*log(p.Na_ex./Na_in)+(2*p.R*p.T*log(K_in./p.K_ext))-(p.F*V/1000)+p.deltau_ATP)./(p.deltau_ATP);
            I_NaK=p.I_NaK_max*((p.K_ext^p.n_NaK_K)./((p.K_ext^p.n_NaK_K)+(p.K_NaK_Kex^p.n_NaK_K))).*((Na_in.^p.n_NaK_Na)./((Na_in.^p.n_NaK_Na)+(p.K_NaK_Nain^p.n_NaK_Na))).*reldeltau_NaK;
            %% NaKCl 
            NaKCl=p.NaKCl_T-NaKCl_PKC;            
            alpha_NaKCl_0=(p.beta_NaKCl_0*p.alpha_NaKCl_4)/(p.beta_NaKCl_4);
            Q=(p.alpha_NaKCl_4+alpha_NaKCl_0*((p.L_NaKCl_Na*p.L_NaKCl_K*p.L_NaKCl_Cl^2)./(p.Na_ex*p.K_ext*p.Cl_ex^2)))./(p.beta_NaKCl_4+p.beta_NaKCl_0*((p.L_NaKCl_Na*p.L_NaKCl_K*p.L_NaKCl_Cl^2)./(Na_in.*K_in.*Cl_in.^2)));
            D1=1+(p.L_NaKCl_Cl./p.Cl_ex)+((p.L_NaKCl_Cl*p.L_NaKCl_K)./(p.Cl_ex*p.K_ext))+(((p.L_NaKCl_Cl^2)*p.L_NaKCl_K)./((p.Cl_ex^2)*p.K_ext))+(((p.L_NaKCl_Cl^2)*p.L_NaKCl_Na*p.L_NaKCl_K)/((p.Cl_ex^2)*p.Na_ex*p.K_ext));
            D2=1+(p.L_NaKCl_Na./Na_in)+((p.L_NaKCl_Cl*p.L_NaKCl_Na)./(Cl_in.*Na_in))+((p.L_NaKCl_Cl*p.L_NaKCl_K*p.L_NaKCl_Na)./(Cl_in.*K_in.*Na_in))+(((p.L_NaKCl_Cl^2)*p.L_NaKCl_Na*p.L_NaKCl_K)./((Cl_in.^2).*Na_in.*K_in));
            J_4_NaKCl=(NaKCl./(D1+Q.*D2)).*(p.alpha_NaKCl_4-p.beta_NaKCl_4.*Q);
            J_4_NaKCl_PKC=((NaKCl_PKC*p.eps_NaKCl_PKC)./(D1+Q.*D2)).*(p.alpha_NaKCl_4-p.beta_NaKCl_4.*Q);
            J_4_NaKCl_ALL=J_4_NaKCl+J_4_NaKCl_PKC;
            J_x_NaKCl_ALL=J_4_NaKCl;
            J_NaK_Cl=(2e15/p.NA)*J_4_NaKCl_ALL;
            J_KCl_Na=(1e15/p.NA)*J_4_NaKCl_ALL;
            J_NaCl_K=(1e15/p.NA)*J_4_NaKCl_ALL;
            I_NaK_Cl=-p.z_Cl*p.F*J_NaK_Cl/1000;
            %% NCX
            beta_NCX_0=p.alpha_NCX_0*p.beta_NCX_4/p.alpha_NCX_4;
            Z_NCX_act=Ca_in.^p.n_NCX_act./(p.K_NCX_act_Ca^p.n_NCX_act+Ca_in.^p.n_NCX_act);
            QQ=(p.alpha_NCX_4+p.alpha_NCX_0.*(p.K_NCX_Ca*p.K_NCX_Na^3.*(exp((p.F*V)/(2000*p.R*p.T)))./(Ca_in*p.Na_ex^3)))./(p.beta_NCX_4+beta_NCX_0.*(p.K_NCX_Ca*p.K_NCX_Na^3.*(exp((-p.F*V)/(2000*p.R*p.T)))./(p.Ca_ex*Na_in.^3)));
            DD1=1+(p.K_NCX_Ca./Ca_in)+((3*p.K_NCX_Na*p.K_NCX_Ca)./(p.Na_ex*Ca_in))+(((3*p.K_NCX_Na^2)*p.K_NCX_Ca)./((p.Na_ex^2).*Ca_in))+(((p.K_NCX_Na^3)*p.K_NCX_Ca)./((p.Na_ex^3).*Ca_in));
            DD2=1+(p.K_NCX_Ca/p.Ca_ex)+((3*p.K_NCX_Na*p.K_NCX_Ca)./(p.Ca_ex.*Na_in))+((3*p.K_NCX_Na^2*p.K_NCX_Ca)./(p.Ca_ex.*Na_in.^2))+(((p.K_NCX_Na^3)*p.K_NCX_Ca)./((p.Ca_ex).*Na_in.^3));
            S4=(p.alpha_NCX_4-QQ*p.beta_NCX_4).*((Z_NCX_act*p.NCX_T)./(DD1+QQ.*DD2));
            I_NCX=(-1e12*p.F/p.NA)*S4;
            %% PMCA

            reldeltau_PMCA=((-p.F*V/1000)+(p.R*p.T*log(p.Ca_ex./Ca_in))+p.deltau_ATP)/(p.deltau_ATP);
            I_PMCA=p.I_PMCA_max*((Ca_in.^p.n_PMCA_Ca)./((Ca_in.^p.n_PMCA_Ca)+(p.K_PMCA_Ca^p.n_PMCA_Ca))).*reldeltau_PMCA;
            %% SERCA
            SERCA=p.SERCA_T;
            relSERCA=((Ca_in.^p.n_SERCA)./((Ca_in.^p.n_SERCA)+(p.K_SERCA_Ca^p.n_SERCA))).*(1./(1+((p.Ca_SRcen)./(p.K_SERCA_inh_Ca_SRcen)).^p.n_SERCA));
            relSERCA_P=((Ca_in.^p.n_SERCA)./((Ca_in.^p.n_SERCA)+(p.K_SERCA_P_Ca^p.n_SERCA))).*(1./(1+((p.Ca_SRcen)./(p.K_SERCA_inh_Ca_SRcen)).^p.n_SERCA));
            reldeltau_SERCA=((2*p.R*p.T*log(p.Ca_SRcen./Ca_in))+p.deltau_ATP)/(p.deltau_ATP);
            J_SERCA=(-1e3/p.F)*p.I_SERCA_max*relSERCA.*reldeltau_SERCA;
            J_SERCA_P=(-1e3/p.F)*p.I_SERCA_max*relSERCA_P.*reldeltau_SERCA;
            J_SERCA_ALL=(SERCA.*J_SERCA)./p.SERCA_T;
            I_SERCA_ALL=-p.F/1000*(J_SERCA_ALL);
%% ENZYMES


            %% cGMP PDE
            relPDE_cGMP=(cGMP.^p.n_PDE)./((cGMP.^p.n_PDE)+(p.K_PDE_cGMP^p.n_PDE));
            relPDE_cGMP_P=(cGMP.^p.n_PDE)./((cGMP.^p.n_PDE)+(p.K_PDE_cGMP_P^p.n_PDE));
            %% 
            cAMP=0.054;
            relPKA_1=(cAMP.^p.n_PKA_cAMP)./((cAMP.^p.n_PKA_cAMP)+(p.K_PKA_cAMP^p.n_PKA_cAMP));
            relPKC_11=(p.DAG./(p.DAG+p.K_PKC_DAG)).*((Ca_in.^p.n_PKC_Ca)./((Ca_in.^p.n_PKC_Ca)+(p.K_PKC_Ca^p.n_PKC_Ca)));
            relPKCeps_DAG=(p.DAG./(p.DAG+p.K_PKCeps_DAG));
            relPKD_DAG=(p.DAG./(p.DAG+p.K_PKD_DAG));
            relPKG_1=(cGMP.^p.n_PKG_cGMP)./((cGMP.^p.n_PKG_cGMP)+(p.K_PKG_cGMP^p.n_PKG_cGMP));

 %% wall mechanism

            K_1 = p.gamma_cross * Ca_in.^p.n_cross;
            K_6 = K_1;
            K_2 = 58.1395 * p.k_mlcp_b + 58.1395 * p.k_mlcp_c * (NO)*1e0;
            K_5 = K_2;
            M = 1 - AM - AMp - Mp;
            F_r = RC_PA.*(AMp + AM);
            E_Pass=((R.*P*133.322./h)./R_PA_x')*p.E_passive+p.E_passive;
            E = (E_Pass + F_r .* (p.E_active - E_Pass));
            R_0 = R_P + F_r .* (p.alpha - 1) .* R_P;


%% Differential Equations
            PDE_cGMP=p.PDE_cGMP_T-PDE_cGMP_P;
            u(idx.cGMP, i) = u(idx.cGMP ,i-1)+ dt*( p.K_GC_act_GTP*GC1_act+p.K_GC_basal_GTP*(p.GC_T-GC1_act)...
                +(-1)*(p.K_PDE_cGMP_small*cGMP.*((PDE_cGMP.*relPDE_cGMP)+(PDE_cGMP_P.*relPDE_cGMP_P))));
            GC0_P=p.GC_T-GC0-GC1-GC1_act-GC2;
            u(idx.GC0, i) = u(idx.GC0 ,i-1)+ dt*((-p.K_GC0_NO_on*GC0.*NO)+(p.K_GC1_NO_off*GC1)+p.K_PP_GC*GC0_P);
            u(idx.GC1, i) = u(idx.GC1 ,i-1)+ dt*( (p.K_GC0_NO_on*GC0.*NO)-(p.K_GC1_NO_off*GC1)-(p.K_GC1_act*GC1)+(p.K_GC1_deact*GC1_act)-(p.K_GC1_NO_on*NO.*GC1)+(p.K_GC2_NO_off*GC2));
            u(idx.GC2, i) = u(idx.GC2 ,i-1)+ dt*( (p.K_GC1_NO_on*NO.*GC1)-(p.K_GC2_NO_off*GC2)-(p.K_GC2_act*GC2));
            u(idx.GC1_act, i) = u(idx.GC1_act ,i-1)+ dt*((p.K_GC1_act*GC1)-(p.K_GC1_deact*GC1_act)+(p.K_GC2_act*GC2)-(p.K_PKG_GC*relPKG_1.*GC1_act));            

            u(idx.V, i) = u(idx.V ,i-1)+ dt*( (-1e-9/p.C_m)*(I_gj+I_KATP_ALL+I_NSCne_ALL+I_BK+I_CaV_ALL+I_Kv+I_kir+I_PMCA+I_NCX+I_NaK+I_TRPM4_ALL+I_ClA...
                +I_ALL_leak+I_gj_SMC));
            u(idx.Cl_in, i) = u(idx.Cl_in ,i-1)+ dt*( (-1e-9/(p.z_Cl*p.F*p.VOL_cell))*(I_ClA+I_Cl_leak)+(1e-12/p.VOL_cell)*J_NaK_Cl);
            u(idx.BUF_1, i) = u(idx.BUF_1 ,i-1)+ dt*(p.K_BUF_on*Ca_in.*(p.BUF_T-BUF_1)-p.K_BUF_off*BUF_1);

            u(idx.CSQ_SRcen1, i) = u(idx.CSQ_SRcen1 ,i-1)+ dt*(p.K_CSQ_on*p.Ca_SRcen.*(p.CSQ_SRcen_T-CSQ_SRcen1)-p.K_CSQ_off*CSQ_SRcen1);
            u(idx.Na_in, i) = u(idx.Na_in ,i-1)+ dt*((-1e-9/(p.z_Na*p.F*p.VOL_cell))*(I_Na_NSCne+I_Na_TRPM4_ALL+3*I_NCX+3*I_NaK+I_Na_leak)+...
                (1e-12/p.VOL_cell)*J_KCl_Na);
            u(idx.K_in, i) = u(idx.K_in ,i-1)+ dt*( (-1e-9/(p.z_K*p.F*p.VOL_cell))*(I_gj+I_KATP_ALL+I_BK+I_K_NSCne+I_Kv+I_kir-2*I_NaK+I_K_leak)+(1e-12/p.VOL_cell)*J_NaCl_K);
            u(idx.Ca_in, i) = u(idx.Ca_in ,i-1)+ dt*((-1e-6/(p.z_Ca*p.F*p.VOL_cell))*(I_CaV_ALL +I_Ca_NSCne+I_Ca_leak)+(-1e-6/(p.F*p.VOL_cell))*...
                (I_PMCA-I_NCX)+(1e-9/p.VOL_cell)*(J_SERCA_ALL+J_juns_cyt_ALL)-(p.K_BUF_on*Ca_in.*(p.BUF_T-BUF_1)-p.K_BUF_off*BUF_1)...
                +(1e-9/p.VOL_cell)*(J_CaNSCstr_cyt));

            u(idx.p.Ca_SRcen, i) = u(idx.p.Ca_SRcen ,i-1);
            u(idx.Ca_NSCstr, i) = u(idx.Ca_NSCstr ,i-1)+ dt*((1e-9/p.VOL_NSCstr).*(J_Ca_NSCstr_ALL6-J_CaNSCstr_cyt)+(-p.K_BUF_on*Ca_NSCstr.*(p.BUF_NSCstr_T-BUF_NSCstr1)+...
                p.K_BUF_off*BUF_NSCstr1));
            u(idx.BUF_NSCstr1, i) = u(idx.BUF_NSCstr1 ,i-1)+ dt*(p.K_BUF_on*Ca_NSCstr.*(p.BUF_NSCstr_T-BUF_NSCstr1)-p.K_BUF_off*BUF_NSCstr1);
            u(idx.PDE_cGMP_P, i) = u(idx.PDE_cGMP_P ,i-1)+ dt*((p.K_PKA_PDE*relPKA_1+p.K_PKG_PDE*relPKG_1).*PDE_cGMP-p.K_PP_PDE*PDE_cGMP_P);
            PP_NaKCl=p.PP_NaKCl_T-PP_P_NaKCl;
            u(idx.NaKCl_PKC, i) = u(idx.NaKCl_PKC ,i-1)+ dt*(p.K_PKC_NaKCl*relPKC_11.*NaKCl-p.K_PP_NaKCl*(PP_NaKCl./p.PP_NaKCl_T).*NaKCl_PKC-p.K_PP_P_NaKCl*(PP_P_NaKCl./p.PP_NaKCl_T).*NaKCl_PKC);
            u(idx.PP_P_NaKCl, i) = u(idx.PP_P_NaKCl ,i-1)+ dt*(p.K_PKA_PP_NaKCl*relPKA_1.*PP_NaKCl+p.K_PKG_PP_NaKCl*relPKG_1.*PP_NaKCl-p.K_PP_PP_NaKCl*PP_P_NaKCl);
            u(idx.w_i, i) = u(idx.w_i, i-1)+dt*(p.lambda_i * (K_act_i - w_i));
            u(idx.Mp, i) = u(idx.Mp, i-1)+ dt*(p.wallMech * ( p.K_4 * AMp + K_1 .* M - (K_2 + p.K_3) .* Mp ));
            u(idx.AMp, i) = u(idx.AMp, i-1)+ dt*(p.wallMech  * ( p.K_3 * Mp + K_6 .* AM - (p.K_4 + K_5) .* AMp ));
            u(idx.AM, i) = u(idx.AM, i-1)+ dt*(p.wallMech  * ( K_5 .* AMp - (p.K_7 + K_6) .* AM ));
            u(idx.R, i) =  u(idx.R, i-1)+ dt*((R_P / p.eta) .* ( R .* (P*133.322) ./ h - E .* (R - R_0) ./ R_0));

            u(idx.Ca_juns1,i)=u(idx.Ca_juns1,i-1)+ dt*((1e-9./(p.VOL_juns(1))).*(-I_CaV3-J_juns_cyt1+J_Ca_GlM));
            u(idx.Ve,i) = u(idx.Ve,i-1)+ dt*((-1e-9./p.Ce) .* (I_kire-I_gj+I_gj_EC));

            x(id.Glu,i)=Glu;
            x(id.WT,i)=R.*P./h;
            x(id.NO,i)=NO;

                %%
    end
    figure;
    plot(t(10:i-1),u(idx.R(12),10:i-1)./u(idx.R(12),20/dt))
    hold on
    plot(t(10:i-1),u(idx.R(7),10:i-1)./u(idx.R(7),20/dt))
    hold on
    plot(t(10:i-1),u(idx.R(2),10:i-1)./u(idx.R(2),20/dt))
%     save_variables(ABNP(j), u(idx.R,:), u(idx.Ca_in,:), u(idx.V,:), t, PP,tp)
%     R_all(:,:,j) = downsample(u(idx.R,:)', downsample_factor)';
%     V_all(:,:,j) = downsample(u(idx.V,:)', downsample_factor)';
%     Ca_all(:,:,j) = downsample(u(idx.Ca_in,:)', downsample_factor)';
%     
%     Q_Cap_L1(j)=mean(qq_L1(:,end/4:end));
%     Q_Cap_L2(j)=mean(qq_L2(:,end/4:end));
%     Q_Cap_L3(j)=mean(qq_L3(:,end/4:end));
%     Q_Cap_L4(j)=mean(qq_L4(:,end/4:end));
%     R_PAs(:,j)=mean(u(idx.R,end/4:end),2);
%     WSS_PAs(:,j)=mean(WS(:,end/4:end),2);
%     WSS_mean_PAs(j)=mean(WSS_PAs(:,j));
%     WSS_std_PAs(j)=std(WSS_PAs(:,j));
%     WT_PAs(:,j)=mean(x(id.WT,end/4:end),2);
%     P_PAs(:,j)=mean(u(idx.R,end/4:end),2);
%     V_PAs(:,j)=mean(u(idx.V,end/4:end),2);
%     Ca_PAs(:,j)=mean(u(idx.Ca_in,end/4:end),2);
    
end
function id = out_indices(NOC)

    id.Glu=1:NOC;
    id.NO=NOC+1:2*NOC;
    id.WT=2*NOC+1:3*NOC;
end

function idx = indices(NOC,NOV)
    % Index of parameters needing inital conditions 
    idx.cGMP=1:NOV:NOC*NOV;
    idx.GC0 =2:NOV:NOC*NOV+1;
    idx.GC1 = 3:NOV:NOC*NOV+2;
    idx.GC2 = 4:NOV:NOC*NOV+3;
    idx.GC1_act = 5:NOV:NOC*NOV+4;
    idx.V = 6:NOV:NOC*NOV+5;
    idx.Cl_in =7:NOV:NOC*NOV+6;
    idx.BUF_1 = 8:NOV:NOC*NOV+7;
    idx.CSQ_SRcen1 = 9:NOV:NOC*NOV+8;
    idx.Na_in = 10:NOV:NOC*NOV+9;
    idx.K_in = 11:NOV:NOC*NOV+10;
    idx.Ca_in = 12:NOV:NOC*NOV+11;
    idx.p.Ca_SRcen = 13:NOV:NOC*NOV+12;
    idx.Ca_NSCstr = 14:NOV:NOC*NOV+13;
    idx.BUF_NSCstr1 = 15:NOV:NOC*NOV+14;
    idx.PDE_cGMP_P = 16:NOV:NOC*NOV+15;
    idx.NaKCl_PKC = 17:NOV:NOC*NOV+16;
    idx.PP_P_NaKCl = 18:NOV:NOC*NOV+17;
    idx.Mp=19:NOV:NOC*NOV+18;
    idx.AMp=20:NOV:NOC*NOV+19;
    idx.AM=21:NOV:NOC*NOV+20;
    idx.R=22:NOV:NOC*NOV+21;
    idx.Ca_juns1=23:NOV:NOC*NOV+22;
    idx.w_i=24:NOV:NOC*NOV+23;
    idx.Ve=25:NOV:NOC*NOV+24;

end

function u = initial_conditions(u,idx,R_0_passive)

    u(idx.GC0,1) =1e-2; 
    u(idx.GC1,1) =0; 
    u(idx.GC2,1) =0; 
    u(idx.GC1_act,1) =0; 
    u(idx.cGMP,1) =0; 
    u(idx.V,1) =-60; 
    u(idx.Cl_in,1) =50; 
    u(idx.BUF_1,1) =2.6;
    u(idx.CSQ_SRcen1,1) =33 ; 
    u(idx.Na_in,1) =7; 
    u(idx.K_in,1) =140; 
    u(idx.Ca_in,1) =0.1; 
    u(idx.p.Ca_SRcen,1) =100; 
    u(idx.Ca_NSCstr,1) =0.1 ; 
    u(idx.BUF_NSCstr1,1) =2.6; 
    u(idx.PDE_cGMP_P,1) =0;    
    u(idx.PP_P_NaKCl,1) =0; 
    u(idx.Mp,1)=0.0842;
    u(idx.AMp,1)=0.0622;
    u(idx.AM,1)=0.2746;
    u(idx.R,1)=R_0_passive;   
    u(idx.Ca_juns1,1)=0.1;
    u(idx.w_i,1)=0.2206;
    u(idx.Ve,1)=-30;
end


function params = parse_inputs()
    parser = inputParser();
    
    
    parser.addParameter('beta',96487/(1000*8.314*310));
    parser.addParameter('F',96487);
    
    parser.addParameter('t_2_2', 22);
    
    parser.addParameter('t_1', 20); 
    parser.addParameter('t_2', 22); 
    parser.addParameter('t_3', 22);


    
    parser.addParameter('Ke_o0', 3);
    parser.addParameter('Ke_o1', 4);
    parser.addParameter('Ke_o2', 5);
    parser.addParameter('Ke_o3', 3);
    parser.addParameter('Ke_o4', 3);
  
    

    parser.addParameter('k_alphaA1', 10);
    parser.addParameter('k_alphaA2', 10);
   
    
    parser.addParameter('NO0', 0); 
    parser.addParameter('NO1', 0); 
    parser.addParameter('k_NO1', 5);
    parser.addParameter('k_NO2', 0.2); %NO decay rate
    
    
    %%
     
  
    
    %BK
   
    
    parser.addParameter('beta_i', 1); %uM^2
    parser.addParameter('v_Ca3_i', -27); %mV
    parser.addParameter('R_K_i', 12); %mV
    parser.addParameter('G_K_i', 0.6e1); %uM mV^-1 s^-1
    parser.addParameter('v_K_i', -94); %mV
    parser.addParameter('lambda_i', 45); % s^-1
    parser.addParameter('perm_BK', 4e-13);
    parser.addParameter('Ca_GlR_T', 1.53);
    parser.addParameter('z_K', 1);
    
   
    % Kir
    parser.addParameter('delta_V_kir', 20); 
    parser.addParameter('G_kir', 5.5); 
    parser.addParameter('n_kir', 0.5);
    parser.addParameter('k_kir', 11);
    parser.addParameter('a_kir', -0.04);
    parser.addParameter('x_kir', -45);
    
    % EC Kir

    parser.addParameter('delta_V_kire', 4.65); 
    parser.addParameter('G_kirbare', 4.8); 
    parser.addParameter('n_kire', 0.16);
    parser.addParameter('k_kire', 12);

    
    %CaV
    parser.addParameter('CaV_T', 3e3);
    parser.addParameter('CaV3_T', 0.06); 
    parser.addParameter('V_CaV_act', 6.2);
    parser.addParameter('zai_CaV_act', 9.3);
    parser.addParameter('K_CaV_inh', 1.2);
    parser.addParameter('perm_CaV', 1.21e-13);
    parser.addParameter('z_Ca', 2);
    parser.addParameter('V_CaV_PKC_act', 1.2);
    parser.addParameter('V_CaV_PKG_act', 1.12e1);
    parser.addParameter('K_PKC_CaV', 2e1);
    parser.addParameter('K_PKG_CaV', 4e1);
    parser.addParameter('K_PP_CaV', 2e1);
    
  
    
    %CLA
    parser.addParameter('zeta_ClA', 2);
    parser.addParameter('V_ClA', 0);
    parser.addParameter('K_ClA_Ca_act_min', 0.3);
    parser.addParameter('K_ClA_Ca_act_max', 4.1);
    parser.addParameter('n_ClA', 2.6);
    parser.addParameter('K_ClA_Ca_inh', 4e1);
    parser.addParameter('z_Cl', -1);
    parser.addParameter('perm_ClA', 1.6e-14);
    parser.addParameter('ClA_T', 2.5e3);
    
    %G.q
    parser.addParameter('G_qT', 6e3);
    parser.addParameter('K_Gq0_GDP_alpha_AR_on', 2e-5);
    parser.addParameter('K_Gq1_GDP_alpha_AR_off', 1); 
    parser.addParameter('K_Gq1_GDP_GTP_exch', 1);
    parser.addParameter('K_Gq_GTPase', 5e-2);
    
    %G.s
    parser.addParameter('G_sT', 6e3); 
    parser.addParameter('K_Gs_GTPase', 5e-2); 
    parser.addParameter('K_Gs0_GDP_beta_AR_on', 1e-5); 
    parser.addParameter('K_Gs1_GDP_beta_AR_off', 5e-1); 
    parser.addParameter('K_Gs1_GDP_GTP_exch', 1); 


    
    %KATP
    parser.addParameter('KATP_T', 3e2); 
    parser.addParameter('perm_KATP', 7e-14);
    parser.addParameter('P_open_KATP', 5e-3);
    parser.addParameter('P_open_KATP_PKA', 1.2e-1);
    parser.addParameter('K_PKA_KATP', 2);
    parser.addParameter('K_PKCeps_KATP', 2);
    parser.addParameter('K_PP_KATP_PKA', 2e-1);
    parser.addParameter('K_PP_KATP_PKCeps', 2e-1);
    
    %Kv
    parser.addParameter('V_Kv', 6);
    parser.addParameter('zeta_Kv', 14);
    parser.addParameter('perm_Kv', 2.5e-14);
    parser.addParameter('Kv_T', 1.2e3);
    
    %Leaks
    parser.addParameter('Cl_leak_T', 2.2e1);
    parser.addParameter('perm_Cl_leak', 1e-15);
    
    parser.addParameter('Na_leak_T', 0);
    parser.addParameter('perm_Na_leak', 1e-15);
    
    parser.addParameter('K_leak_T', 2e1);
    parser.addParameter('perm_K_leak', 1e-15);
    
    parser.addParameter('Ca_leak_T', 1.3e1);
    parser.addParameter('perm_Ca_leak', 1e-14);
    
    %NaK
    parser.addParameter('R', 8.314);
    parser.addParameter('deltau_ATP', -5e4);
    parser.addParameter('I_NaK_max', 6e1);
    parser.addParameter('n_NaK_K', 1.1);
    parser.addParameter('K_NaK_Kex', 1.6);
    parser.addParameter('K_NaK_Nain', 2.2e1);
    parser.addParameter('n_NaK_Na', 2.5);
    parser.addParameter('z_Na', 1);
    
    %NaKCl
    parser.addParameter('PP_NaKCl_T', 1e4); 
    parser.addParameter('PP_NaKCl_0', 1e4);
    parser.addParameter('NaKCl_T', 1e4); 
    parser.addParameter('beta_NaKCl_0', 5e4);
    parser.addParameter('alpha_NaKCl_4', 5e4);
    parser.addParameter('beta_NaKCl_4', 4e4);
    parser.addParameter('L_NaKCl_Na', 3.2e1);
    parser.addParameter('L_NaKCl_K', 2.7e1);
    parser.addParameter('L_NaKCl_Cl', 6.3e1);
    parser.addParameter('eps_NaKCl_PKC', 2);
    parser.addParameter('NA', 6.022e23);
    parser.addParameter('K_PKC_NaKCl', 2e1);
    parser.addParameter('K_PP_NaKCl', 2e1);
    parser.addParameter('K_PP_P_NaKCl', 4e1);
    parser.addParameter('K_PKA_PP_NaKCl', 2e1);
    parser.addParameter('K_PKG_PP_NaKCl', 2e1);
    parser.addParameter('K_PP_PP_NaKCl', 2e1);
    
    %NCX
    parser.addParameter('alpha_NCX_0', 3e4); 
    parser.addParameter('beta_NCX_4', 5e4);
    parser.addParameter('alpha_NCX_4', 3e4);
    parser.addParameter('n_NCX_act', 2);
    parser.addParameter('K_NCX_act_Ca', 1e-1);
    parser.addParameter('K_NCX_Ca', 8e-2);
    parser.addParameter('K_NCX_Na', 1e2);
    parser.addParameter('NCX_T', 1.6e3);
    

    
    
    %NSCne
    parser.addParameter('n_NSCne', 2);
    parser.addParameter('K_NSCne_DAG', 1.2);
    parser.addParameter('K_NSCne_Ca_act', 2e-1);
    parser.addParameter('K_NSCne_Ca_inh', 1);
    parser.addParameter('NSCne_T', 3.1e1);
    parser.addParameter('perm_Na_NSCne', 3.3e-14);
    parser.addParameter('perm_K_NSCne', 3.3e-14);
    parser.addParameter('perm_Ca_NSCne', 1.5e-13);
    
    %NSCstr
    parser.addParameter('NSC_str_T', 2.5e2); 
    parser.addParameter('n_NSCstr_DAG', 2);
    parser.addParameter('K_NSCstr_DAG', 1.38);
    parser.addParameter('n_NSCstr_IP3', 2);
    parser.addParameter('K_NSCstr_IP3', 1.38);
    parser.addParameter('n_NSCstr_bp', 1);
    parser.addParameter('bp50', 3e2);
    parser.addParameter('n_NSCstr_Ca', 1);
    parser.addParameter('K_NSCstr_Ca_act', 1.7);
    parser.addParameter('K_NSCstr_PKC_Ca_act', 1);
    parser.addParameter('K_NSCstr_PKG_Ca_act', 4);
    parser.addParameter('K_NSCstr_Ca_inh', 1e2);
    parser.addParameter('NSCstr_T', 2.5e2);
    parser.addParameter('perm_Na_NSCstr', 3.3e-14);
    parser.addParameter('perm_Ca_NSCstr', 1.5e-13);
    parser.addParameter('perm_K_NSCstr', 3.3e-14);
    parser.addParameter('K_PKC_NSCstr', 3);
    parser.addParameter('K_PKG_NSCstr', 2e1);
    parser.addParameter('K_PP_NSCstr', 1e1);
    
    %TRPM4
    parser.addParameter('TRPM4_T', 50e2);
    parser.addParameter('n_TRPM4_Ca', 1);
    parser.addParameter('K_TRPM4_Ca_act', 10);
    parser.addParameter('perm_Na_TRPM4', 3.3e-14);
    parser.addParameter('n_TRPM4_DAG', 2);
    parser.addParameter('K_TRPM4_DAG', 1.38);
    parser.addParameter('n_TRPM4_IP3', 2);
    parser.addParameter('K_TRPM4_IP3', 1.38);

    
    %PDE
    parser.addParameter('K_PKA_PDE', 1e1); 
    parser.addParameter('K_PKG_PDE', 1e1); 
    parser.addParameter('K_PP_PDE', 1e1); 
    parser.addParameter('n_PDE', 1); 
    parser.addParameter('PDE_cAMP_T', 1e-2); 
    parser.addParameter('K_PDE_cAMP', 8);
    parser.addParameter('K_PDE_cAMP_P', 3); 
    parser.addParameter('K_PDE_cAMP_small', 1e1);
    parser.addParameter('PDE_cGMP_T', 1e-2); 
    parser.addParameter('K_PDE_cGMP', 8);
    parser.addParameter('K_PDE_cGMP_P', 3);
    parser.addParameter('K_PDE_cGMP_small', 1e1);
    
    %PKA_PKC_PKG_PKCe
    parser.addParameter('n_PKA_cAMP', 2);
    parser.addParameter('K_PKA_cAMP', 4);
    parser.addParameter('K_PKC_DAG', 5e-1);
    parser.addParameter('n_PKC_Ca', 2);
    parser.addParameter('K_PKC_Ca', 2e-1);
    parser.addParameter('K_PKCeps_DAG', 1.7);
    parser.addParameter('K_PKD_DAG', 1.7);
    parser.addParameter('n_PKG_cGMP', 2);
    parser.addParameter('K_PKG_cGMP', 3);
    
    
    %PMCA
    parser.addParameter('n_PMCA_Ca', 3);
    parser.addParameter('K_PMCA_Ca', 4e-1);
    parser.addParameter('I_PMCA_max', 1.8e1);
    

    
    %SERCA
    parser.addParameter('SERCA_T', 1e3); 
    parser.addParameter('n_SERCA', 1.6); 
    parser.addParameter('K_SERCA_Ca', 2e-1);
    parser.addParameter('K_SERCA_inh_Ca_SRcen', 1e2);
    parser.addParameter('K_SERCA_P_Ca', 1e-1);
    parser.addParameter('I_SERCA_max', 8e-1);
    parser.addParameter('K_PKA_SERCA', 2e1);
    parser.addParameter('K_PKC_SERCA', 2e1);
    parser.addParameter('K_PKG_SERCA', 2e1);
    parser.addParameter('K_PP_SERCA', 2e1);
    
    %sGC
    parser.addParameter('GC_T', 1e-2); 
    parser.addParameter('K_GC_act_GTP', 1e1); 
    parser.addParameter('K_GC_basal_GTP', 1e-1); 
    parser.addParameter('K_GC0_NO_on', 3e2); 
    parser.addParameter('K_GC1_NO_off', 2e2); 
    parser.addParameter('K_GC1_act', 1); 
    parser.addParameter('K_GC1_deact', 1e-1); 
    parser.addParameter('K_GC1_NO_on', 2e2); 
    parser.addParameter('K_GC2_NO_off', 2e2);
    parser.addParameter('K_GC2_act', 1e1); 
    parser.addParameter('K_PKG_GC', 1e1);
    parser.addParameter('K_PP_GC', 1e1);
    
    %CAPAC & VOLS
    
    parser.addParameter('C_m', 2.5e-11); 
    parser.addParameter('VOL_cell', 1e-12); 
    parser.addParameter('VOL_SRper', 5.6e-16);
    parser.addParameter('VOL_SRper1', 5.6e-16/5);
    parser.addParameter('VOL_jun', 0.6e-16); 
    parser.addParameter('VOL_juns', (1.3e-16/5));
    parser.addParameter('VOL_SRcen', 7e-14); 
    parser.addParameter('VOL_NSCstr', 6.3e-17); 
    
    % BUFFERS
    parser.addParameter('BUF_T', 3e2); 
    parser.addParameter('K_BUF_on', 2.2e1); 
    parser.addParameter('K_BUF_off', 7.7e1);
    parser.addParameter('BUF_jun_T', 3e2); 
    parser.addParameter('BUF_NSCstr_T', 3e2);
    parser.addParameter('CSQ_SRcen_T', 3e2);
    parser.addParameter('CSQ_SRper_T', 3e2); 
    parser.addParameter('K_CSQ_on', 2);
    parser.addParameter('K_CSQ_off', 1.6e3);
    
    %wall Mechanism 
    parser.addParameter('wallMech', 5);
    
    % Contraction Equation Constants
    parser.addParameter('K_3', 2.4); % s^-1
    parser.addParameter('K_4', 0.1); % s^-1
    parser.addParameter('K_7', 0.1); % s^-1
    parser.addParameter('gamma_cross', 35); %uM^-3 s^-1
    parser.addParameter('n_cross', 3); % fraction constant of the phosphorylation crossbridge
    
    % Mechanical Equation Constants
    parser.addParameter('eta', 1e4); %Pa s
    parser.addParameter('R_0_passive', 14.5e-6); % m
    parser.addParameter('trans_p', 1000); % Pa  transmural pressure
    parser.addParameter('E_passive', 66e3); % Pa
    parser.addParameter('E_active', 433e3); % Pa
    parser.addParameter('alpha', 0.6); % [-]
    parser.addParameter('k_mlcp_b', 0.0086); % [s^-1]
    parser.addParameter('k_mlcp_c', 0.0327); % [s^-1]
    
    % intracell flux
    parser.addParameter('lambda_jun_cyt', 8e-13); 
    parser.addParameter('lambda_juns_cyt', 8e-13/5);
    parser.addParameter('lambda_cen_per', 1.6e-14);
    parser.addParameter('lambda_cen_per1', 1.6e-14/5);
    parser.addParameter('lambda_NSCstr_cyt', 4e-13);

    %ions and T
    parser.addParameter('RT_F', 8314.0*293/96485.0);
    parser.addParameter('T', 3.1e2);
    parser.addParameter('DAG', 0.17);
    parser.addParameter('Na_ex', 1.43e2);
    parser.addParameter('K_ext', 4); 
    parser.addParameter('Ke_i', 150);
    parser.addParameter('Ca_ex', 1.6e3);
    parser.addParameter('Ca_SRcen', 100);
    parser.addParameter('Cl_ex', 1.27e2);
    parser.addParameter('R_gj', 1e0); 
    parser.addParameter('R_gj_SMC', 1e-1); 
    parser.addParameter('Ce', 8e-12);

    
    parser.parse();
    params = parser.Results;
end