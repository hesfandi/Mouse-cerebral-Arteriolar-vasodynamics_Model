function params = parse_inputs(n,NOC)
    parser = inputParser();
    
    parser.addParameter('beta',96487/(1000*8.314*310));
    parser.addParameter('F',96487);
    %% effectors
%     parser.addParameter('t_1', 60); 
%     parser.addParameter('t_2', 63); 
%     parser.addParameter('t_3', 90);
%     parser.addParameter('t_4', 93); 
%     parser.addParameter('t_5', 120); 
%     parser.addParameter('t_6', 123); 
%     parser.addParameter('t_7', 150); 
%     parser.addParameter('t_8', 153); 
%     parser.addParameter('t_9', 180); 
%     parser.addParameter('t_10', 183); 
%     parser.addParameter('t_11', 210); 
%     parser.addParameter('t_12', 213); 
%     parser.addParameter('t_13', 240); 
%     parser.addParameter('t_14', 243); 
%     parser.addParameter('t_15', 270); 
%     parser.addParameter('t_16', 273); 
%     parser.addParameter('t_17', 300); 
    
    
    parser.addParameter('t_1', 20); 
    parser.addParameter('t_2', 23); 
    parser.addParameter('t_3', 40);
    parser.addParameter('t_4', 41); 
    parser.addParameter('t_5', 50); 
    parser.addParameter('t_6', 2.1); 
    parser.addParameter('t_7', 3.5); 
    parser.addParameter('t_8', 4); 
    
    
    parser.addParameter('EC0', -35); 
    parser.addParameter('EC1', -70);
    
    
    parser.addParameter('P0', 1); 
    parser.addParameter('P1', 40);
    parser.addParameter('P2', 50);
    parser.addParameter('P3', 30);
    parser.addParameter('P4', 40);
    parser.addParameter('P5', 50);
    parser.addParameter('P6', 60);
    parser.addParameter('P7', 70);
    parser.addParameter('P8', 80);
    parser.addParameter('P9', 90);
    parser.addParameter('P10', 100);
  
    
    parser.addParameter('alphaA0', 0); 
    parser.addParameter('alphaA1', 0); 
    parser.addParameter('k_alphaA1', 10);
    parser.addParameter('k_alphaA2', 10);
    
    parser.addParameter('ATP0', 0); 
    parser.addParameter('ATP1', 0); 
    parser.addParameter('k_ATP1', 2);
    parser.addParameter('k_ATP2', 2);
    parser.addParameter('t_pulse_init', 100);
    parser.addParameter('deltat_pulse_on', 1);
    parser.addParameter('deltat_pulse_off', 9);
    parser.addParameter('n_cycle', 20);
    
    parser.addParameter('betaA0', 0); 
    parser.addParameter('betaA1', 0); 
    parser.addParameter('k_betaA1', 2);
    parser.addParameter('k_betaA2', 2);
    
    parser.addParameter('BP0', 0); 
    parser.addParameter('BP1', 0); 
    parser.addParameter('k_BP1', 2);
    parser.addParameter('k_BP2', 2);
    parser.addParameter('t_BP1', 10); 
    parser.addParameter('t_BP2', 15); 
    parser.addParameter('t_BP3', 300); 
    
    parser.addParameter('EET0', 0); 
    parser.addParameter('EET1', 0); 
    parser.addParameter('k_EET1', 2);
    parser.addParameter('k_EET2', 2);
    parser.addParameter('t_EET1', 40); 
    parser.addParameter('t_EET2', 41); 
    parser.addParameter('t_EET3', 55); 
    
    parser.addParameter('NO0', 0); 
    parser.addParameter('NO1', 7); 
    parser.addParameter('k_NO1', 10);
    parser.addParameter('k_NO2', 10);
    
    
    %%
    % alphaAR
    parser.addParameter('K_alpha_AR_on', 1); 
    parser.addParameter('K_alpha_AR_off', 1e1); 
    parser.addParameter('alpha_AR_T', 6e3); 
    parser.addParameter('K_alpha_AR_P', 2e-2); 
    parser.addParameter('K_alpha_AR_endo', 2e-2); 
    parser.addParameter('K_alpha_AR_recyc', 1e-2); 
    
    %AC
    parser.addParameter('AC_T', 6e3); 
    parser.addParameter('K_PKA_AC', 1e-2); 
    parser.addParameter('K_PP_AC', 1e-2); 
    parser.addParameter('K_AC0_Gs_on', 4e-5); 
    parser.addParameter('K_AC1_Gs_off', 5e-1);
    parser.addParameter('K_AC1_ATP', 2e2);
    parser.addParameter('Nu', 6.022e23*1e-12);
    parser.addParameter('K_AC0_ATP', 1e-1);
     
  
    % beta_AR        
    parser.addParameter('beta_AR_T', 6e3); 
    parser.addParameter('K_beta_AR_on', 2e-1); 
    parser.addParameter('K_beta_AR_off', 4); 
    parser.addParameter('K_beta_AR_P', 4e-3); 
    parser.addParameter('K_beta_AR_endo', 4e-3); 
    parser.addParameter('K_beta_AR_recyc', 4e-3); 
    
    %BK
    parser.addParameter('BK_T', 0.75e3/n); 
    parser.addParameter('z_BK_J', 5.7e-1);
    parser.addParameter('V_BK_closed', 8e1);
    parser.addParameter('V_BK_open', -3.4e1);
    parser.addParameter('K_BK_closed', 4.72);
    parser.addParameter('K_BK_open', 8e-1);
    parser.addParameter('L_BK', 2.5e-6);
    parser.addParameter('L_BK_PKA', 7.5e-6);
    parser.addParameter('L_BK_PKC', 8.33e-7);
    parser.addParameter('L_BK_PKG', 7.5e-6);
    parser.addParameter('z_BK_L', 4.1e-1);
    parser.addParameter('perm_BK', 4e-13);
    parser.addParameter('z_K', 1);
    parser.addParameter('K_PKA_BK', 2e1);
    parser.addParameter('K_PKC_BK', 2e1);
    parser.addParameter('K_PKG_BK', 4e1);
    parser.addParameter('K_PP_BK', 2e1);
    
    % Kir
    parser.addParameter('delta_V_kir', 28.89); 
    parser.addParameter('G_kir', 0.5); 
    parser.addParameter('n_kir', 0.5);
    parser.addParameter('k_kir', 10);
    parser.addParameter('a_kir', -0.04);
    parser.addParameter('x_kir', -45);
    
    %CaV
    parser.addParameter('CaV_T', 3e3); 
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
    
    %CaV3.2 
    parser.addParameter('CaV3_T', 0.07e0); 
    parser.addParameter('V_Cav3_b', -28); 
    parser.addParameter('Cav3_b', 6.1);
    parser.addParameter('Cav_coef_b', 1.068);
    parser.addParameter('V_Cav3_b1',-16.3);
    parser.addParameter('V_coef_Cav3_b', 30);
    parser.addParameter('V_Cav3_g', -60);
    parser.addParameter('Cav3_g', 6.6);
    parser.addParameter('Cav_coef_g', 0.015);
    parser.addParameter('V_Cav3_g1', -71.7);
    parser.addParameter('V_coef_Cav3_g1', 83.33);
    parser.addParameter('V_coef_Cav3_g2', 15.38);
    
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

    %IP3R
    parser.addParameter('IP3R_T', 2e3);
    parser.addParameter('IP3R_0', 2e3);
    parser.addParameter('n_IP3R_IP3', 3);
    parser.addParameter('K_IP3R_IP3', 3.5e-1);
    parser.addParameter('K_IP3R_IRAG_PKG_IP3', 8e-1);
    parser.addParameter('K_IP3R_Ca_max', 1);
    parser.addParameter('K_IP3R_Ca_min', 1e-1);
    parser.addParameter('K_IP3R_Ca_inh', 3e-1);
    parser.addParameter('n_IP3R_Ca_inh', 2.5);
    parser.addParameter('n_IP3R_CaSR', 2.5);
    parser.addParameter('K_IP3R_CaSR', 1.4e2);
    parser.addParameter('n_IP3R_Ca_act', 3);
    parser.addParameter('perm_IP3R', 1.1e-13);
    parser.addParameter('K_PKG_IP3R_IRAG', 2e1);
    parser.addParameter('K_PP_IP3R_IRAG_PKG', 2e1);
    
    parser.addParameter('IP3R1_T', 2e3);
    parser.addParameter('K_on_IP3', 1.4e-6);
    parser.addParameter('K_inh_IP3', 100);
    parser.addParameter('K_IP3', 0.12);
%     parser.addParameter('IP3_sense', 2e1);
    parser.addParameter('K_act_IP3', 0.1);

    
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
    parser.addParameter('V_Kv', -2e1);
    parser.addParameter('zeta_Kv', 7);
    parser.addParameter('perm_Kv', 2.5e-14);
    parser.addParameter('Kv_T', 2.7e2);
    
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
    
    %NSCeet
    parser.addParameter('V_NSCeet_min', -1e1);
    parser.addParameter('V_NSCeet_max', 50);
    parser.addParameter('K_NSCeet_EET', 5e-2);
    parser.addParameter('n_NSCeet', 1);
    parser.addParameter('zeta_NSCeet', 1e1);
    parser.addParameter('NSCeet', 4e1);
    parser.addParameter('perm_Na_NSCeet', 1.3e-13);
    parser.addParameter('perm_Ca_NSCeet', 9e-13);
    parser.addParameter('perm_K_NSCeet', 1.3e-13);
    
    
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
    parser.addParameter('NSC_str_T', 2.61e2); 
    parser.addParameter('n_NSCstr_DAG', 2);
    parser.addParameter('K_NSCstr_DAG', 1.38);
    parser.addParameter('n_NSCstr_IP3', 2);
    parser.addParameter('K_NSCstr_IP3', 1.38);
    parser.addParameter('n_NSCstr_bp', 1);
    parser.addParameter('bp50', 3e2);
    parser.addParameter('n_NSCstr_Ca', 1);
    parser.addParameter('K_NSCstr_Ca_act', 1.7);
    parser.addParameter('K_TRPM4_Ca_act', 10);
    parser.addParameter('K_NSCstr_PKC_Ca_act', 1);
    parser.addParameter('K_NSCstr_PKG_Ca_act', 4);
    parser.addParameter('K_NSCstr_Ca_inh', 1e2);
    parser.addParameter('NSCstr_T', 2.61e2);
    parser.addParameter('perm_Na_NSCstr', 3.3e-14);
    parser.addParameter('perm_Ca_NSCstr', 1.5e-13);
    parser.addParameter('perm_K_NSCstr', 3.3e-14);
    parser.addParameter('K_PKC_NSCstr', 3);
    parser.addParameter('K_PKG_NSCstr', 2e1);
    parser.addParameter('K_PP_NSCstr', 1e1);
    
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
    parser.addParameter('K_PDE_cGMP_small', 1e3);
    
    %PKA_PKC_PKG_PKCe
    parser.addParameter('n_PKA_cAMP', 2);
    parser.addParameter('K_PKA_cAMP', 4);
    parser.addParameter('K_PKC_DAG', 5e-1);
    parser.addParameter('n_PKC_Ca', 2);
    parser.addParameter('K_PKC_Ca', 2e-1);
    parser.addParameter('K_PKCeps_DAG', 1.7);
    parser.addParameter('K_PKD_DAG', 1.7);
    parser.addParameter('n_PKG_cGMP', 1);
    parser.addParameter('K_PKG_cGMP', 0.4);
    parser.addParameter('K_PKG_TRPC6', 0.7);
    %PLC
    parser.addParameter('PLC_T', 6e3); 
    parser.addParameter('K_PLC_PIP_on', 1e-6); 
    parser.addParameter('K_PLC_PIP_off', 1e2);
    parser.addParameter('K_PLC_PIP_hyd', 2e1); 
    parser.addParameter('K_PLC_Gq0_GTP_on', 2e-4);
    parser.addParameter('K_PLC_Gq0_GTP_off', 2e-2); 
    parser.addParameter('K_PLC_Gq0_GTP_hyd', 1e-1);
    parser.addParameter('K_PLC_Gq0_GTP_PIP_on', 3e-6); 
    parser.addParameter('K_PLC_Gq0_GTP_PIP_off', 1e1);
    parser.addParameter('K_PLC_Gq0_GTP_PIP_hyd', 8e1); 
    parser.addParameter('K_PKA_PLC_Gq0_GTP', 1e1); 
    parser.addParameter('K_PKC_PLC_Gq0_GTP', 1e1);
    parser.addParameter('K_PLC_P_Gq0_GTP_on', 5e-6); 
    parser.addParameter('K_PLC_P_Gq0_GTP_off', 1);
    parser.addParameter('K_PP_PLC', 2); 
    parser.addParameter('PIP_T', 1e7);
    parser.addParameter('K_met_DAG', 0.09); 
    parser.addParameter('K_met_IP3', 0.09);
    
    %PMCA
    parser.addParameter('n_PMCA_Ca', 3);
    parser.addParameter('K_PMCA_Ca', 4e-1);
    parser.addParameter('I_PMCA_max', 1.8e1);
    
    %P2XR
    parser.addParameter('P2XR_T', 1e2);
    parser.addParameter('perm_NaP2XR', 4.3e-14);
    parser.addParameter('perm_KP2XR', 4.3e-14);
    parser.addParameter('perm_CaP2XR', 1.3e-13);
    parser.addParameter('K_P2XR_0_on', 1e1);
    parser.addParameter('K_P2XR_3_off', 1e2);
    parser.addParameter('n_P2XR', 2);
    parser.addParameter('K_P2XR_3_act', 1.2);
    parser.addParameter('K_P2XR_3_deact', 1e1);
    parser.addParameter('K_P2XR_3_act_desens', 1);
    parser.addParameter('K_P2XD_3_resens', 1e-2);
    parser.addParameter('K_P2XR_0_desens', 1e-1);
    parser.addParameter('K_P2XD_0_resens', 1);
    parser.addParameter('K_P2XD_0_on', 1e1);
    parser.addParameter('K_P2XD_3_off', 1e-2);
    parser.addParameter('K_P2XD_0_PKD_on', 1);
    parser.addParameter('K_P2XD_3_PKD_off', 1);
    parser.addParameter('K_PKD_P2XD_3', 1);
    parser.addParameter('K_PP_P2XD_0_PKD', 1);
    
    %RyRs    
    parser.addParameter('RyR_T', 3e3);
    parser.addParameter('RyRs_T', 3e3/n); % 
    parser.addParameter('perm_RyR', 2e-13);
    parser.addParameter('K_RyR_Ca_SR', 8e2);
    parser.addParameter('n_RyR_Ca_SR', 5e-1);
    parser.addParameter('K_RyR_Ca_jun_max', 1e1);
    parser.addParameter('K_RyR_Ca_jun_min', 4);
    parser.addParameter('n_RyR_Ca_jun', 2);
    parser.addParameter('K_RyR_Ca_inh', 3);
    parser.addParameter('n_RyR_Ca_inh', 3);
    parser.addParameter('K_PKA_RyR', 4e1);
    parser.addParameter('K_PKC_RyR', 4e1);
    parser.addParameter('K_PP_RyR', 4e1);
    parser.addParameter('f_RyR_PKA_min', 4e-1);
    parser.addParameter('f_RyR_PKC_max', 5);
    
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
    parser.addParameter('K_GC_act_GTP', 1e2); 
    parser.addParameter('K_GC_basal_GTP', 1e-1); 
    parser.addParameter('K_GC0_NO_on', 9e2); 
    parser.addParameter('K_GC1_NO_off', 2e2); 
    parser.addParameter('K_GC1_act', 10); 
    parser.addParameter('K_GC1_deact', 1e-1); 
    parser.addParameter('K_GC1_NO_on', 2e3); 
    parser.addParameter('K_GC2_NO_off', 2e2);
    parser.addParameter('K_GC2_act', 1e3); 
    parser.addParameter('K_PKG_GC', 1e3);
    parser.addParameter('K_PP_GC', 1e3);
    
    %CAPAC & VOLS
    
    parser.addParameter('C_m', 1.5e-11); 
    parser.addParameter('VOL_cell', 1e-12); 
    parser.addParameter('VOL_SRper', 5.6e-16);
    parser.addParameter('VOL_SRper1', 5.6e-16/n);
    parser.addParameter('VOL_jun', 0.6e-16); 
    parser.addParameter('VOL_juns', (1.3e-16/n)+0.5*(1.3e-16/n)*(-1+2*rand(n*NOC,1)));
%     parser.addParameter('VOL_juns', (1.3e-16/n)*ones(n,1));
    parser.addParameter('VOL_SRcen', 7e-14); 
    parser.addParameter('VOL_NSCstr', 6.3e-17); 
    parser.addParameter('R_gj', 1e0); 
    
    
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
    parser.addParameter('wallMech', 8);
    
    % Contraction Equation Constants
    parser.addParameter('K_3', 0.4); % s^-1
    parser.addParameter('K_4', 0.1); % s^-1
    parser.addParameter('K_7', 0.1); % s^-1
    parser.addParameter('gamma_cross', 25); %uM^-3 s^-1
    parser.addParameter('n_cross', 3); % fraction constant of the phosphorylation crossbridge
    
    % Mechanical Equation Constants
    parser.addParameter('eta', 1e4); %Pa s
    parser.addParameter('R_0_passive', 16e-6); % m
    parser.addParameter('trans_p', 1000); % Pa  transmural pressure
    parser.addParameter('E_passive', 66e3); % Pa
    parser.addParameter('E_active', 433e3); % Pa
    parser.addParameter('alpha', 0.6); % [-]
    parser.addParameter('k_mlcp_b', 0.0086); % [s^-1]
    parser.addParameter('k_mlcp_c', 0.0327); % [s^-1]
    
    % intracell flux
    parser.addParameter('lambda_jun_cyt', 8e-13); 
    parser.addParameter('lambda_juns_cyt', 8e-13/n);
%     parser.addParameter('lambda_juns_cyt', (8e-13/n)+0.1*(8e-13/n)*(-1+2*rand(n,1)));
    parser.addParameter('lambda_cen_per', 1.6e-14);
    parser.addParameter('lambda_cen_per1', 1.6e-14/n);
    parser.addParameter('lambda_NSCstr_cyt', 4e-13);

    %ions and T
    parser.addParameter('T', 3.1e2);
    parser.addParameter('Na_ex', 1.43e2);
    parser.addParameter('K_ex', 5.9); 
    parser.addParameter('Ca_ex', 1.6e3);
    parser.addParameter('Cl_ex', 1.27e2);
    
    parser.parse();
    params = parser.Results;
end