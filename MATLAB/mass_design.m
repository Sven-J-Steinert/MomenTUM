clc;
clear;

%==============================================
% mission timeline / deltaV maneuvers
%==============================================

fuel_margine = 1.02; % [1]
sc_duration = 730; % [days] around phobos

% (0) ariane 62 ejection with c3_ariane
% (1) trans-mars-injection to match c3_required
%deltaV1 = sqrt(c3_required) - sqrt(c3_ariane);
deltaV1 = 516.8679 * fuel_margine; % [m/s]
% (2) TCM
deltaV2 = 1.3947 * fuel_margine; % [m/s]
% (3) MOI
deltaV3 = 850.9806 * fuel_margine; % [m/s]
% (4) Match Plane
deltaV4 = 14.5501 * fuel_margine; % [m/s]
% (5) Periapsis raise
deltaV5 = 4.5147 * fuel_margine; % [m/s]
% (6) Circularize
deltaV6 = 713.7529 * fuel_margine; % [m/s]
% (7) Spiral down
deltaV7 = 831.5261 * fuel_margine; % [m/s]
% (8) maintainence @ Phobos
deltaV8 = sc_duration*1.16/8.23 * fuel_margine; % [m/s]
% (9) End-of-life deorbit to 40km altitude periapsis
deltaV9 = 708.3948 * fuel_margine; % [m/s]


sum_deltaV = deltaV1 + deltaV2 + deltaV3 + deltaV4 + deltaV5 + deltaV6 + deltaV7 + deltaV8 + deltaV9; % [m/s]

%=============================================
% Engines
%=============================================

%RCS-Thruster S10-26
dp_bi = 1.5; % [bar] pressure loss
dp_RCS = 1.5; % [bar] pressure loss
p_tank_bi = 20; % [bar]
F_RCS = ((12.5-6)/(23-10)) * (p_tank_bi-dp_RCS) + 1; % [N] thrust depending on tank pressure
m_engine_RCS = 12 * 0.6825; % [kg] redundancy
ISP_RCS = 292; % [s]
OF_RCS = 1.65; % [1]
m_dot_RCS = (((4.2-2.3)/(23-10)) * (p_tank_bi-dp_RCS) + 0.838462) / 1000; % [kg/s] mass flow depending on pressure

%S400-15
F_bi = ((450-340)/(18.5-12.5)) * (p_tank_bi-dp_bi) + 110.833; % [N]
m_engine_bi = 5.46; % [kg]
ISP_bi = (321 * F_bi + 2* ISP_RCS * F_RCS) / (F_bi + 2 * F_RCS);  % [s] mixed ISP
OF_bi = 1.65; % [1]
m_dot_bi = F_bi/(ISP_bi * 9.81); % [kg/s]

%PPS 1350
mode = 1.1;
m_engine_el = 5.04; % [kg]

%PPS 1350-G

if(mode == 1.1)
    
    F_el = 0.090; % [N]
    p_PWR_electric = 1500; % [W]
    ISP_el = 1800; % [s]

elseif(mode == 1.2)
    
    F_el = 0.140; % [N]
    p_PWR_electric = 2500; % [W]
    ISP_el = 1800; % [s]

    
%PPS 1350-E 

elseif(mode == 2.1)
    
    F_el = 0.100; % [N]
    p_PWR_electric = 2500; % [W]
    ISP_el = 1600; % [s]

elseif(mode == 2.2)
    
    F_el = 0.300; % [N]
    p_PWR_electric = 5000; % [W]
    ISP_el = 1900; % [s]

end    
    

%==============================================
% Loopstarter
%==============================================
m_wet = 1; % [kg]
m_wet_1 = 2; % [kg]
m_wet_2 = 1; % [kg]
m_prop = 1; % [kg]
%m_Noble = 1; % [kg]
t_1 = 1; % [s]
t_2 = 1; % [s]
t_3 = 1; % [s]
t_4 = 1; % [s]
t_5 = 1; % [s]
t_6 = 1; % [s]

%==============================================
% power budget
%==============================================

S_Earth = 1361; % [W/m^2]
S_Mars = S_Earth/(1.524^2); % [W/m^2]
eta_solar = 0.29; % [1]
p_PWR_Base = 550; % [W]
p_PWR_tot = p_PWR_Base + p_PWR_electric; % [W]
p_PWR_Earth_tot = p_PWR_tot * 1.524^2 ;
A_solar = p_PWR_tot / (S_Mars * eta_solar) ; % [m^2]


while abs(m_wet_1 - m_wet_2) > 10^-3

    %==============================================
    % propulsion budget
    %==============================================
    ISP_1 = ISP_bi; % [s] R-4D-15 HiPAT™ 445 N
    ISP_2 = ISP_bi; % [s]
    ISP_3 = ISP_bi; % [s]
    ISP_4 = ISP_bi; % [s]
    ISP_5 = ISP_bi; % [s]
    ISP_6 = ISP_el; % [s]
    ISP_7 = ISP_el; % [s] 
    ISP_8 = ISP_el; % [s]
    ISP_9 = ISP_el; % [s]

    m_1_RCS = 2* m_dot_RCS * t_1; % [kg]
    m_2_RCS = 2* m_dot_RCS * t_2; % [kg]
    m_3_RCS = 2* m_dot_RCS * t_3; % [kg]
    m_4_RCS = 2* m_dot_RCS * t_4; % [kg]
    m_5_RCS = 2* m_dot_RCS * t_5; % [kg]
    %m_prop_RCS = (m_1_RCS + m_2_RCS + m_3_RCS + m_4_RCS + m_5_RCS); % [kg] with margine
    %m_MMH_RCS = 1/(OF_RCS+1) * m_prop_RCS;
    %m_MON_3_RCS = m_prop_RCS - m_MMH_RCS;

    roh_MMH = 880; % [kg/m^3]
    roh_MON_3 = 1433.12; % [kg/m^3]
    %m_MMH = 1/(OF_bi+1) * m_dot_bi * (t_1 + t_2 + t_3 + t_4 + t_5); % [kg]
    m_MMH = 1/(OF_bi+1) * 751.275; % [kg]
    %m_MON_3 = m_dot_bi * (t_1 + t_2 + t_3 + t_4 + t_5) - m_MMH; % [kg]
    m_MON_3 = OF_bi * m_MMH; % [kg]
    V_MMH = (m_MMH / roh_MMH)*(1+0.02); % [m^3] pressure regulated
    V_MON_3 = (m_MON_3 / roh_MON_3)*(1+0.02); % [m^3] pressure regulated
    r_i_MMH = (V_MMH / ((4/3)*pi)) ^ (1/3); % [m]
    r_i_MON_3 = (V_MON_3 / ((4/3)*pi)) ^ (1/3); % [m]
    sigma_material = 498*10^6; % [N/m^2] Titan
    roh_material = 4450; % [kg/m^3] Titan
    %m_tank_MMH = 2*pi*r_i_MMH^3*p_tank_bi*10^5*roh_material/sigma_material * (1+0.993); % [kg]
    m_tank_MMH = 22.05; % [kg]
    %m_tank_MON_3 = 2*pi*r_i_MON_3^3*p_tank_bi*10^5*roh_material/sigma_material * (1+0.993); % [kg]
    m_tank_MON_3 = 22.05; % [kg]
    
    roh_Noble = 1723; % [kg/m^3] "Noble density at 187 bar at 65°C" Xenon Tank – COTS
    m_Noble = 166; % [kg]
    V_Noble = (m_Noble / roh_Noble)*(1+0.02); % [m^3] pressure regulated
    r_i_Noble = (V_Noble / ((4/3)*pi)) ^ (1/3); % [m]
    p_tank_el = 187*10^5; % [Pa] = [N/m^2]
    sigma_material_el = 733*10^6; % [N/m^2] CFRP
    roh_material_el = 1800; % [kg/m^3] CFRP
    %m_tank_Noble = 2*pi*r_i_Noble^3*p_tank_el*roh_material_el/sigma_material_el * (1+0.985); % [kg] for Titanium liner
    m_tank_Noble = 20.16; % [kg]

    R_He = 2077.1; % [J/(kg*K)]
    T_He = 293.15; % [K]
    %m_He = (p_tank_bi*10^5) * V_MMH / (R_He *T_He) + (p_tank_bi*10^5) * V_MON_3 / (R_He *T_He); % [kg]
    m_He = 3.5251; % [kg]
    p_tank_He = 310*10^5; % [Pa] = [N/m^2]
    roh_He = p_tank_He / (R_He * T_He); % [kg/m^3] "He density at 310 bar at 20°C" Xenon Tank – COTS 
    V_He = (m_He / roh_He)*(1+0.02); % [m^3] 
    r_i_He = (V_He / ((4/3)*pi)) ^ (1/3); % [m]
    sigma_material_He = 733*10^6; % [N/m^2] CFRP
    roh_material_He = 1800; % [kg/m^3] CFRP
    %m_tank_He = 2*pi*r_i_He^3*p_tank_He*roh_material_He/sigma_material_He * (1+0.752); % [kg] for Titanium liner
    m_tank_He = 8.925; % [kg]

   
    %m_CPROP = (m_engine_bi + m_engine_RCS + m_tank_MMH + m_tank_MON_3 + m_tank_He) * (1+0.171); % [kg] with 26% margine for valves and helium tank and so on
    m_CPROP = m_engine_bi + m_engine_RCS + 1* m_tank_MMH + 2 * m_tank_MON_3 + 2 * m_tank_He + 29.29848; % [kg]
    %m_EPROP = (m_engine_el + m_tank_Noble) * (1+1.05); % [kg]
    m_EPROP = m_engine_el + m_tank_Noble + 20.16; % [kg]
    m_PROP = m_EPROP + m_CPROP; % [kg]
    
    %==============================================
    % mass budget
    %==============================================
    m_INS = 50; % [kg]
    m_AOGNC = 40; % [kg]
    m_COM = 25; % [kg]
    m_CDH = 15; % [kg]
    m_MEC = 55; % [kg]
    m_PWR_el_only = 0.0758 * p_PWR_electric * 1.524^2; % [kg]
    m_PWR = 0.0758 * p_PWR_Earth_tot + 37.741; % [kg]
    m_TC = 15; % [kg]
    m_STR = 0.0224*m_wet+27.823; % [kg]
    m_Probe = 20; % [kg]
    m_dry_no_Harness = m_INS + m_AOGNC + m_COM + m_CDH + m_MEC + m_PWR + m_TC + m_STR + m_Probe + m_PROP; % [kg]
    m_Harness = m_dry_no_Harness * 0.1; % [kg]
    m_dry = (m_dry_no_Harness + m_Harness) * (1+0.2); % [kg] with 20% margine
    m_pres = m_He * 1.02; % [kg] with margine
    m_wet_1 = m_dry + m_prop; % [kg]
    m_wet = m_wet_1 + m_pres; % [kg]
    m_adapter = 85; % [kg]
    m_launch = m_wet + m_adapter; % [kg]

   
    %==============================================
    % propellant budget
    %==============================================
    m_final_9 = m_dry - m_Probe; % [kg]
    m_0_9 = m_final_9 * exp((deltaV9/(9.81*ISP_9))); % [kg]
    m_prop_9 = m_0_9 - m_final_9; % [kg]
    m_dot_prop9 = F_el/(ISP_9 * 9.81); % [kg/s]
    t_9 = (m_prop_9/m_dot_prop9); %[days]
    
    m_final_8 = m_0_9; % [kg]
    m_0_8 = m_final_8 * exp((deltaV8/(9.81*ISP_8))); % [kg]
    m_prop_8 = m_0_8 - m_final_8; % [kg]
    m_dot_prop8 = F_el/(ISP_8 * 9.81); % [kg/s]
    t_8 = (m_prop_8/m_dot_prop8); %[days]

    m_final_7 = m_0_8 + m_Probe; % [kg]
    m_0_7 = m_final_7 * exp((deltaV7/(9.81*ISP_7))); % [kg]
    m_prop_7 = m_0_7 - m_final_7; % [kg]
    m_dot_prop7 = F_el/(ISP_7 * 9.81); % [kg/s]
    t_7 = (m_prop_7/m_dot_prop7); %[days]

    m_final_6 = m_0_7; % [kg]
    m_0_6 = m_final_6 * exp((deltaV6/(9.81*ISP_6))); % [kg]
    m_prop_6 = m_0_6 - m_final_6; % [kg]
    m_dot_prop6 = F_el/(ISP_6 * 9.81); % [kg/s]
    t_6 = (m_prop_6/m_dot_prop6); %[days]

    m_final_5 = m_0_6; % [kg]
    m_0_5 = m_final_5 * exp((deltaV5/(9.81*ISP_5))); % [kg]
    m_prop_5 = m_0_5 - m_final_5; % [kg]
    m_dot_prop5 = F_bi/(ISP_5 * 9.81); % [kg/s]
    t_5 = (m_prop_5/m_dot_prop5); %[days]

    m_final_4 = m_0_5; % [kg]
    m_0_4 = m_final_4 * exp((deltaV4/(9.81*ISP_4))); % [kg]
    m_prop_4 = m_0_4 - m_final_4; % [kg]
    m_dot_prop4 = F_bi/(ISP_4 * 9.81); % [kg/s]
    t_4 = (m_prop_4/m_dot_prop4); %[days]

    m_final_3 = m_0_4; % [kg]
    m_0_3 = m_final_3 * exp((deltaV3/(9.81*ISP_3))); % [kg]
    m_prop_3 = m_0_3 - m_final_3; % [kg]
    m_dot_prop3 = F_bi/(ISP_3 * 9.81); % [kg/s]
    t_3 = (m_prop_3/m_dot_prop3); %[days]

    m_final_2 = m_0_3; % [kg]
    m_0_2 = m_final_2 * exp((deltaV2/(9.81*ISP_2))); % [kg]
    m_prop_2 = m_0_2 - m_final_2; % [kg]
    m_dot_prop2 = F_bi/(ISP_2 * 9.81); % [kg/s]
    t_2 = (m_prop_2/m_dot_prop2); %[days]

    m_final_1 = m_0_2; % [kg]
    m_0_1 = m_final_1 * exp((deltaV1/(9.81*ISP_1))); % [kg]
    m_prop_1 = m_0_1 - m_final_1; % [kg]
    m_dot_prop1 = F_bi/(ISP_1 * 9.81); % [kg/s]
    t_1 = (m_prop_1/m_dot_prop1); %[days]

    m_wet_2 = m_0_1; % [kg]
    m_prop = m_prop_9 + m_prop_8 + m_prop_7 + m_prop_6 + m_prop_5 + m_prop_4 + m_prop_3 + m_prop_2 + m_prop_1; % [kg]
    %m_Noble = m_prop_9 + m_prop_8 + m_prop_7 + m_prop_6;

end

m_prop_chem = m_prop_5 + m_prop_4 + m_prop_3 + m_prop_2 + m_prop_1; % [kg]
m_prop_el = m_prop_9 + m_prop_8 + m_prop_7 + m_prop_6; % [kg]
m_prop_RCS = 0.05 * m_prop_chem; % [kg]

disp(["Launch mass [kg]",m_launch, "< 2600 kg"]);
disp(["time limit station keeping", deltaV8/t_8 *86400,">= 0.1409 dV/day"])
disp(["D_V1 to D_V9 [m/s]",deltaV1, deltaV2, deltaV3, deltaV4, deltaV5, deltaV6, deltaV7, deltaV8, deltaV9, "Summ_Delta_V [m/s]", sum_deltaV]);
disp(["Dry mass [kg]", m_dry]);
disp(["PROP system mass dry [kg]", m_PROP]);
disp(["CPROP system mass dry [kg]", m_CPROP]);
disp(["EPROP system mass dry [kg]", m_EPROP]);
disp(["MMH tank mass [kg]", m_tank_MMH, "MON_3 tank mass [kg]", m_tank_MON_3, "Noble tank mass [kg]", m_tank_Noble, "He tank mass [kg]", m_tank_He]);
disp(["MMH volume [l]", V_MMH*1000, "MON_3 volume [l]", V_MON_3*1000,"Noble volume [l]", V_Noble*1000, "He volume [l]", V_He*1000]);
disp(["MMH mass [kg]", m_MMH, "MON_3 mass [kg]", m_MON_3, "Noble mass [kg]", m_Noble, "He mass [kg]", m_He]);
disp(["m_0 [kg]",m_0_1, m_0_2, m_0_3, m_0_4, m_0_5, m_0_6, m_0_7, m_0_8, m_0_9]);
disp(["m_PWR_el_only [kg]",m_PWR_el_only]);
disp(["Burn time [days]", t_1/(24*3600), t_2/(24*3600), t_3/(24*3600), t_4/(24*3600), t_5/(24*3600), t_6/(24*3600), t_7/(24*3600), t_8/(24*3600), t_9/(24*3600)]);
disp(["Solar array size [m^2]", A_solar]);
disp(["RCS prop mass [kg]", m_prop_RCS]);
disp(["ISP", ISP_bi]);
disp(["m_prop_8 [kg]", m_prop_8]);
disp(["cprop [kg]", m_prop_chem]);
disp(["eprop [kg]", m_prop_el]);
disp(["F_bi [N]", F_bi]);
disp(["F_RCS [N]", F_RCS]);
