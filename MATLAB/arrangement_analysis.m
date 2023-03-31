clc;
clear;

%%FRAME in an equatorial orbit%%

% X = axis of spacecraft with engines, oriented to Mars/Phobos gravicenter with the
% lowest moment of inertia

%Y = solar panels left and right of the X-axis 

%Z = axis of constant rotation with the highest inertia

%Satellite box
l_sat = [4; 1.3529; 1.3529]; % [m] dimensions of half the satellite box (box from the tank outline + wthat RCS needs)
V_sat = l_sat(1) * l_sat(2)*2 * l_sat(3)*2; % [m^3]

%dry_mass
m_dry = 1141.5902; % [kg]

%RCS calculation

F_chem = 449.9997; % [N]
F_RCS = 10.25; % [N]
l_engine_chem_offset = [0; 0.05; 0]; % [m]
l_engine_chem_missalignment = [0; 0; 1]; % [deg]
l_engine_chem_missalignment = l_engine_chem_missalignment * pi / 180; % [rad]


%solar panel
A_solar = 13; % [m^2]
m_solar = A_solar/2 * 5.67; % [kg]

%INS
m_INS = 50; % [kg]

%COM
m_COM = 25; % [kg]

%CPROP
m_engine_chem = 5.46; % [kg]
m_CPROP_dry = 29.295; % [kg]
m_MMH_tank = 22.05; % [kg]
m_MON_tank = 22.05; % [kg]
m_He_tank = 8.925; % [kg]
l_engine_chem = [0; 0; 0]; % [m]
l_CPROP_dry = [0.5763 ; 0 ; 0]; % [m]
l_MMH_tank_1 = [1.274 ; 0 ; 0]; % [m]
l_MMH_tank_2 = [0 ; 0 ; 0]; % [m]
l_MON_tank_1 = [1.274 ; 0.8 ; 0]; % [m]
l_MON_tank_2 = [1.274 ; 0.8 ; 0]; % [m]
l_He_tank_1 = [2.347 ; 0.4+0.04 ; 0]; % [m]
l_He_tank_2 = [2.347 ; -0.4-0.04 ; 0]; % [m]

geo_MMH_tank = [0.753; 0.928; 0.331]; % [m, m, m^3] diameter and length and volume of tank
geo_MON_tank = [0.753; 0.631; 0.198]; % [m, m, m^3] diameter and length and volume of tank
geo_He_tank = [0.432; 0.467; 0.04]; % [m, m, m^3] diameter and length and volume of tank
geo_Xe_tank = [0.405; 1.072; 0.120]; % [m, m, m^3] diameter and length and volume of tank


m_CPROP = m_engine_chem + m_CPROP_dry + m_MMH_tank + 2 * m_MON_tank + 2 * m_He_tank; % [kg]

%EPROP
m_engige_el = 5.04; % [kg]
m_EPROP_dry = 20.16; % [kg]
m_Xe_tank = 20.16; % [kg]
l_engine_el = [4 ; 0 ; 0]; % [m]
l_EPROP_dry = [3.5 ; 0 ; 0]; % [m]
l_Xe_tank_1 = [2.5825 ; 0 ; 0]; % [m]
l_Xe_tank_2 = [0 ; 0 ; 0]; % [m]

m_EPROP = m_engige_el + m_EPROP_dry + m_Xe_tank; % [kg]

%Propellant
OF = 1.65; % [1]
m_cprop_tot = 751.275; % [kg]
m_cprop_maneuver = [309.2921; 0.7675; 411.5063; 6.3687;  0.4490]; % [kg] propellant needed in the chemical maneuvers
m_eprop_tot = 166; % [kg]
m_eprop_maneuver = [51.805485270; 56.0652395514; 6.9816; 47.72525170401]; % [kg] propellant needed in the electrical maneuvers
m_MMH = 1/(OF+1) * [m_cprop_tot; m_cprop_tot-m_cprop_maneuver(1); m_cprop_tot-m_cprop_maneuver(1)-m_cprop_maneuver(2); 
                    m_cprop_tot-m_cprop_maneuver(1)-m_cprop_maneuver(2)-m_cprop_maneuver(3); 
                    m_cprop_tot-m_cprop_maneuver(1)-m_cprop_maneuver(2)-m_cprop_maneuver(3)-m_cprop_maneuver(4); 
                    m_cprop_tot-m_cprop_maneuver(1)-m_cprop_maneuver(2)-m_cprop_maneuver(3)-m_cprop_maneuver(4)-m_cprop_maneuver(5); 
                    m_cprop_tot-m_cprop_maneuver(1)-m_cprop_maneuver(2)-m_cprop_maneuver(3)-m_cprop_maneuver(4)-m_cprop_maneuver(5); 
                    m_cprop_tot-m_cprop_maneuver(1)-m_cprop_maneuver(2)-m_cprop_maneuver(3)-m_cprop_maneuver(4)-m_cprop_maneuver(5); 
                    m_cprop_tot-m_cprop_maneuver(1)-m_cprop_maneuver(2)-m_cprop_maneuver(3)-m_cprop_maneuver(4)-m_cprop_maneuver(5); 
                    m_cprop_tot-m_cprop_maneuver(1)-m_cprop_maneuver(2)-m_cprop_maneuver(3)-m_cprop_maneuver(4)-m_cprop_maneuver(5)]; % [kg]
m_MON = OF * m_MMH; % [kg]
m_Xe = [m_eprop_tot; m_eprop_tot; m_eprop_tot; m_eprop_tot; m_eprop_tot; m_eprop_tot; m_eprop_tot-m_eprop_maneuver(1); 
        m_eprop_tot-m_eprop_maneuver(1)-m_eprop_maneuver(2); m_eprop_tot-m_eprop_maneuver(1)-m_eprop_maneuver(2)-m_eprop_maneuver(3); 
        m_eprop_tot-m_eprop_maneuver(1)-m_eprop_maneuver(2)-m_eprop_maneuver(3)-m_eprop_maneuver(4)]; % [kg]


ratio_MMH = m_MMH ./ m_MMH(1); % [1]
ratio_MON = m_MON ./ m_MON(1); % [1]

m_He_maneuver_MMH = [0.028757; 0.41319; 0.0010257; 0.53192; 0.0078748; 0.00058064]; % [kg] first value is the ullage mass of Helium already in the tank befor the burn
m_He_maneuver_MON = [0.22998; 0.42268; 0.0010491; 0.54376; 0.0086394; 0.00058064]; % [kg] first value is the ullage mass of Helium already in the tank befor the burn

m_He_tot = 3.5251; % [kg] all Helium in the system
m_He_down_MMH = [m_He_maneuver_MMH(1); m_He_maneuver_MMH(1) + m_He_maneuver_MMH(2); m_He_maneuver_MMH(1) + m_He_maneuver_MMH(2) + m_He_maneuver_MMH(3); 
                 m_He_maneuver_MMH(1) + m_He_maneuver_MMH(2) + m_He_maneuver_MMH(3) + m_He_maneuver_MMH(4); 
                 m_He_maneuver_MMH(1) + m_He_maneuver_MMH(2) + m_He_maneuver_MMH(3) + m_He_maneuver_MMH(4) + m_He_maneuver_MMH(5);
                 m_He_maneuver_MMH(1) + m_He_maneuver_MMH(2) + m_He_maneuver_MMH(3) + m_He_maneuver_MMH(4) + m_He_maneuver_MMH(5) + m_He_maneuver_MMH(6);
                 m_He_maneuver_MMH(1) + m_He_maneuver_MMH(2) + m_He_maneuver_MMH(3) + m_He_maneuver_MMH(4) + m_He_maneuver_MMH(5) + m_He_maneuver_MMH(6);
                 m_He_maneuver_MMH(1) + m_He_maneuver_MMH(2) + m_He_maneuver_MMH(3) + m_He_maneuver_MMH(4) + m_He_maneuver_MMH(5) + m_He_maneuver_MMH(6);
                 m_He_maneuver_MMH(1) + m_He_maneuver_MMH(2) + m_He_maneuver_MMH(3) + m_He_maneuver_MMH(4) + m_He_maneuver_MMH(5) + m_He_maneuver_MMH(6);
                 m_He_maneuver_MMH(1) + m_He_maneuver_MMH(2) + m_He_maneuver_MMH(3) + m_He_maneuver_MMH(4) + m_He_maneuver_MMH(5) + m_He_maneuver_MMH(6)]; % [kg] Helium in the lower tank for MMH pessurization

m_He_down_MON = [m_He_maneuver_MON(1); m_He_maneuver_MON(1) + m_He_maneuver_MON(2); m_He_maneuver_MON(1) + m_He_maneuver_MON(2) + m_He_maneuver_MON(3);
                 m_He_maneuver_MON(1) + m_He_maneuver_MON(2) + m_He_maneuver_MON(3) + m_He_maneuver_MON(4);
                 m_He_maneuver_MON(1) + m_He_maneuver_MON(2) + m_He_maneuver_MON(3) + m_He_maneuver_MON(4) + m_He_maneuver_MON(5);
                 m_He_maneuver_MON(1) + m_He_maneuver_MON(2) + m_He_maneuver_MON(3) + m_He_maneuver_MON(4) + m_He_maneuver_MON(5) + m_He_maneuver_MON(6);
                 m_He_maneuver_MON(1) + m_He_maneuver_MON(2) + m_He_maneuver_MON(3) + m_He_maneuver_MON(4) + m_He_maneuver_MON(5) + m_He_maneuver_MON(6);
                 m_He_maneuver_MON(1) + m_He_maneuver_MON(2) + m_He_maneuver_MON(3) + m_He_maneuver_MON(4) + m_He_maneuver_MON(5) + m_He_maneuver_MON(6);
                 m_He_maneuver_MON(1) + m_He_maneuver_MON(2) + m_He_maneuver_MON(3) + m_He_maneuver_MON(4) + m_He_maneuver_MON(5) + m_He_maneuver_MON(6);
                 m_He_maneuver_MON(1) + m_He_maneuver_MON(2) + m_He_maneuver_MON(3) + m_He_maneuver_MON(4) + m_He_maneuver_MON(5) + m_He_maneuver_MON(6)]; % [kg] Helium in the lower tank for MON pessurization

m_He_up = m_He_tot - m_He_down_MMH - m_He_down_MON; % [kg]

m_He_tank_1 = m_He_up/2; % [kg]
m_He_tank_2 = m_He_up/2; % [kg]
m_He_MMH_tank_1 = m_He_down_MMH; % [kg]
m_He_MMH_tank_2 = 0; % [kg]
m_He_MON_tank_1 = m_He_down_MON/2; % [kg]
m_He_MON_tank_2 = m_He_down_MON/2; % [kg]
m_MMH_tank_1 = m_MMH; % [kg]
m_MMH_tank_2 = 0; % [kg]
m_MON_tank_1 = m_MON/2; % [kg]
m_MON_tank_2 = m_MON/2; % [kg]
m_Xe_tank_1 = m_Xe; % [kg]
m_Xe_tank_2 = 0; % [kg]

m_cprop = m_MMH + m_MON + m_He_tot; % [kg]
m_eprop = m_Xe; % [kg]


xl_MMH = l_MMH_tank_1(1) - geo_MMH_tank(2) /2 + geo_MMH_tank(2)/2 .* ratio_MMH; % [m] absolute position of MMH in tank 1 and 2 in x
xl_MON = l_MON_tank_1(1) - geo_MON_tank(2) /2 + geo_MON_tank(2)/2 .* ratio_MON; % [m] absolute position of MMH in tank 1 and 2 in x


l_He_up_tank_1 = l_He_tank_1; % [m]
l_He_up_tank_2 = l_He_tank_2; % [m]

xl_He_down_MMH = l_MMH_tank_1(1) + geo_MMH_tank(2) /2 - geo_MMH_tank(2)/2 .* (1-ratio_MMH); % [m]
xl_He_down_MON = l_MON_tank_1(1) + geo_MON_tank(2) /2 - geo_MON_tank(2)/2 .* (1-ratio_MON); % [m]

l_Xe_1 = l_Xe_tank_1; % [m]
l_Xe_2 = 0; % [m]


m_tot = m_CPROP + m_EPROP + m_cprop + m_eprop; % [kg]
m_inertia = m_dry - m_CPROP - m_EPROP  - 2 * m_solar - m_INS - m_COM; % [kg] for later inertia calculation, only the inner block for the cuboid
V_inertia = m_inertia / 1330; % [m^3] density of a 1U cubesat



%center of gravity calculation for all maneuvers for the prop system

cgx_pre = (m_engine_chem * l_engine_chem(1) + m_CPROP_dry * l_CPROP_dry(1) +  m_MMH_tank * l_MMH_tank_1(1) ...
       + m_MON_tank * l_MON_tank_1(1) + m_MON_tank * l_MON_tank_2(1) + m_He_tank * l_He_tank_1(1) ...
       + m_He_tank * l_He_tank_2(1) + m_engige_el * l_engine_el(1) + m_EPROP_dry * l_EPROP_dry(1) ...
       + m_Xe_tank * l_Xe_tank_1(1) + m_He_tank_1 .* l_He_up_tank_1(1) ...
       + m_He_tank_2 .* l_He_up_tank_2(1) + m_He_MMH_tank_1 .* xl_He_down_MMH + m_He_MMH_tank_2 .* xl_He_down_MMH ...
       + m_He_MON_tank_1 .* xl_He_down_MON + m_He_MON_tank_2 .* xl_He_down_MON + m_MMH_tank_1 .* xl_MMH + m_MMH_tank_2 .* xl_MMH ...
       + m_MON_tank_1 .* xl_MON + m_MON_tank_2 .* xl_MON + m_Xe_tank_1 .* l_Xe_1(1) + m_Xe_tank_2 .* l_Xe_2(1) ) ...
       ./ m_tot; % [m]

%position optimized for station keeping


l_solar_1 = [cgx_pre(8); 4.6029; 0]; % [m]
l_solar_2 = [cgx_pre(8); -l_solar_1(2); 0]; % [m]
l_INS = [cgx_pre(8); 0; 0.83]; % [m]
l_COM = [cgx_pre(8); 0; -l_INS(3) * m_INS/m_COM]; % [m]
l_inertia = [cgx_pre(8); 0; 0]; % [m]

cgx = (cgx_pre .* m_tot + m_solar * l_solar_1(1) + m_solar * l_solar_1(1) +... %cg with all systems involved
        m_COM * l_COM(1) + m_INS * l_INS(1) + m_inertia * l_inertia(1)) ...
       ./(m_dry + m_cprop + m_eprop); % [m]

cg = zeros(3,10); % [m] [cgx0, cgx1, cgx2, ...; cgy0, cgy1, cgy2, ...; cgz0, cgz1, cgz2, ...]
cg(1,:) = cgx; % [m]
%cg(2,:) = 0.06; % [m] for sensitivity analysis
%cg(3,:) = 0.06; % [m] for sensitivity analysis

l_MMH_1 = zeros(3,10); % [m]
l_MMH_1(1,:) = xl_MMH; % [m]
l_MMH_1(2,:) = l_MMH_tank_1(2); % [m]
l_MMH_1(3,:) = l_MMH_tank_1(3); % [m]

l_MMH_2 = zeros(3,10); % [m]
l_MMH_2(1,:) = 0; % [m]
l_MMH_2(2,:) = 0; % [m]
l_MMH_2(3,:) = 0; % [m]


l_MON_1 = zeros(3,10); % [m]
l_MON_1(1,:) = xl_MON; % [m]
l_MON_1(2,:) = l_MON_tank_1(2); % [m]
l_MON_1(3,:) = l_MON_tank_1(3); % [m]

l_MON_2 = zeros(3,10); % [m]
l_MON_2(1,:) = xl_MON; % [m]
l_MON_2(2,:) = l_MON_tank_2(2); % [m]
l_MON_2(3,:) = l_MON_tank_2(3); % [m]

l_He_down_MMH_1 = zeros(3,10); % [m]
l_He_down_MMH_1(1,:) = xl_He_down_MMH; % [m]
l_He_down_MMH_1(2,:) = l_MMH_tank_1(2); % [m]
l_He_down_MMH_1(3,:) = l_MMH_tank_1(3); % [m]

l_He_down_MMH_2 = zeros(3,10); % [m]
l_He_down_MMH_2(1,:) = 0; % [m]
l_He_down_MMH_2(2,:) = 0; % [m]
l_He_down_MMH_2(3,:) = 0; % [m]


l_He_down_MON_1 = zeros(3,10); % [m]
l_He_down_MON_1(1,:) = xl_He_down_MON; % [m]
l_He_down_MON_1(2,:) = l_MON_tank_1(2); % [m]
l_He_down_MON_1(3,:) = l_MON_tank_1(3); % [m]

l_He_down_MON_2 = zeros(3,10); % [m]
l_He_down_MON_2(1,:) = xl_He_down_MON; % [m]
l_He_down_MON_2(2,:) = l_MON_tank_2(2); % [m]
l_He_down_MON_2(3,:) = l_MON_tank_2(3); % [m]



%moment of inertia

m_inertia = m_inertia * ones(1,1,10); % [kg]
m_inertia(8) = m_inertia(8) - 20; % [kg] accounting for the dry mass without the probe
m_inertia(9) = m_inertia(9) - 20; % [kg] accounting for the dry mass without the probe
m_inertia(10) = m_inertia(10) - 20; % [kg] accounting for the dry mass without the probe
 


I_sat = 1/12 .* m_inertia .*  [(2*l_sat(2))^2 + (2*l_sat(3))^2; (2*l_sat(1))^2 + (2*l_sat(3))^2; (2*l_sat(1))^2 + (2*l_sat(2))^2]; % [kgm^2]
I_sat_cg = I_sat + m_inertia .* (cg - cg(:,8)).^2; % [kgm^2]

L = zeros(3,28,10); %matrix to put all parts in

for i = 1:10 %dimension 1 = data of part at certain timestep, dimension 2 = parts, dimension 3 = time
L(:,1,i) = l_engine_chem(:,1); % [m]
L(:,2,i) = l_CPROP_dry(:,1); % [m]
L(:,3,i) = l_MMH_tank_1(:,1); % [m]
L(:,4,i) = l_MMH_tank_2(:,1); % [m]
L(:,5,i) = l_MON_tank_1(:,1); % [m]
L(:,6,i) = l_MON_tank_2(:,1); % [m]
L(:,7,i) = l_He_tank_1(:,1); % [m]
L(:,8,i) = l_He_tank_2(:,1); % [m]
L(:,9,i) = l_engine_el(:,1); % [m]
L(:,10,i) = l_EPROP_dry(:,1); % [m]
L(:,11,i) = l_Xe_tank_1(:,1); % [m]
L(:,12,i) = l_Xe_tank_2(:,1); % [m]
L(:,13,i) = l_He_up_tank_1(:,1); % [m]
L(:,14,i) = l_He_up_tank_2(:,1); % [m]
L(:,15,i) = l_He_down_MMH_1(:,i); % [m]
L(:,16,i) = l_He_down_MMH_2(:,i); % [m]
L(:,17,i) = l_He_down_MON_1(:,i); % [m]
L(:,18,i) = l_He_down_MON_2(:,i); % [m]
L(:,19,i) = l_MMH_1(:,i); % [m]
L(:,20,i) = l_MMH_2(:,i); % [m]
L(:,21,i) = l_MON_1(:,i); % [m]
L(:,22,i) = l_MON_2(:,i); % [m]
L(:,23,i) = l_Xe_1(:,1); % [m]
L(:,24,i) = l_Xe_2(:,1); % [m]
L(:,25,i) = l_solar_1(:,1); % [m]
L(:,26,i) = l_solar_2(:,1); % [m]
L(:,27,i) = l_INS(:,1); % [m]
L(:,28,i) = l_COM(:,1); % [m]   
end

cg_inertia = zeros(3,1,10); %transformed cg vector to matrix for inertia calculation
for i = 1:10
    cg_inertia(:,1,i) = cg(:,i); % [m]
end

L_1 = zeros(3,28,10); %relative coordinate calculation

for j = 1:28
    for i = 1:10
        L_1(:,j,i) = L(:,j,i)-cg_inertia(:,1,i); % [m]
    end
end


L_inertia = zeros(3,28,10); %squared and added values for inertia calculation
for j = 1:10
    for i = 1:28   
        L_inertia(1,i,j) = L_1(2,i,j)^2 + L_1(3,i,j)^2; % [m^2]
        L_inertia(2,i,j) = L_1(1,i,j)^2 + L_1(3,i,j)^2; % [m^2]
        L_inertia(3,i,j) = L_1(1,i,j)^2 + L_1(2,i,j)^2; % [m^2]
        
    end
end

L_inertia(:,4,:) = 0; % [m]
L_inertia(:,16,:) = 0; % [m]
L_inertia(:,20,:) = 0; % [m]
L_inertia(:,12,:) = 0; % [m]
L_inertia(:,24,:) = 0; % [m]


I_tot = I_sat + m_engine_chem .* L_inertia(:,1,:) + m_CPROP_dry .* L_inertia(:,2,:) + m_MMH_tank .* L_inertia(:,3,:) ...
        + m_MMH_tank .* L_inertia(:,4,:) + m_MON_tank .* L_inertia(:,5,:) + m_MON_tank .* L_inertia(:,6,:) ...
        + m_He_tank .* L_inertia(:,7,:) + m_He_tank .* L_inertia(:,8,:) + m_engige_el .* L_inertia(:,9,:) ...
        + m_EPROP_dry .* L_inertia(:,10,:).^2 + m_Xe_tank .* L_inertia(:,11,:) + m_Xe_tank .* L_inertia(:,12,:) ...
        + m_He_tank_1' .* L_inertia(:,13,:) + m_He_tank_2' .* L_inertia(:,14,:) ...
        + m_He_MMH_tank_1' .* L_inertia(:,15,:) + m_He_MMH_tank_2' .* L_inertia(:,16,:) ...
        + m_He_MON_tank_1' .* L_inertia(:,17,:) + m_He_MON_tank_2' .* L_inertia(:,18,:) ...
        + m_MMH_tank_1' .* L_inertia(:,19,:) + m_MMH_tank_2' .* L_inertia(:,20,:) ...
        + m_MON_tank_1' .* L_inertia(:,21,:) + m_MON_tank_2' .* L_inertia(:,22,:) ... 
        + m_Xe_tank_1' .* L_inertia(:,23,:) + m_Xe_tank_2' .* L_inertia(:,24,:) ...
        + m_solar .* L_inertia(:,25,:) + m_solar .* L_inertia(:,26,:) ...
        + m_INS .* L_inertia(:,27,:) + m_COM .* L_inertia(:,28,:); % [kgm^2]


Inertia_cond = zeros(1,3);

    if I_tot(1,:,1) < I_tot(2,:,1) 
        Inertia_cond(1) = 1;
    else 
        Inertia_cond(1) = 0;
    end

    if I_tot(2,:,1) < I_tot(3,:,1) 
        Inertia_cond(2) = 1;
    else 
        Inertia_cond(2) = 0;
    end

    if I_tot(1,:,1) < I_tot(3,:,1) 
        Inertia_cond(3) = 1;
    else 
        Inertia_cond(3) = 0;
    end


%RCS
off = 0.5596; %offset for RCS calculation or additional distance to the cuboids corners in Y and Z
l_RCS_1 = [cgx_pre(8); l_sat(2)+off; l_sat(3)+off]; % [m]
l_RCS_2 = [cgx_pre(8); l_sat(2)+off; l_sat(3)+off]; % [m]

T_offset = F_chem * l_engine_chem_offset; % [Nm]
l_engine_chem_missalignment(1) = sin(l_engine_chem_missalignment(1)); % [1]
l_engine_chem_missalignment(2) = sin(l_engine_chem_missalignment(2)); % [1]
l_engine_chem_missalignment(3) = sin(l_engine_chem_missalignment(3)); % [1]

a = F_chem * l_engine_chem_missalignment; % [N]
A = repmat(a,1,10); % [N]

b = (cg - l_engine_chem); % [m]
T_missalignment = cross(A,b); % [Nm]

T_tot = (T_offset+T_missalignment); % [Nm]
T_tot = abs(T_tot); % [Nm]
T_RCS = abs(F_RCS * l_RCS_1) + abs(F_RCS * l_RCS_2); % [Nm]
T_RCS_matrix = repmat(T_RCS,1,10); % [Nm]

T_tot = T_tot(:,1:6); % [Nm]
T_RCS_matrix = T_RCS_matrix(:,1:6); % [Nm]


RCS_cond = zeros(1,3);

for i = 1:3
    if T_RCS_matrix(i,:) >= T_tot(i,:)
        RCS_cond(i) = 1;
    else 
        RCS_cond(i) = 0;
    end
end


%output

M = zeros(3,8);
M(:,1) = l_MMH_tank_1; % [m]
M(:,2) = l_MMH_tank_2; % [m]
M(:,3) = l_MON_tank_1; % [m]
M(:,4) = l_MON_tank_2; % [m]
M(:,5) = l_He_tank_1; % [m]
M(:,6) = l_He_tank_2; % [m]
M(:,7) = l_Xe_tank_1; % [m]
M(:,8) = l_Xe_tank_2; % [m]


disp("Center of gravity in [x;y;z] over 10 time steps [0,1,2,...] [m]");
disp(cg);
disp("Total moment of inertia around the main axis in [x;y;z] over 10 timesteps [0,1,2,...] [kgm^2]");
disp(I_tot(:,:,1));
disp("Geometrical location of parts in [x;y;z] over the MMH tank 1 % 2, MON tank 1 % 2, He tank 1 % 2 and Xe tank 1 % 2 [m]");
disp(M);
disp("Inertia analysis Iz > Iy > Ix over all timesteps with first Ix < Iy, second Iy < Iz, thrid Ix < Iz");

for i = 1:3
    if Inertia_cond(i) == 1
        disp(["Inertia concition in", i, "matched", "1 = X, 2 = Y, 3 = Z"]);
    else 
        disp(["Inertia concition in", i, "dismatched", "1 = X, 2 = Y, 3 = Z"]);
    end
end

disp("Acting torques that need to be coutered by RCS in [x;y;z] over 6 timesteps [0,1,2,...] [Nm]");
disp(T_tot);
disp("Torques RCS can provide in [x;y;z] over 6 timesteps [0,1,2,...] [Nm]");
disp(T_RCS_matrix);
disp("RCS analysis");

for i = 1:3
    if RCS_cond(i) == 1
        disp(["RCS concition in", i, "matched", "1 = X, 2 = Y, 3 = Z"]);
    else 
        disp(["RCS concition in", i, "dismatched", "1 = X, 2 = Y, 3 = Z"]);
    end
end

disp("Volume analysis");
    if V_sat >= V_inertia
        disp("SC is big enough for the mass to fit within it");
    else 
        disp("SC is not big enough for the mass to fit within it");
    end
    
    if l_sat(2) *2 <= 4.6
        disp("SC is small enough for the faring diameter");
    else 
        disp("SC is not small enough the faring diameter");
    end

disp("Solar panel dimensions [m]");
disp(["length", 2* (l_solar_1(2)-l_sat(2)), "width", A_solar/(2* 2* (l_solar_1(2)-l_sat(2)))]);