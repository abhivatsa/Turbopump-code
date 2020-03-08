clc;
close all;
clear all;

%% ------------ shaft discretization for turbine impeller and inducer ----------------

do=0.001*[8 8 8 8 10 10 10 10 10 11 12 11 10.5 11 11 12 12 12 14 14 18 18 18 18 14 14 15.5 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 16 14 11 16 17 17]; %outer dia. of shaft
di=0.001*[3.50 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; % inner dia. of shaft
le=0.001*[2.6 2.4 2 3 3 3 3 3 3 3 3 3.5 2.5 3 3 3 3 2 3 2 3 3 3 2 2 3 1.50 3 3 3 3 3 3 3 3 3 3 3 3 3 2 2 3 3 3 2 2 4.5 3 4 4 4.7]; % length of each sub-element

E=211*10^9; % elastic modulus of steel
rho=7810; % density of steel
% v=0.3; % poisson ratio of steel
G=81.2*10^9; % shear modulus of steel

% G=E/(2*(1+v));
v=E/(2*G)-1;

tot_num_ele=length(le);
tot_node=tot_num_ele+1;

[M_glo_mat_shaft_bend,M_glo_mat_shaft_rotary,K_glo_mat_shaft_bend,G_glo_mat_shaft,K_glo_mat_skew] =shaft_glo_matrix(do,di,le,E,rho,v,G,tot_node);


%% ------------------- Disc Properties ----------------------

disc_node=[4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 53]; 
md=[0.052/20 0.052/20 0.052/20 0.052/20 0.052/20 0.052/20 0.052/20 0.052/20 0.052/20 0.052/20 (0.052/20+0.115/10) (0.052/20+0.115/10) (0.052/20+0.115/10) (0.052/20+0.115/10) (0.052/20+0.115/10) (0.052/20+0.115/10) (0.052/20+0.115/10) (0.052/20+0.115/10) (0.052/20+0.115/10) (0.052/20+0.115/10) 0.005 0.605];
Ip=1e-6* [0 0 0 0 0 0 0 0 0 3.883 0 0 0 0 0 29.737 0 0 0 0 0 779.45];
Id=[0 0 0 0 0 0 0 0 0 1.645e-5 0 0 0 0 0 2.127e-5 0 0 0 0 0 3.99e-4];

[M_glo_mat_disc,G_glo_mat_disc]=Disc_glo_matrix(tot_node,disc_node,md,Ip,Id);

%% ------------------ Bearing properties ------------------

bearing_node=[29 43];

%% ------------------ internal Damping coeff ----------------

nh=0;
nv=0.00001;

%% --------------------- Seal properties ---------------------

% delp=0.2*[0; 22.813; 27.983; 37.923; 148.9]*10^5;
% speed=[0; 25998; 28932; 33324; 61300];
% 
% p_seal=polyfit(speed,delp,3);
% 
rho_fluid=750; % density of fluid
% mu_fluid=(0.85e-6)*750; % dynamic viscosity of fluid
% H= 0.025; % seal radial clearence
% Clearence=5*10^(-4); % nominal seal radial clearence
% len_seal=9*10^(-3); % seal length
% v0=0; %initial swirl
% D=24*10^(-3); % diameter
% zeta=0.1; %inlet poressure loss coefficient
% seal_node=26;

K_seal_speed=[];
C_seal_speed=[];
M_seal_speed=[];

%% ------------------ Alford force torque ------------------

% T_data=[0; 3.82; 4.83; 6.67; 16.77];
% speed=[0; 25998; 28932; 33324; 61300];
% 
% p=polyfit(speed,T_data,3);
% 
% speed_rpm=0:1000:100000;
% T_vec=polyval(p,speed_rpm).';
% 
% turbine_node=53;

%% ------------------ unbalance response --------------------
% [Force/Moment, node, value, phase]
unb_mass=[1 40 0.395*10^(-6) 0];

%% ------------------ hydraulic unbalance -------------------

imp_dia= 0.055;
exit_width=0.00311;
imp_node=19;
hyd_force_array=[imp_node, rho_fluid, imp_dia, exit_width];

%% ------------------ Solver --------------------
ctr=0;
for N=1000:1000:100000% speed in RPM
    
    ctr=ctr+1;
    
    w=2*pi*N/60;  % speed in rad/sec
    
%     delp=polyval(p_seal,N);
%     T=T_vec(ctr+1);

    % ------------- bearing matrix calculation -----------
    [K_glo_mat_bearing,C_glo_mat_bearing]=bearing_glo_matrix(bearing_node,tot_node,N);
    
    % ------------- hydraulic effect calculation -----------
%     [K_glo_mat_imp,C_glo_mat_imp,M_glo_mat_imp]=impeller_diffuser(N,rho_fluid,imp_dia,exit_width,imp_node,tot_node);
    
    % ------------- alford force effect calculation -----------
%    [K_alford]=alford_force(turbine_node,tot_node,T);
    
    % ------------- seal effect calculation -----------
%    [K_glo_mat_seal,C_glo_mat_seal,M_glo_mat_seal,K_seal,k_seal,C_seal,c_seal,m_seal]=seal(N,rho_fluid,mu_fluid,Clearence,len_seal,v0,D,zeta,delp,seal_node,tot_node);

    % -------------- Assembly of Global Matrices ---------------
    M_glo_mat= M_glo_mat_shaft_bend + M_glo_mat_disc + M_glo_mat_shaft_rotary;%+M_glo_mat_seal+M_glo_mat_imp;
    K_glo_mat= (1+nh)/(sqrt(1+nh^2))*K_glo_mat_shaft_bend + K_glo_mat_bearing + (nv*w+nh/(sqrt(1+nh^2)))*K_glo_mat_skew;%+K_glo_mat_seal+K_glo_mat_imp;%+K_alford;
    G_glo_mat= G_glo_mat_shaft + G_glo_mat_disc;
    C_glo_mat= w*G_glo_mat+ C_glo_mat_bearing+nv*K_glo_mat_shaft_bend;%+C_glo_mat_seal+C_glo_mat_imp;
    
    % ------------ Constraint update in Global Matrix ------------
    % ----------- if some node is fixed -----------
    constraint_vec=[1:212];
    count_row=0;

    for x=constraint_vec
        count_row=count_row+1;
        count_column=0;
        for y=constraint_vec
            count_column=count_column+1;
            Global_M(count_row,count_column)=M_glo_mat(x,y);
            Global_K(count_row,count_column)=K_glo_mat(x,y);
            Global_C(count_row,count_column)=C_glo_mat(x,y);
        end
    end
    
    
    % ------------- State-Space formulation ---------------
    len_zero=length(Global_M);
    
    A=[zeros(len_zero) eye(len_zero); -inv(Global_M)*Global_K -inv(Global_M)*Global_C];
    
    % ----------- Computation of Eigen values & Eigen Vectors ---------
    [eig_vec_right,eig_val,eig_vec_left]=eig(A);
    eig_val1=diag(eig_val);

    % ------------ Sorting of eigen values and vectors ------------
    [eig_val2,eig_val_index]=sort(eig_val1);
    
    for j=1:length(eig_val2)
        eig_vec1(:,j)=eig_vec_right(:,eig_val_index(j));
    end
    
    % ------------- Imaginary and real part --------------
    eig_val_imag=imag(eig_val2);
    eig_val_real=real(eig_val2);
    
    eig_vec2=eig_vec1(1:len_zero,:);
    
    %     Eigen values, Eigen vectors & Speed
    nat_freq_update(:,ctr)=eig_val_imag;
    nat_freq=nat_freq_update(1:12,:);
    eig_vec_final(:,:,ctr)=eig_vec2;
    freq_rpm(ctr,1)=N;
    
    % ----------- Stability_analysis -----------
    Stab_anal(ctr,:)=eig_val_real(1:12);
    Stab_anal_seal(ctr)=max(eig_val_real);
    

% K_seal_speed(:,:,ctr)=K_glo_mat_seal;
% C_seal_speed(:,:,ctr)=C_glo_mat_seal;
% M_seal_speed(:,:,ctr)=M_glo_mat_seal;

end
nf=nat_freq_update/(2*pi);

for i=1:length(Stab_anal_seal)
    if Stab_anal_seal(i)<0
        unstab_speed=freq_rpm(i);
    end
end

L=[0 cumsum(le)];

for i=2:length(freq_rpm)
    for j=1:length(nat_freq(:,1))
        u=eig_vec_final(1:4:end,j,i);
        v=eig_vec_final(2:4:end,j,i);
        kappa(:,j,i)=whirl(u,v);
    end
end

% ----- calculation of steady state response due to various forces ------
response=force_response(K_glo_mat_shaft_bend,K_glo_mat_skew,G_glo_mat,(M_glo_mat_shaft_bend + M_glo_mat_disc + M_glo_mat_shaft_rotary),K_seal_speed, C_seal_speed, M_seal_speed,[100:100:freq_rpm(end,1)],unb_mass,nv,nh,hyd_force_array,bearing_node,tot_node);
 
%% ------------------------ Campbell diagram -----------------------------
figure(1)
campbell(nat_freq,freq_rpm,kappa)

% N_Hz=(freq_rpm/60)*[4];
% plot(freq_rpm,N_Hz,'b--');
% campbell(nat_freq,freq_rpm)

%% ------------------------ Stability analysis map --------------------------
figure(2)
plot(freq_rpm,zeros(length(freq_rpm),1),'+','DisplayName','Zero line')
hold on
plot(freq_rpm,Stab_anal(:,2),'x','DisplayName','Mode1')
hold on
plot(freq_rpm,Stab_anal(:,4),'v','DisplayName','Mode2')
hold on
plot(freq_rpm,Stab_anal(:,6),'*','DisplayName','Mode3')
hold on
plot(freq_rpm,Stab_anal(:,8),'+','DisplayName','Mode4')
% hold on
% plot(freq_rpm,Stab_anal(:,10),'x','DisplayName','Mode5')
hold off
legend('show')
xlabel('Speed in RPM')
ylabel('Real part of Eigenvalues')
title('Stability analysis')
set(gca,'fontsize',24)

%% ------------------------ Mode shapes -----------------------
 
L=[0 cumsum(le)];

figure(3)
subplot(2,3,1)
mode_shapes(eig_vec_final(:,2,10),L,nat_freq_update(2,10))
subplot(2,3,2)
mode_shapes(eig_vec_final(:,4,10),L,nat_freq_update(4,10)) 
subplot(2,3,3)
mode_shapes(eig_vec_final(:,6,10),L,nat_freq_update(6,10))
subplot(2,3,4)
mode_shapes(eig_vec_final(:,8,10),L,nat_freq_update(8,10))
subplot(2,3,5)
mode_shapes(eig_vec_final(:,10,10),L,nat_freq_update(10,10))
subplot(2,3,6)
mode_shapes(eig_vec_final(:,12,10),L,nat_freq_update(12,10))

%% ---------------------- orbit plot ------------------------- 
figure(4)
subplot(2,3,1)
orbit_plot(eig_vec_final(:,2,10),[19 53],'Mode 1',nat_freq_update(2,10))
subplot(2,3,2)
orbit_plot(eig_vec_final(:,4,10),[19 53],'Mode 2',nat_freq_update(4,10))
subplot(2,3,3)
orbit_plot(eig_vec_final(:,6,10),[19 53],'Mode 3',nat_freq_update(6,10))
subplot(2,3,4)
orbit_plot(eig_vec_final(:,8,10),[19 53],'Mode 4',nat_freq_update(8,10))
subplot(2,3,5)
orbit_plot(eig_vec_final(:,10,10),[19 53],'Mode 5',nat_freq_update(10,10))
subplot(2,3,6)
orbit_plot(eig_vec_final(:,12,10),[19 53],'Mode 6',nat_freq_update(12,10))

%% ---------------------------- Unbalance response  --------------------------

figure(5)
plot_FRF([100:100:freq_rpm(end,1)],response,[19 1; 53 1],unstab_speed)
figure(6)
plot_FRF_with_backward_whirl([100:100:freq_rpm(end,1)],response,53)


%% ---------------------------  Coast-up response --------------------------
[resp,N_rpm_coast,Time]=coastup(K_glo_mat_shaft_bend,K_glo_mat_skew,G_glo_mat,(M_glo_mat_shaft_bend + M_glo_mat_disc + M_glo_mat_shaft_rotary),K_seal_speed, C_seal_speed, M_seal_speed,unb_mass,nv,nh,hyd_force_array,bearing_node,tot_node,[0 0 100],25000);
figure(7)
plot(N_rpm_coast,resp(:,210).')
xlabel('Speed (in RPM)')
ylabel('Magnitude (in y-direction)')
xlim([0 30000])
set(gca,'FontSize',24)

%% --------------------- Waterfall diagram ---------------------
[freq_wf,speed_wf,recep_wf] = waterfall_diag(resp(:,210),N_rpm_coast,Time);
figure(8)
waterfall(freq_wf,speed_wf,recep_wf(:,:))
xlabel('Frequency (Hz)')
ylabel('Speed (RPM)')
zlabel('Amplitude (m)')
xlim([0 400])
ylim([0 24000])
set(gca,'fontsize',24)


% %% Bent response
% 
% % (K_glo_mat_shaft_bend + K_glo_mat_bearing),G_glo_mat,C_glo_mat_bearing,(M_glo_mat_shaft_bend + M_glo_mat_disc + M_glo_mat_shaft_rotary),K_seal_speed, C_seal_speed, M_seal_speed
% 
% % bend_meas_node=[3 11 24 36 53];
% % bend_defl=0.001*[0.156*(exp(1i*237/180*pi)) 0.029*(exp(1i*14.91/180*pi)) 0.075*(exp(1i*358.33/180*pi)) 0.023*(exp(1i*125.04/180*pi)) 0.114*(exp(1i*125.04/180*pi))];
% % response_bend=bend_response(K_glo_mat_shaft_bend,K_glo_mat_bearing,K_glo_mat_skew,G_glo_mat,C_glo_mat_bearing,(M_glo_mat_shaft_bend + M_glo_mat_disc + M_glo_mat_shaft_rotary),[10:10:freq_rpm(end,1)],bend_meas_node.',bend_defl.',nv,nh);
% % figure(8)
% % plot_FRF([10:10:freq_rpm(end,1)],response_bend,[19 1; 53 1],unstab_speed)
% % figure(9)
% % plot_FRF_with_backward_whirl([10:10:freq_rpm(end,1)],response_bend,53)

% % %% operating deflection shape
% % figure(10)
% % operating_deflection_shapes(response(:,201),L,20000)
% % figure(11)
% % operating_deflection_shapes(response_bend(:,2001),L,20000)
