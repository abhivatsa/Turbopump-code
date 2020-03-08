function [resp,N_rpm,Time]=coastup(K_glo_mat_shaft_bend,K_glo_mat_skew,G,M,K_seal_speed,C_seal_speed,M_seal_speed,force_array,nv,nh,hyd_force_array,bearing_node,tot_node,alpha,speed_range)
b0=zeros(length(M),1);
hyd_force=zeros(length(M),1);
vane_pass_force=zeros(length(M),1);

for i=1:length(force_array(:,1))
    force=force_array(i,:);
    node_pos=force(2);
    if force(i,1) == 1
        b0(4*node_pos-3)=force(3)*exp(1i*force(4));
        b0(4*node_pos-2)=-1i*force(3)*exp(1i*force(4));
    end
    if force(i,1) == 2
        b0(4*node_pos-1)=1i*force(3)*exp(1i*force(4));
        b0(4*node_pos)=force(3)*exp(1i*force(4));
    end
end

phi=alpha(1);
phi_dot=alpha(2);
phi_double_dot=alpha(3);
   
N=0;
ctr=0;
ti=0;

M_new=M;
C_new=nv*K_glo_mat_shaft_bend;

while N<speed_range
    ctr=ctr+1;
    tf=ti+0.001;
    w=2*pi*N/60;
    
    [K_glo_mat_bearing,C_glo_mat_bearing]=bearing_glo_matrix(bearing_node,tot_node,N);
    K_new = (1+nh)/(sqrt(1+nh^2))*K_glo_mat_shaft_bend + K_glo_mat_bearing + (nv*w+nh/(sqrt(1+nh^2)))*K_glo_mat_skew;
    C_new = nv*K_glo_mat_shaft_bend + C_glo_mat_bearing;
    
%     if N>0
%         N_round_rpm=round(N,-2);
%         count_round=round(N_round_rpm/100);
%         if count_round>0
%             k_seal_2d(:,:)=K_seal_speed(:,:,count_round);
%             c_seal_2d(:,:)=C_seal_speed(:,:,count_round);
%             m_seal_2d(:,:)=M_seal_speed(:,:,count_round);
%             K_new=(1+nh)/(sqrt(1+nh^2))*K_glo_mat_shaft_bend + K_glo_mat_bearing + (nv*w+nh/(sqrt(1+nh^2)))*K_glo_mat_skew+k_seal_2d;
%             M_new=M+m_seal_2d;
%             C_new=C_glo_mat_bearing+nv*K_glo_mat_shaft_bend+c_seal_2d;
%         end
%     end
%     
%     if isempty(hyd_force_array) == 0
%             imp_node=hyd_force_array(1);
%             rho_fluid=hyd_force_array(2);
%             imp_dia=hyd_force_array(3);
%             exit_width=hyd_force_array(4);
%             
%             head=6.1316*10^(-7)*N^2-0.0046*N+1.7320;
%             
%             Kr=0.015;
%             mag_hyd_force=Kr*rho_fluid*9.81*head*imp_dia*exit_width;
%             mag_vane_pass_force=0.01*(1.7112*10^(-7)*N^2-6.9655*10^(-4)*N-0.1744);
% 
%             hyd_force(4*imp_node-3)=mag_hyd_force;
%             hyd_force(4*imp_node-2)=-1i*mag_hyd_force;
%             vane_pass_force(4*imp_node-3)=mag_vane_pass_force;
%             vane_pass_force(4*imp_node-2)=-1i*mag_vane_pass_force;
%             
%             [K_glo_mat_imp,C_glo_mat_imp,M_glo_mat_imp]=impeller_diffuser(N,rho_fluid,imp_dia,exit_width,imp_node,tot_node);
%             K_new = K_new + K_glo_mat_imp;
%             C_new = C_new + C_glo_mat_imp;
%             M_new = M_new + M_glo_mat_imp;
%     end
    
    A1=[zeros(length(M_new)) eye(length(M_new)); -M_new\K_new -M_new\(C_new+phi_dot*G)];
    B1=[zeros(length(M_new)); M_new];
    
    [red_A,red_B,trans_mat,inv_trans_mat]=serep_modified(A1,B1);
    
    if ctr==1
        ini_cond=zeros(length(red_A),1);
    end
    
    options=odeset;
    [t,x]=ode45(@(t,x) solver_ode45(t,x,red_A,red_B,b0,hyd_force,vane_pass_force,phi,phi_dot,phi_double_dot),[ti tf],ini_cond,options);
    phi=phi+phi_dot*(tf-ti)+0.5*phi_double_dot*(tf-ti)^2;
    phi_dot=phi_dot+phi_double_dot*(tf-ti);
    ti=tf;
    ini_cond=x(end,:);
    Time(ctr)=t(end);
    resp(ctr,:)=inv_trans_mat*x(end,:).';
    N_rpm(ctr)=60*phi_dot/(2*pi);
    N=60*phi_dot/(2*pi)
end
end

