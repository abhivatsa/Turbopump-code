function [response]=force_response(K_glo_mat_shaft_bend,K_glo_mat_skew,G,M,K_seal_speed,C_seal_speed,M_seal_speed,N,force_array,nv,nh,hyd_force_array,bearing_node,tot_node,rho_fluid)

w=2*pi*N/60;
b0=zeros(length(M),1);
hyd_force=zeros(length(M),1);
vane_pass_force=zeros(length(M),1);
% bend_force=zeros(length(M),1);

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

% if isempty(bend_meas_node) == 0
%     [bend_force]=force_bend(K_sh,bend_meas_node,bend_defl);
% end

response=zeros(length(M),length(w));


M_new=M;
C_new=nv*K_glo_mat_shaft_bend;

for j=1:length(w)
    N=60*w(j)/(2*pi);
    [K_glo_mat_bearing,C_glo_mat_bearing]=bearing_glo_matrix(bearing_node,tot_node,N);
    K_new = (1+nh)/(sqrt(1+nh^2))*K_glo_mat_shaft_bend + K_glo_mat_bearing + (nv*w(j)+nh/(sqrt(1+nh^2)))*K_glo_mat_skew;
    C_new = nv*K_glo_mat_shaft_bend + C_glo_mat_bearing;
%     if N>0
%         N_round_rpm=round(N,-2);
%         count_round=round(N_round_rpm/100);
%         if count_round>0
%             N_round_rpm
%             count_round
%             k_seal(:,:)=K_seal_speed(:,:,count_round);
%             c_seal(:,:)=C_seal_speed(:,:,count_round);
%             m_seal(:,:)=M_seal_speed(:,:,count_round);
%             K_new=(1+nh)/(sqrt(1+nh^2))*K_glo_mat_shaft_bend + K_glo_mat_bearing + (nv*w(j)+nh/(sqrt(1+nh^2)))*K_glo_mat_skew+k_seal;
%             M_new=M+m_seal;
%             C_new=C_glo_mat_bearing+nv*K_glo_mat_shaft_bend+c_seal;
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
      
    unbalance=(w(j))^2*b0;
%     force=bend_force+unbalance+hyd_force;
%     force=unbalance+hyd_force;
    force = unbalance;
    
    response(:,j)=((K_new-(w(j))^2*M_new)+1i*w(j)*(w(j)*G+C_new))\(force);
%     response1(:,j)=((K_new-(w(j))^2*M_new)+1i*w(j)*(w(j)*G+C_new))\(force);
%     response2(:,j)=((K_new-16*(w(j))^2*M_new)+4*1i*w(j)*(w(j)*G+C_new))\(vane_pass_force);
%     response(:,j)=response1(:,j)+response2(:,j);
    
end

end
