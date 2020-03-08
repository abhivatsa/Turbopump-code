function [response_bend]=bend_response(K_glo_mat_shaft_bend,K_glo_mat_bearing,K_glo_mat_skew,G,Cb,M,N,bend_meas_node,bend_defl,nv,nh)
w=2*pi*N/60;
resp_master=[];
master_dof=[];
for count=1:length(bend_meas_node)
    master_dof=[master_dof 4*bend_meas_node(count)-3 4*bend_meas_node(count)-2];
    resp_master=[resp_master; bend_defl(count,1)*[1;-1i]];
end
slave_dof=1:length(K_glo_mat_shaft_bend);
slave_dof(master_dof)=[];
Kss=K_glo_mat_shaft_bend(slave_dof,slave_dof);
Ksm=K_glo_mat_shaft_bend(slave_dof,master_dof);
resp_slave=-Kss\(Ksm*resp_master);
resp_bend(master_dof)=resp_master;
resp_bend(slave_dof)=resp_slave;
bend_force=K_glo_mat_shaft_bend*resp_bend.';

response_bend=zeros(length(M),length(w));

for j=1:length(w)
    K= (1+nh)/(sqrt(1+nh^2))*K_glo_mat_shaft_bend + K_glo_mat_bearing + (nv*w(j)+nh/(sqrt(1+nh^2)))*K_glo_mat_skew;
    C= Cb+nv*K_glo_mat_shaft_bend;
    response_bend(:,j)=((K-(w(j))^2*M)+1i*w(j)*(w(j)*G+C))\bend_force;
end
end