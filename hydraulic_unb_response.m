function [hyd_response]=hydraulic_unb_response(K,G,Cb,M,K_seal_speed,C_seal_speed,M_seal_speed,N,rho_fluid,head,imp_dia,exit_width,imp_node)
w=2*pi*N/60;
b0=zeros(length(M),1);
Kr=0.015;
hyd_force=Kr*rho_fluid*9.81*head*imp_dia*exit_width;

b0(4*imp_node-3)=hyd_force;
b0(4*imp_node-2)=-1i*hyd_force;

hyd_response=zeros(length(M),length(w));

for j=1:length(w)
    N=60*w(j)/(2*pi);
    if N>0
        N_round_rpm=round(N,-3);
        count_round=round(N_round_rpm/1000);
        if count_round>0
            k_seal(:,:)=K_seal_speed(:,:,count_round);
            c_seal(:,:)=C_seal_speed(:,:,count_round);
            m_seal(:,:)=M_seal_speed(:,:,count_round);
            K_new=K+k_seal;
            M_new=M+m_seal;
            Cb_new=Cb+c_seal;
        end
    end
    hyd_response(:,j)=((K_new-(w(j))^2*M_new)+1i*w(j)*(w(j)*G+Cb_new))\(b0);
end

end
