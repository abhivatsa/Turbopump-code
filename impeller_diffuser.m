function [K_glo_mat_imp,C_glo_mat_imp,M_glo_mat_imp]=impeller_diffuser(N,rho_fluid,imp_dia,exit_width,imp_node,tot_node)
K_glo_mat_imp=zeros(4*tot_node);
C_glo_mat_imp=zeros(4*tot_node);
M_glo_mat_imp=zeros(4*tot_node);

w=2*pi*N/60;

% Dimensionless coefficients
K_star = -4.2;
k_star = 5.1;
C_star = 4.6;
c_star = 13.5;
M_star = 11.0;
m_star = 4.0;

A=pi*(imp_dia/2)^2*exit_width*rho_fluid;

K = K_star*(A*w^2);
k = k_star*(A*w^2);
C = C_star*(A*w);
c = c_star*(A*w);
M = M_star*A;
m = m_star*A;

K_mat=[K k; -k K];
C_mat=[C c; -c C];
M_mat=[M m; -m M];

m=imp_node-1;
for ele_hor_pos=1:2
    for ele_ver_pos=1:2
        K_glo_mat_imp(ele_hor_pos+4*m,ele_ver_pos+4*m)=K_glo_mat_imp(ele_hor_pos+4*m,ele_ver_pos+4*m)+K_mat(ele_hor_pos,ele_ver_pos);
        C_glo_mat_imp(ele_hor_pos+4*m,ele_ver_pos+4*m)=C_glo_mat_imp(ele_hor_pos+4*m,ele_ver_pos+4*m)+C_mat(ele_hor_pos,ele_ver_pos);
        M_glo_mat_imp(ele_hor_pos+4*m,ele_ver_pos+4*m)=M_glo_mat_imp(ele_hor_pos+4*m,ele_ver_pos+4*m)+M_mat(ele_hor_pos,ele_ver_pos);
    end
end

end
