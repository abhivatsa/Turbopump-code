function [K_glo_mat_bearing,C_glo_mat_bearing]=bearing_glo_matrix(bearing_node,tot_node,N)

mr=1.027; % mass of rotor
d_c_b1=37.27; % Distsnce between centre of mass and Bearing-1
d_b1_b2=40; % Distsnce between Bearing-1 and Bearing-2
d_c_b2=d_b1_b2-d_c_b1; % Distsnce between centre of mass and Bearing-2
d_i_b1=25.5; % Distsnce between impeller and Bearing-1
d_i_b2=65.5; % Distsnce between impeller and Bearing-1
z=8; % no. of balls
d=0.0065; % dia. of balls
F_sta_1=mr*9.81*d_c_b2/d_b1_b2; % Static load on Bearing-1
F_sta_2=mr*9.81*d_c_b1/d_b1_b2; % Static load on Bearing-2

z=8; % no. of balls

% Fh=1.7112*10^(-7)*N^2-6.9655*10^(-4)*N-0.1744;
% F_sta_1= abs(Fh*d_i_b2/d_b1_b2 + mr*9.81*d_c_b2/d_b1_b2); % Static load on Bearing-1
% F_sta_2= abs(Fh*d_i_b1/d_b1_b2 - mr*9.81*d_c_b1/d_b1_b2); % Static load on Bearing-2


F=[F_sta_1 F_sta_2]; % Force matrix

K_glo_mat_bearing=zeros(4*tot_node);
C_glo_mat_bearing=zeros(4*tot_node);

for i=1:length(bearing_node)
    Kyy=1.3*z^(2/3)*d^(1/3)*(F(i))^(1/3)*10^7;
%     Kxx=0.46*Kyy;
    Kxx=Kyy;
    Kxy=0;
    Kyx=0;
  
    Cxx=0;
    Cyy=0;
    
    K=[Kxx Kxy; Kyx Kyy];
    C=[Cxx 0; 0 Cyy];
    m=bearing_node(i)-1;
    for ele_hor_pos=1:2
        for ele_ver_pos=1:2
            K_glo_mat_bearing(ele_hor_pos+4*m,ele_ver_pos+4*m)=K_glo_mat_bearing(ele_hor_pos+4*m,ele_ver_pos+4*m)+K(ele_hor_pos,ele_ver_pos);
            C_glo_mat_bearing(ele_hor_pos+4*m,ele_ver_pos+4*m)=C_glo_mat_bearing(ele_hor_pos+4*m,ele_ver_pos+4*m)+C(ele_hor_pos,ele_ver_pos);
        end
    end
end
end
