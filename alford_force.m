function [K_alford]=alford_force(turbine_node,tot_node,T)

beta=5;
D=0.092;
L=0.021;

K_alford=zeros(4*tot_node);
Ksw=beta*T/(D*L);

for i=1:length(turbine_node)
    Kyy=0;
    Kxx=0;
    Kxy=Ksw;
    Kyx=-Ksw;
    
    K=[Kxx Kxy; Kyx Kyy];
    m=turbine_node(i)-1;
    for ele_hor_pos=1:2
        for ele_ver_pos=1:2
            K_alford(ele_hor_pos+4*m,ele_ver_pos+4*m)=K_alford(ele_hor_pos+4*m,ele_ver_pos+4*m)+K(ele_hor_pos,ele_ver_pos);
        end
    end
end
end


