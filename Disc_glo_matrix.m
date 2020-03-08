function [M_glo_mat_disc,G_glo_mat_disc]=Disc_glo_matrix(tot_node,disc_node,md,Ip,Id)

M_glo_mat_disc=zeros(4*tot_node);
G_glo_mat_disc=zeros(4*tot_node);

% Disc_glo_matrix(tot_node,disc_node,md,Ip,Id)
% Disc_glo_matrix(rho,tot_node,disc_node,do_disc,di_disc,thi_disc)


for i=1:length(disc_node)

%     Ad(i)=pi*((do_disc(i))^2-(di_disc(i))^2)/4; % area of ith disc 
%     md(i)=Ad(i)*rho*thi_disc(i); % mass of ith disc
%     Ip(i)=(md(i)/8)*((do_disc(i))^2+(di_disc(i))^2);
%     Id(i)=(md(i)/12)*(3*((do_disc(i))^2+(di_disc(i))^2)/4 + (thi_disc(i))^2);

    Md=[md(i) 0 0 0 ; 0 md(i) 0 0 ;0 0 Id(i) 0; 0 0 0 Id(i)];
    Gd=[0 0 0 0; 0 0 0 0; 0 0 0 Ip(i); 0 0 -Ip(i) 0];
    
    m=disc_node(i)-1;
    for row_no=1:4
        for column_no=1:4
            M_glo_mat_disc(row_no+4*m,column_no+4*m)=M_glo_mat_disc(row_no+4*m,column_no+4*m)+Md(row_no,column_no);
            G_glo_mat_disc(row_no+4*m,column_no+4*m)=G_glo_mat_disc(row_no+4*m,column_no+4*m)+Gd(row_no,column_no);
        end
    end
end
end




% disc_node=[34]; % disc node number
% discmass=[0.605];
% discpolarinertia=1e-6* [779.45];
% discdiametralinertia=[3.99e-4];

% disc_node=[9 13 33]; % disc node number
% discmass=[0.0456 0.106 0.575];
% discpolarinertia=1e-6* [3.11579  27.4023 741.845];
% discdiametralinertia=[1.4830e-5 1.96e-5 3.79e-4];

% disc_node=[9 13 33]; % disc node number
% discmass=[0.052 0.115 0.605];
% discpolarinertia=1e-6* [3.875  29.758 779.45];
% discdiametralinertia=[1.6284e-5 2.128e-5 3.99e-4];

% dia_disc=[0.28 0.35]; % dia. of rotor disc
% thi_disc=[0.07 0.07]; % thickness of rotor disc
% 
% E=211*10^9; % elastic modulus of steel
% pho=7810; % density of steel
% % v=0.27; % poisson ratio of steel
% G=81.2*10^9; % shear modulus of steel
% 
% % G=E/(2*(1+v));
% v=E/(2*G)-1;

%     Ad=pi*((dia_disc(i))^2-(do(disc_node(i)-1))^2)/4; % area of ith disc 
%     md=Ad*pho*thi_disc(i); % mass of ith disc
%     %Id=(pi*pho*thi_disc(i)/12)*(3*((dia_disc(i))^4-(do(node))^4)/16+((thi_disc(i))^2)*((dia_disc(i))^2-(do(node))^2)/4);% diameteral moment of inertia of disc1
%     %Ip=pho*pi*((dia_disc(i))^4-(do(node))^4)*thi_disc(i)/32; %polar moment of inertia of disc1I
%     Ip=(md/8)*((dia_disc(i))^2+(do(disc_node(i)-1))^2);
%     Id=(md/12)*(3*((dia_disc(i))^2+(do(disc_node(i)-1))^2)/4 + (thi_disc(i))^2);