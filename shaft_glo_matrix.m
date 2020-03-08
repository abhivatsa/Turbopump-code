function [M_glo_mat_shaft_bend,M_glo_mat_shaft_rotary,K_glo_mat_shaft_bend,G_glo_mat_shaft,K_glo_mat_skew] =shaft_glo_matrix(do,di,le,E,rho,v,G,tot_node)

M_glo_mat_shaft_bend=zeros(4*tot_node);
K_glo_mat_shaft_bend=zeros(4*tot_node);
G_glo_mat_shaft=zeros(4*tot_node);
K_glo_mat_skew=zeros(4*tot_node);
M_glo_mat_shaft_rotary=zeros(4*tot_node);

glo_pos_ini=0;

for ele_num=1:1:tot_node-1
%% Shaft Properties

    A=pi*((do(ele_num))^2-(di(ele_num))^2)/4; % area of shaft
    Ie=pi*((do(ele_num))^4-(di(ele_num))^4)/64; % diameteral area moment of inertia of shaft 
    Ip=pi*((do(ele_num))^4-(di(ele_num))^4)/32; %  % polar area moment of inertia of shaft 
    u=di(ele_num)/do(ele_num); %ratio of inner shaft radius to outer shaft radius
     ke=(6*(1+v)*(1+u^2)^2)/((7+6*v)*(1+u^2)^2 + (20+12*v)*u^2);
%     ke=(6*((1+v)^2)*(1+u^2)^2)/((7+12*v+4*v^2)*(1+u^2)^2 + 4*(5+6*v+2*v^2)*u^2);% shear constant
    phi = 12*E*Ie/(ke*G*A*(le(ele_num))^2);

%% Timoshenko Beam 8*8 Elemental Stiffness Matrix

% Stiffness Element Bending Matrix

    K_ele_mat_bend=zeros(8);
    
    K_ele_mat_bend(1,1)=12;  K_ele_mat_bend(1,4)=6*le(ele_num);  K_ele_mat_bend(1,5)=-12;  K_ele_mat_bend(1,8)=6*le(ele_num);
    K_ele_mat_bend(2,2)=12;  K_ele_mat_bend(2,3)=-6*le(ele_num);  K_ele_mat_bend(2,6)=-12;  K_ele_mat_bend(2,7)=-6*le(ele_num);
    
    K_ele_mat_bend(3,2)=-6*le(ele_num); K_ele_mat_bend(3,3)=(4+phi)*(le(ele_num))^2;  K_ele_mat_bend(3,6)=6*le(ele_num);  K_ele_mat_bend(3,7)=(2-phi)*(le(ele_num))^2;
    K_ele_mat_bend(4,1)=6*le(ele_num); K_ele_mat_bend(4,4)=(4+phi)*(le(ele_num))^2;  K_ele_mat_bend(4,5)=-6*le(ele_num);  K_ele_mat_bend(4,8)=(2-phi)*(le(ele_num))^2;
    
    K_ele_mat_bend(5,1)=-12;  K_ele_mat_bend(5,4)=-6*le(ele_num);  K_ele_mat_bend(5,5)=12;  K_ele_mat_bend(5,8)=-6*le(ele_num);
    K_ele_mat_bend(6,2)=-12;  K_ele_mat_bend(6,3)=6*le(ele_num);  K_ele_mat_bend(6,6)=12;  K_ele_mat_bend(6,7)=6*le(ele_num);
    
    K_ele_mat_bend(7,2)=-6*le(ele_num); K_ele_mat_bend(7,3)=(2-phi)*(le(ele_num))^2; K_ele_mat_bend(7,6)=6*le(ele_num); K_ele_mat_bend(7,7)=(4+phi)*(le(ele_num))^2;
    K_ele_mat_bend(8,1)=6*le(ele_num); K_ele_mat_bend(8,4)=(2-phi)*(le(ele_num))^2; K_ele_mat_bend(8,5)=-6*le(ele_num); K_ele_mat_bend(8,8)=(4+phi)*(le(ele_num))^2;
    
    K_ele_mat_bend=(E*Ie/((1+phi)*(le(ele_num))^3))*K_ele_mat_bend;
    
   
% Skew symmetric Stiffness Matrix
    
    K_ele_mat_skew = zeros(8);
    
    K_ele_mat_skew(1,2)=12; K_ele_mat_skew(1,3)=-6*le(ele_num); K_ele_mat_skew(1,6)=-12; K_ele_mat_skew(1,7)=-6*le(ele_num);
    K_ele_mat_skew(2,1)=-12; K_ele_mat_skew(2,4)=-6*le(ele_num); K_ele_mat_skew(2,5)=12; K_ele_mat_skew(2,8)=-6*le(ele_num);
    
    K_ele_mat_skew(3,1)=6*le(ele_num); K_ele_mat_skew(3,4)=(4+phi)*(le(ele_num))^2; K_ele_mat_skew(3,5)=-6*le(ele_num); K_ele_mat_skew(3,8)=(2-phi)*(le(ele_num))^2;
    K_ele_mat_skew(4,2)=6*le(ele_num); K_ele_mat_skew(4,3)=-(4+phi)*(le(ele_num))^2; K_ele_mat_skew(4,6)=-6*le(ele_num); K_ele_mat_skew(4,7)=-(2-phi)*(le(ele_num))^2;
    
    K_ele_mat_skew(5,2)=-12; K_ele_mat_skew(5,3)=6*le(ele_num); K_ele_mat_skew(5,6)=12; K_ele_mat_skew(5,7)=6*le(ele_num);
    K_ele_mat_skew(6,1)=12; K_ele_mat_skew(6,4)=6*le(ele_num); K_ele_mat_skew(6,5)=-12; K_ele_mat_skew(6,8)=6*le(ele_num);
    
    K_ele_mat_skew(7,1)=6*le(ele_num); K_ele_mat_skew(7,4)=(2-phi)*(le(ele_num))^2; K_ele_mat_skew(7,5)=-6*le(ele_num); K_ele_mat_skew(7,8)=(4+phi)*(le(ele_num))^2;
    K_ele_mat_skew(8,2)=6*le(ele_num); K_ele_mat_skew(8,3)=-(2-phi)*(le(ele_num))^2; K_ele_mat_skew(8,6)=-6*le(ele_num); K_ele_mat_skew(8,7)=-(4+phi)*(le(ele_num))^2;
    
    K_ele_mat_skew=(E*Ie/((1+phi)*(le(ele_num))^3))*K_ele_mat_skew;
    
%% Timoshenko Beam 8*8 Mass Matrix


% Mass element bending matrix

    M_ele_mat_bend=zeros(8);
    M_ele_mat_rotary=zeros(8);

    m1=312+588*phi+280*phi^2;  m2=(44+77*phi+35*phi^2)*le(ele_num);  m3=108+252*phi+140*phi^2;
    m4=-(26+63*phi+35*phi^2)*le(ele_num);  m5=(8+14*phi+7*phi^2)*(le(ele_num))^2;  m6=-(6+14*phi+7*phi^2)*(le(ele_num))^2;
    m7=36;  m8=(3-15*phi)*le(ele_num);  m9=(4+5*phi+10*phi^2)*(le(ele_num))^2;  m10=(-1-5*phi+5*phi^2)*(le(ele_num))^2;

    M_ele_mat_bend(1,1)=m1;  M_ele_mat_bend(1,4)=m2;  M_ele_mat_bend(1,5)=m3;  M_ele_mat_bend(1,8)=m4;
    M_ele_mat_bend(2,2)=m1;  M_ele_mat_bend(2,3)=-m2;  M_ele_mat_bend(2,6)=m3;  M_ele_mat_bend(2,7)=-m4;
    
    M_ele_mat_bend(3,2)=-m2; M_ele_mat_bend(3,3)=m5;  M_ele_mat_bend(3,6)=m4;  M_ele_mat_bend(3,7)=m6;
    M_ele_mat_bend(4,1)=m2; M_ele_mat_bend(4,4)=m5;  M_ele_mat_bend(4,5)=-m4;  M_ele_mat_bend(4,8)=m6;
    
    M_ele_mat_bend(5,1)=m3;  M_ele_mat_bend(5,4)=-m4;  M_ele_mat_bend(5,5)=m1;  M_ele_mat_bend(5,8)=-m2;
    M_ele_mat_bend(6,2)=m3;  M_ele_mat_bend(6,3)=m4;  M_ele_mat_bend(6,6)=m1;  M_ele_mat_bend(6,7)=m2;
    
    M_ele_mat_bend(7,2)=-m4; M_ele_mat_bend(7,3)=m6; M_ele_mat_bend(7,6)=m2; M_ele_mat_bend(7,7)=m5;
    M_ele_mat_bend(8,1)=m4; M_ele_mat_bend(8,4)=m6; M_ele_mat_bend(8,5)=-m2; M_ele_mat_bend(8,8)=m5;
    
    M_ele_mat_bend=(rho*A*le(ele_num)/(840*(1+phi)^2))*M_ele_mat_bend;
    
    % Effect of rotary inertia
    
    M_ele_mat_rotary(1,1)=m7;  M_ele_mat_rotary(1,4)=m8;  M_ele_mat_rotary(1,5)=-m7;  M_ele_mat_rotary(1,8)=m8;
    M_ele_mat_rotary(2,2)=m7;  M_ele_mat_rotary(2,3)=-m8;  M_ele_mat_rotary(2,6)=-m7;  M_ele_mat_rotary(2,7)=-m8;
    
    M_ele_mat_rotary(3,2)=-m8; M_ele_mat_rotary(3,3)=m9;  M_ele_mat_rotary(3,6)=m8;  M_ele_mat_rotary(3,7)=m10;
    M_ele_mat_rotary(4,1)=m8; M_ele_mat_rotary(4,4)=m9;  M_ele_mat_rotary(4,5)=-m8;  M_ele_mat_rotary(4,8)=m10;
    
    M_ele_mat_rotary(5,1)=-m7;  M_ele_mat_rotary(5,4)=-m8;  M_ele_mat_rotary(5,5)=m7;  M_ele_mat_rotary(5,8)=-m8;
    M_ele_mat_rotary(6,2)=-m7;  M_ele_mat_rotary(6,3)=m8;  M_ele_mat_rotary(6,6)=m7;  M_ele_mat_rotary(6,7)=m8;
    
    M_ele_mat_rotary(7,2)=-m8; M_ele_mat_rotary(7,3)=m10; M_ele_mat_rotary(7,6)=m8; M_ele_mat_rotary(7,7)=m9;
    M_ele_mat_rotary(8,1)=m8; M_ele_mat_rotary(8,4)=m10; M_ele_mat_rotary(8,5)=-m8; M_ele_mat_rotary(8,8)=m9;
    
    M_ele_mat_rotary=(rho*Ie/(30*(1+phi)^2*le(ele_num)))*M_ele_mat_rotary;
    
%% Timoshenko Beam 12*12 Gyroscopic Matrix

    G_ele_mat_shaft=zeros(8);
    
    g1=36;  g2=(3-15*phi)*le(ele_num);  g3=(4+5*phi+10*phi^2)*(le(ele_num))^2;
    g4=(-1-5*phi+5*phi^2)*(le(ele_num))^2;
    
    G_ele_mat_shaft(1,2)=g1;  G_ele_mat_shaft(1,3)=-g2;  G_ele_mat_shaft(1,6)=-g1;  G_ele_mat_shaft(1,7)=-g2;
    G_ele_mat_shaft(2,1)=-g1;  G_ele_mat_shaft(2,4)=-g2;  G_ele_mat_shaft(2,5)=g1;  G_ele_mat_shaft(2,8)=-g2;
    
    G_ele_mat_shaft(3,1)=g2; G_ele_mat_shaft(3,4)=g3; G_ele_mat_shaft(3,5)=-g2; G_ele_mat_shaft(3,8)=g4;
    G_ele_mat_shaft(4,2)=g2; G_ele_mat_shaft(4,3)=-g3; G_ele_mat_shaft(4,6)=-g2; G_ele_mat_shaft(4,7)=-g4;
    
    G_ele_mat_shaft(5,2)=-g1;  G_ele_mat_shaft(5,3)=g2;  G_ele_mat_shaft(5,6)=g1;  G_ele_mat_shaft(5,7)=g2;
    G_ele_mat_shaft(6,1)=g1;  G_ele_mat_shaft(6,4)=g2;  G_ele_mat_shaft(6,5)=-g1;  G_ele_mat_shaft(6,8)=g2;
    
    G_ele_mat_shaft(7,1)=g2; G_ele_mat_shaft(7,4)=g4; G_ele_mat_shaft(7,5)=-g2; G_ele_mat_shaft(7,8)=g3;
    G_ele_mat_shaft(8,2)=g2; G_ele_mat_shaft(8,3)=-g4; G_ele_mat_shaft(8,6)=-g2; G_ele_mat_shaft(8,7)=-g3;
    
    G_ele_mat_shaft=(rho*Ie/(15*(1+phi)^2*le(ele_num)))*G_ele_mat_shaft;
    
    
%% Global Matrix

    for row_no=1:8
        for column_no=1:8
            M_glo_mat_shaft_bend(row_no+glo_pos_ini,column_no+glo_pos_ini)=M_glo_mat_shaft_bend(row_no+glo_pos_ini,column_no+glo_pos_ini)+M_ele_mat_bend(row_no,column_no);
            M_glo_mat_shaft_rotary(row_no+glo_pos_ini,column_no+glo_pos_ini)=M_glo_mat_shaft_rotary(row_no+glo_pos_ini,column_no+glo_pos_ini)+M_ele_mat_rotary(row_no,column_no);
            K_glo_mat_shaft_bend(row_no+glo_pos_ini,column_no+glo_pos_ini)=K_glo_mat_shaft_bend(row_no+glo_pos_ini,column_no+glo_pos_ini)+K_ele_mat_bend(row_no,column_no);
            G_glo_mat_shaft(row_no+glo_pos_ini,column_no+glo_pos_ini)=G_glo_mat_shaft(row_no+glo_pos_ini,column_no+glo_pos_ini)+G_ele_mat_shaft(row_no,column_no);
            K_glo_mat_skew(row_no+glo_pos_ini,column_no+glo_pos_ini)=K_glo_mat_skew(row_no+glo_pos_ini,column_no+glo_pos_ini)+K_ele_mat_skew(row_no,column_no);
        end
    end
    glo_pos_ini=glo_pos_ini+4;

end
end
