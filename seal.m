function [K_glo_mat_seal,C_glo_mat_seal,M_glo_mat_seal,K,k,C,c,m]=seal(N,rho_fluid,mu_fluid,Clearence,len_seal,v0,D,zeta,delp,seal_node,tot_node)

w=2*pi*N/60;

M_glo_mat_seal=zeros(4*tot_node);
C_glo_mat_seal=zeros(4*tot_node);
K_glo_mat_seal=zeros(4*tot_node);
R=D/2;
n0=0.066; % Co-efficients for hir's turbulence equation
m0=-0.25; % Co-efficients for hir's turbulence equation
err=100;
V=0;
while err>10
    V=V+0.0001;
    Ra0=rho_fluid*V*Clearence/mu_fluid;
    Rc0=rho_fluid*(R*w)*Clearence/mu_fluid;
    b=Ra0/Rc0; % dimensionless co-efficient
    lmda=n0*(Ra0^m0)*(1+1/(4*b^2))^((1+m0)/2);
    sigma=lmda*len_seal/Clearence;
    beta0=1/(1+(16/7*b)^2); % dimensionless co-efficient
    beta=1/(1+4*b^2);
    a=sigma*(1+beta0*(1+m0)); % dimensionless co-efficient
    B=1+4*(b^2)*beta*(1+m0); % dimensionless co-efficient
    delp_new=(1+zeta+2*sigma)*rho_fluid*V^2/2;
    err=abs(delp-delp_new);
end

T=len_seal/V;
E=(1+zeta)/(2*(1+zeta+B*sigma)); % dimensionless co-efficient
r=2*sigma^2/(1+zeta+2*sigma); % dimensionless co-efficient


K=(2*sigma^2)/(1+zeta+2*sigma)*(E*(1-m0)-(w*T)^2/(4*sigma)*(0.5*(1/6+E)+2*v0/a*((E+1/a^2)*(1-exp(-a))-(0.5+1/a)*exp(-a))))*pi*R*delp/lmda;
k=w*T*sigma^2/(1+zeta+2*sigma)*(E/sigma+B/2*(1/6+E)+2*v0/a*(E*B+(1/sigma-B/a)*((1-exp(-a))*(E+1/2+1/a)-1)))*pi*R*delp/lmda;
C=2*sigma^2/(1+zeta+2*sigma)*(E/sigma+B/2*(1/6+E))*T*pi*R*delp/lmda;
c=2*sigma*(w*T)/(1+zeta+2*sigma)*(1/2*(1/6+E)+v0/a*((1-exp(-a))*(E+1/2+1/a^2)-(1/2+exp(-a)/a)))*T*pi*R*delp/lmda;
m=sigma*(1/6+E)/(1+zeta+2*sigma)*T^2*pi*R*delp/lmda;

K_ele_mat_seal=[K k; -k K];
C_ele_mat_seal=[C c; -c C];
M_ele_mat_seal=[m 0; 0 m];
    for row_no=1:2
        for column_no=1:2
            M_glo_mat_seal(row_no+4*(seal_node-1),column_no+4*(seal_node-1))=M_glo_mat_seal(row_no+4*(seal_node-1),column_no+4*(seal_node-1))+M_ele_mat_seal(row_no,column_no);
            C_glo_mat_seal(row_no+4*(seal_node-1),column_no+4*(seal_node-1))=C_glo_mat_seal(row_no+4*(seal_node-1),column_no+4*(seal_node-1))+C_ele_mat_seal(row_no,column_no);
            K_glo_mat_seal(row_no+4*(seal_node-1),column_no+4*(seal_node-1))=K_glo_mat_seal(row_no+4*(seal_node-1),column_no+4*(seal_node-1))+K_ele_mat_seal(row_no,column_no);
        end
    end

end