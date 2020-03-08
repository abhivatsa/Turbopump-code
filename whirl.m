function [kappa]=whirl(u,v)
ru=abs(u);
rv=abs(v);
nu=angle(u);
nv=angle(v);
for i=1:length(u)
    if ru(i)*rv(i)<1e-16
        kappa(i)=0;
    else
        T=[ru(i)*cos(nu(i)) -ru(i)*sin(nu(i)); rv(i)*cos(nv(i)) -rv(i)*sin(nv(i))];
        H=inv(T*T.');
        lmda=eig(H);
        lmda_update=sort(lmda);
        Kappa_mag=sqrt(lmda_update(1,1)/lmda_update(2,1));
        phase_diff=mod(nv(i)-nu(i),2*pi);
        if phase_diff>0 && phase_diff<pi 
            Kappa_mag=-Kappa_mag;
        end
        kappa(i)=Kappa_mag;
    end
end
end



