function [] = campbell(nat_freq,N,kappa)
nat_freq_Hz=nat_freq/(2*pi);
whirl_dirn=zeros(length(nat_freq_Hz),length(N));
N_Hz=(N/60)*[1];
plot(N,nat_freq_Hz,'k-',N,N_Hz,'b--');
ylim([0 N_Hz(end)])
set(gca,'fontsize',24)
grid;
pos_whirl=[];
neg_whirl=[];
mix_whirl=[];
if nargin>2
for i=1:length(kappa(1,:,1)) % from kappa
    for j=2:length(N)
        modal_kappa=kappa(:,i,j);
        n=0;
        p=0;
        for k=1:length(modal_kappa)
            if modal_kappa(k)<0
                n=n+1;
            end
            if modal_kappa(k)>0
                p=p+1;
            end
        end
        
        if n>0 && p==0
%             whirl_dirn(i,j)=-1;
            neg_whirl=[neg_whirl; N(j) nat_freq_Hz(i,j)];
        end
        if p>0 && n==0
%             whirl_dirn(i,j)=1;
            pos_whirl=[pos_whirl; N(j) nat_freq_Hz(i,j)];
        end
        if p>0 && n>0
%             whirl_dirn(i,j)=0;
            mix_whirl=[mix_whirl; N(j) nat_freq_Hz(i,j)];
        end
    end
end
hold on 
if isempty(mix_whirl)~=1
    camp_diag=plot(pos_whirl(:,1),pos_whirl(:,2),'r.',neg_whirl(:,1),neg_whirl(:,2),'g.',mix_whirl(:,1),mix_whirl(:,2),'k.');
    set(camp_diag(3),'Markersize',10)
    xlabel('Speed in RPM')
    ylabel('Frequency (in Hz)')
%     title('Forward whirl is red, Backward whirl is green & Mixed whirl is black.')
else
    camp_diag=plot(pos_whirl(:,1),pos_whirl(:,2),'r.',neg_whirl(:,1),neg_whirl(:,2),'g.');
    xlabel('Speed in RPM')
    ylabel('Frequency (in Hz)')
%     title('Forward whirl is red & Backward whirl is green.')
end
set(camp_diag(1),'Markersize',10)
set(camp_diag(2),'Markersize',10)
% set(axes.camp_diag,'fontsize',24)
hold off
end
end
