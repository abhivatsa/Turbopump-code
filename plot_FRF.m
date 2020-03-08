function []=plot_FRF(freq_rpm,response,node_dir,unstab_speed)
legend_text=[];
for i =1:length(node_dir(:,1))
    node=node_dir(i,1);
    dir=node_dir(i,2);
    if dir == 1
        node_string=['Node' num2str(node) ',x'];
    end
    if dir == 2
        node_string=['Node' num2str(node) ',y'];
    end
    legend_text=[legend_text; node_string];
    mag_data(:,i)=abs(response((4*node)-4+dir,:));
    phase_data(:,i)=-180*angle(response((4*node)-4+dir,:))/pi;
end
% subplot(2,1,1)
semilogy(freq_rpm,mag_data(:,1));
hold on
semilogy(freq_rpm,mag_data(:,2));
hold on
v=[unstab_speed min(ylim); freq_rpm(end) min(ylim); freq_rpm(end) max(ylim); unstab_speed max(ylim)];
f=[1 2 3 4];
patch('faces',f,'vertices',v,'EdgeColor','None','FaceColor','red')
alpha(0.3)
hold off
legend('Impeller node,x','Turbine node,x')
xlabel('Rotor spin speed (rpm)')
ylabel('Response magnitude (m)')
set(gca,'fontsize',24)
grid
% subplot(2,1,2)
% plot(freq_rpm,phase_data(:,1));
% hold on
% plot(freq_rpm,phase_data(:,2));
% hold off
% % legend(legend_text(1,:),legend_text(2,:))
% legend('Impeller node,x','Turbine node,x')
% xlabel('Rotor spin speed (rev/min)')
% ylabel('Phase (degrees)')
% set(gca,'fontsize',15)
% grid
end