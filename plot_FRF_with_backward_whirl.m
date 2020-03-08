function []=plot_FRF_with_backward_whirl(freq_rpm,response,desired_node)
freq_step=freq_rpm(2)-freq_rpm(1);
legend_text=[];
node_dir=[desired_node 1; desired_node 2];
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
backward_whirl=[];
dont_plot=0;
cnt=0;
for N=1:length(freq_rpm)
    resp_x=response(4*desired_node-3,N);
    resp_y=response(4*desired_node-2,N);
    if whirl(resp_x,resp_y)<0
        cnt=cnt+1;
        backward_whirl(cnt,1)=freq_rpm(N);
    end
end

if isempty(backward_whirl)==1
    dont_plot=1;
end

back_whirl_padding=[0;backward_whirl;0];

% subplot(2,1,1)
semilogy(freq_rpm,mag_data(:,1));
hold on
semilogy(freq_rpm,mag_data(:,2));
hold on
if dont_plot==0
    for N=2:length(back_whirl_padding)-1
        if back_whirl_padding(N)+freq_step==back_whirl_padding(N+1)
            start1=back_whirl_padding(N);
            end1=back_whirl_padding(N+1);
            v=[start1 min(ylim); end1 min(ylim); end1 max(ylim); start1 max(ylim)];
            f=[1 2 3 4];
            patch('faces',f,'vertices',v,'EdgeColor','None','FaceColor','yellow')
        else
            start1=back_whirl_padding(N);
            v=[start1-0.5 min(ylim);start1+0.5 min(ylim); start1+0.5 max(ylim); start1-0.5 max(ylim)];
            f=[1 2 3 4];
            patch('faces',f,'vertices',v,'EdgeColor','None','FaceColor','yellow')
        end
            
    end
end
semilogy(freq_rpm,mag_data(:,1));
hold on
semilogy(freq_rpm,mag_data(:,2));
grid
hold off
xlabel('Rotor spin speed (rpm)')
ylabel('Response magnitude (m)')
set(gca,'fontsize',24)
% legend(legend_text(1,:),legend_text(2,:))
legend('Turbine node,x','Turbine node,y')

% subplot(2,1,2)
% plot(freq_rpm,phase_data(:,1));
% hold on
% plot(freq_rpm,phase_data(:,2));
% hold on
% if dont_plot==0
%     for N=2:length(back_whirl_padding)-1
%         if back_whirl_padding(N)+freq_step==back_whirl_padding(N+1)
%             start1=back_whirl_padding(N);
%             end1=back_whirl_padding(N+1);
%             v=[start1 min(ylim); end1 min(ylim); end1 max(ylim); start1 max(ylim)];
%             f=[1 2 3 4];
%             patch('faces',f,'vertices',v,'EdgeColor','None','FaceColor','yellow')
%         else
%             start1=back_whirl_padding(N);
%             v=[start1-0.5 min(ylim);start1+0.5 min(ylim); start1+0.5 max(ylim); start1-0.5 max(ylim)];
%             f=[1 2 3 4];
%             patch('faces',f,'vertices',v,'EdgeColor','None','FaceColor','yellow')
%         end
%             
%     end
% end
% plot(freq_rpm,phase_data(:,1));
% hold on
% plot(freq_rpm,phase_data(:,2));
% grid
% hold off
% xlabel('Rotor spin speed (rev/min)')
% ylabel('Phase (degrees)')
% % legend(legend_text(1,:),legend_text(2,:))
% legend('Turbine node,x','Turbine node,y')
end