function []=orbit_plot(Mode,orbit_node,titletext,eigen_value)
cnt=0;
for m=1:length(orbit_node)
cnt=cnt+1;
orbit_x=Mode(4*(orbit_node(m)-1)+1);
orbit_y=Mode(4*(orbit_node(m)-1)+2);

wt=(0:1:340)*(pi)/180;

loc_x(:,cnt)=real(orbit_x*exp(1i*sign(eigen_value)*wt));
loc_y(:,cnt)=real(orbit_y*exp(1i*sign(eigen_value)*wt));
end

max_loc_x=max(loc_x);
max_loc_y=max(loc_y);

max_val=1.2*max([max_loc_x max_loc_y]);

plot(loc_x(:,1),loc_y(:,1),'k-',loc_x(1,1),loc_y(1,1),'kx',loc_x(341,1),loc_y(341,1),'kd')
axis square
axis([-max_val max_val -max_val max_val])
axis off

hold on
if nargin >= 3
    title(titletext)
end
for j=2:length(orbit_node)
    plot(loc_x(:,j),loc_y(:,j),'k--',loc_x(1,j),loc_y(1,j),'kx',loc_x(end,j),loc_y(end,j),'kd')
end
hold off
set(gca,'fontsize',24)
end

