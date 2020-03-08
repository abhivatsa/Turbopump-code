function []=mode_shapes(Mode,L,nat_freq)
shaft_stat=zeros(length(L),1);
nat_freq_Hz=nat_freq/(2*pi);
len_zero=length(Mode(:,1));

mode_x = Mode(1:4:len_zero);
mode_y = Mode(2:4:len_zero);
    
[xmax,index_xmax] = max(abs(mode_x));
[ymax,index_ymax] = max(abs(mode_y));
    
if ymax > 0.4*xmax
    Mode = Mode/mode_y(index_ymax);
else
    Mode = Mode/mode_x(index_xmax);
end

Mode_x = Mode(1:4:len_zero);
Mode_y= Mode(2:4:len_zero);
    
wt=(0:1:180)*(2*pi)/180;
for node=1:length(Mode_x(:,1))
    mode_loc_x=(real(Mode_x*exp(1i*wt)))';
    mode_loc_y=(real(Mode_y*exp(1i*wt)))';
end
    
Mode_x=real(Mode_x);
Mode_y=real(Mode_y);
% Normalization of modes
max_xy=max(abs([Mode_x; Mode_y]));
% Mode changing from vector to array Modes
Mode_x=Mode_x/max_xy;
Mode_y=Mode_y/max_xy;
mode_loc_x=mode_loc_x/max_xy;
mode_loc_y=mode_loc_y/max_xy;


node=ones(181,1)*L;

p=plot3(L,Mode_x,Mode_y,'-k',node,mode_loc_x,mode_loc_y,'-b');
set(p(1),'LineWidth',1.5);
set(p(2),'LineWidth',0.25);
title(['Nat Freq = ' num2str(nat_freq_Hz) 'Hz'],'HorizontalAlignment','center')
set(gca,'fontsize',16)
axis off
