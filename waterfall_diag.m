function [freq,speed,recep] = waterfall_diag(resp,N_rpm_coast,Time)
length(Time.')
data_set = zeros(length(Time.'),3);
length(data_set)
data_set(:,2)=resp;
data_set(:,1)=Time.';
data_set(:,3)=N_rpm_coast.';

ctr=0;
for i=1001:500:(length(data_set(:,2))-300)
    ctr=ctr+1
    data_fft=data_set(i-200:i+200,:);

    Time = data_fft(:,1);
    samp_time = Time(2)-Time(1); 
    samp_freq = 1/samp_time;

    len = length(Time);


    measure_data = data_fft(:,2);

    measure_data_var = measure_data - mean(measure_data);

    Four_tran = fft(measure_data_var);
    Four_tran_2=abs(Four_tran/len);
    Four_tran_1 = Four_tran_2(1:(floor(len/2))+1);
    Four_tran_1(2:end-1)=2*Four_tran_1(2:end-1);
    freq = samp_freq*(0:(len/2))/len;
    freq=freq';

    recep(ctr,:)=Four_tran_1.';
    speed(ctr,1)=data_set(i,3);


end
end
% axis square