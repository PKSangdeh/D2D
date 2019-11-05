clear;clear;
clc;
close all;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_ant_ap = 4;
num_valid_sc = 52;
num_subcarrier = 64;
num_symbol_per_frame = 20;
num_frame = 20;
disp('Cellular Communications')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LTF
LTF_DAT_GRP = zeros(10, 52);
LTF_DAT_GRP(1,:) = [-1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1 1 1 -1 -1 -1];
LTF_DAT_GRP(2,:) = [-1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1 1 1 -1 -1 -1];
LTF_DAT_GRP(3,:) = [-1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1 1 1 -1 -1 -1];
LTF_DAT_GRP(4,:) = [-1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1 1 1 -1 -1 -1];
LTF_DAT_GRP(5,:) = [-1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1 1 1 -1 -1 -1];
LTF_DAT_GRP(6,:) = [-1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1 1 1 -1 -1 -1];
LTF_DAT_GRP(7,:) = [-1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1 1 1 -1 -1 -1];
LTF_DAT_GRP(8,:) = [-1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1 1 1 -1 -1 -1];
LTF_DAT_GRP(9,:) = [-1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1 1 1 -1 -1 -1];


%%%% pilot
POS_VALID_SC = [39:64 2:27];
POS_PAYLOAD_SC = [1:5 7:19 21:32 34:46 48:52];
POS_PILOT_SC = [6 20 33 47];
PILOT_DAT_GRP = zeros(4,num_symbol_per_frame);
PILOT_DAT_GRP(1,:) = [-1 -1 -1 -1 1 1 1 -1 1 -1 -1 1 1 1 1 -1 -1 1 -1 1];
PILOT_DAT_GRP(2,:) = [1 -1 -1 1 -1 1 1 1 1 1 1 -1 1 1 1 -1 -1 1 1 1];
PILOT_DAT_GRP(3,:) = [-1 -1 1 1 -1 -1 -1 1 1 1 -1 -1 1 -1 1 1 1 1 1 1];
PILOT_DAT_GRP(4,:) = [1 1 1 -1 1 1 -1 -1 -1 -1 1 -1 -1 -1 -1 -1 -1 1 1 -1];

PREAMABLE = [1,1,-1,1,-1,-1,-1,-1,1,1,1,1,-1,1,-1,1,-1,1,-1,1,1,1,1,1,-1,1,-1,1,1,1,-1,-1,1,-1,1,1,1,1,1,1,1,-1,-1,-1,1,-1,1,-1,-1,-1,-1,1,-1,-1,-1,1,1,1,1,-1,-1,1,1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,-1,-1,-1,1,-1,1,1,1,-1,1,-1,-1,-1,-1,-1,-1,1,1,-1,1,-1,-1,1,-1,1];
PREAMABLE = [PREAMABLE PREAMABLE];
PREAMABLE = PREAMABLE/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate freq data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
time_data_matrix_all_ant = zeros(num_ant_ap, 80, num_symbol_per_frame);
for ant_idx = 1:num_ant_ap
    % generate freq data
    freq_data_matrix = zeros(num_valid_sc, num_symbol_per_frame);
    for sym_idx = 3*num_ant_ap+1:num_symbol_per_frame %2*num_ant_ap+1
        freq_data_matrix(:,sym_idx) = (1/2)*(sign(randn(num_valid_sc,1))+1i*sign(randn(num_valid_sc,1)));
    end
    freq_data_matrix(:,3*ant_idx-2) = transpose(LTF_DAT_GRP(ant_idx,:));
    freq_data_matrix(:,3*ant_idx-1) = transpose(LTF_DAT_GRP(ant_idx,:));
    freq_data_matrix(:,3*ant_idx-0) = transpose(LTF_DAT_GRP(ant_idx,:));
    % Convert to time domain
    freq_data_matrix_tmp = zeros(num_subcarrier, num_symbol_per_frame);
    freq_data_matrix_tmp(POS_VALID_SC,:) = freq_data_matrix;
    time_data_matrix = ifft(freq_data_matrix_tmp)*sqrt(num_subcarrier);
    time_data_matrix = [time_data_matrix(end-15:end,:); time_data_matrix];
    % % add preamble
    time_data_matrix_all_ant(ant_idx,:,:) = time_data_matrix;
end
 



    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build orthogonal tx signal for channel estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_data_matrix_grp = zeros(num_ant_ap, 80*num_symbol_per_frame);
for ant_idx = 1:num_ant_ap
    time_data_matrix_grp(ant_idx,:) = reshape(time_data_matrix_all_ant(ant_idx,:,:), 1, []);
end
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add preamble
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ampl_avg = mean(mean(mean(abs(time_data_matrix_grp))));
preamable_matrix = repmat(num_ant_ap*ampl_avg*PREAMABLE, num_ant_ap, 1);
tx_signal_in_time = [preamable_matrix time_data_matrix_grp];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% how much data to write
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_data_matrix = repmat(tx_signal_in_time, 1, 20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write to the file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:num_ant_ap
    %fileID= fopen(['/tmp/data/cali_ap_signal_tx' num2str(ii-1) '.dat'], 'wb');
    fileID= fopen(['signals/cali_ap_signal_tx' num2str(ii-1) '.dat'], 'wb');
    if (fileID < 0)
        disp('Error: fail to open files!');
        pause;
    end
    data_tmp_c = time_data_matrix(ii,:);
    data_tmp_f = [real(data_tmp_c); imag(data_tmp_c)];
    data_tmp_f = data_tmp_f(:);
    fwrite(fileID,data_tmp_f,'float');
    fclose(fileID);
end



%% Parameters
num_ant_sta = 1;
num_valid_sc = 52;
num_subcarrier = 64;
num_symbol_per_frame = 20;
num_frame = 20;
num_sample_process = 50000;
num_sample_shift = 2e6;
num_move_sample = -1;

%% STF 
STF_POS = [3  7  11  15  19  23  30  34  38  42  46  50];
STF_DAT_GRP = zeros(10, 52);    
STF_DAT_GRP(1,STF_POS) = sqrt(2)*[-1-1i  -1+1i  -1-1i  -1-1i  -1-1i  1-1i  -1-1i  -1-1i  1+1i  1-1i  1+1i  -1+1i];
STF_DAT_GRP(2,STF_POS) = sqrt(2)*[-1+1i  -1+1i  1-1i  1-1i  -1-1i  1-1i  -1+1i  1+1i  1-1i  -1+1i  1-1i  -1+1i];
STF_DAT_GRP(3,STF_POS) = sqrt(2)*[-1+1i  1+1i  -1+1i  1+1i  1-1i  1+1i  -1+1i  -1+1i  1-1i  1-1i  -1+1i  -1+1i];
STF_DAT_GRP(4,STF_POS) = sqrt(2)*[1+1i  -1+1i  -1+1i  -1+1i  -1-1i  1+1i  1+1i  1-1i  1-1i  -1+1i  -1+1i  -1-1i];
STF_DAT_GRP(5,STF_POS) = sqrt(2)*[-1+1i  -1+1i  1+1i  -1+1i  -1-1i  1+1i  1-1i  -1-1i  -1-1i  -1+1i  -1+1i  1-1i];
STF_DAT_GRP(6,STF_POS) = sqrt(2)*[1+1i  1-1i  -1-1i  1+1i  -1+1i  1+1i  1-1i  -1-1i  1+1i  -1-1i  -1+1i  -1+1i];
STF_DAT_GRP(7,STF_POS) = sqrt(2)*[1+1i  1+1i  1+1i  1-1i  -1-1i  1-1i  -1-1i  1-1i  1+1i  1+1i  -1-1i  1-1i];
STF_DAT_GRP(8,STF_POS) = sqrt(2)*[-1+1i  1-1i  1+1i  -1+1i  1+1i  1-1i  -1-1i  -1-1i  -1+1i  1-1i  -1+1i  -1-1i];
STF_DAT_GRP(9,STF_POS) = sqrt(2)*[1-1i  -1+1i  1-1i  -1+1i  1-1i  1+1i  1+1i  -1-1i  -1-1i  -1+1i  1+1i  -1-1i];
STF_DAT_GRP(10,STF_POS) = sqrt(2)*[1-1i  1-1i  -1-1i  1+1i  1+1i  1+1i  1+1i  1-1i  -1-1i  -1+1i  -1+1i  -1-1i];  
%% LTF
LTF_DAT_GRP = zeros(10, 52);
LTF_DAT_GRP(1,:) = [-1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1 1 1 -1 -1 -1];
LTF_DAT_GRP(2,:) = [-1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1 1 1 -1 -1 -1];
LTF_DAT_GRP(3,:) = [-1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1 1 1 -1 -1 -1];
LTF_DAT_GRP(4,:) = [-1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1 1 1 -1 -1 -1];
LTF_DAT_GRP(5,:) = [-1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1 1 1 -1 -1 -1];
LTF_DAT_GRP(6,:) = [-1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1 1 1 -1 -1 -1];
LTF_DAT_GRP(7,:) = [-1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1 1 1 -1 -1 -1];
LTF_DAT_GRP(8,:) = [-1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1 1 1 -1 -1 -1];
LTF_DAT_GRP(9,:) = [-1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1 1 1 -1 -1 -1];
%% Pilot
POS_VALID_SC = [39:64 2:27];
POS_PAYLOAD_SC = [1:5 7:19 21:32 34:46 48:52];
POS_PILOT_SC = [6 20 33 47];
PILOT_DAT_GRP = zeros(4,num_symbol_per_frame);
PILOT_DAT_GRP(1,:) = [-1 -1 -1 -1 1 1 1 -1 1 -1 -1 1 1 1 1 -1 -1 1 -1 1];
PILOT_DAT_GRP(2,:) = [1 -1 -1 1 -1 1 1 1 1 1 1 -1 1 1 1 -1 -1 1 1 1];
PILOT_DAT_GRP(3,:) = [-1 -1 1 1 -1 -1 -1 1 1 1 -1 -1 1 -1 1 1 1 1 1 1];
PILOT_DAT_GRP(4,:) = [1 1 1 -1 1 1 -1 -1 -1 -1 1 -1 -1 -1 -1 -1 -1 1 1 -1];
%% Preamble
PREAMABLE = [1,1,-1,1,-1,-1,-1,-1,1,1,1,1,-1,1,-1,1,-1,1,-1,1,1,1,1,1,-1,1,-1,1,1,1,-1,-1,1,-1,1,1,1,1,1,1,1,-1,-1,-1,1,-1,1,-1,-1,-1,-1,1,-1,-1,-1,1,1,1,1,-1,-1,1,1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,-1,-1,-1,1,-1,1,1,1,-1,1,-1,-1,-1,-1,-1,-1,1,1,-1,1,-1,-1,1,-1,1];
PREAMABLE = [PREAMABLE PREAMABLE];
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ==================== UL Channels ==================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
signal_stream_original = zeros(num_ant_ap, num_sample_process);
for ii = 1:num_ant_ap
    fileID = fopen(['signals/cali_ap_signal_rx' num2str(ii-1) '.dat'], 'rb');
    if (fileID < 0)
        error('Error: fail to open files!');
    end
    frewind(fileID);
    fseek(fileID, num_sample_shift, 'bof');
    data_ant_f = fread(fileID, 2*num_sample_process, 'float');
    data_ant_c = transpose(data_ant_f(1:2:end) + 1i*data_ant_f(2:2:end));
    signal_stream_original(ii,:) = data_ant_c(1:num_sample_process);
    fclose(fileID);
end
signal_stream_original =  signal_stream_original*1000;
signal_stream_combined = sum(signal_stream_original,1);

%% time synchronization
% Cross correlation
data_freq = zeros(1, num_subcarrier);
data_freq(POS_VALID_SC) = LTF_DAT_GRP(1,:);
data_time = ifft(data_freq)*sqrt(num_subcarrier);
ltf_signal = [data_time(end-15:end) data_time];
local_signal = [ltf_signal  ltf_signal];
cross_correlation_arr = zeros(1,num_sample_process);
for ii = 200:10000
    signal_segment = signal_stream_combined(ii:ii+160-1);
    cross_correlation_arr(ii) =  local_signal*signal_segment'/(norm(local_signal)*norm(signal_segment));
end
% figure; hold on; grid on; 
% plot(abs(cross_correlation_arr),'k*-');
% xlabel('sample index');
% ylabel('correlation');
% title('Cross correlation (UL)');
% box on;
% axis([0 length(signal_stream_combined) 0 1]);

%% Capture frames
[~, begofframe] = max(abs(cross_correlation_arr(1:1600)));
begofframe = begofframe(1) + 1600 - 160 - num_move_sample;
signal_stream_buf = signal_stream_original(:, begofframe:end);
signal_stream_combined = signal_stream_combined(:, begofframe:end);


%% frequency synchronization
freq_offset1 = 0;
freq_offset2 = 0;
for ii = 0:1600:length(signal_stream_combined)-1600
    % CP
    for jj = 1:20
        freq_offset1 = freq_offset1 + signal_stream_combined(ii+80*(jj-1)+1:ii+80*(jj-1)+16)*signal_stream_combined(ii+80*(jj-1)+65:ii+80*(jj-1)+80)';
    end
    % LTF
    freq_offset2 = freq_offset2 + signal_stream_combined(ii+177:ii+240)*signal_stream_combined(ii+257:ii+320)';
end
angle_os1 = angle(freq_offset1)/64;
angle_os2 = angle(freq_offset2)/80;
freq_comp = exp(1i*angle_os1*(1:length(signal_stream_combined)));
freq_comp_fram = repmat(freq_comp, size(signal_stream_buf,1), 1);
signal_stream_buf = signal_stream_buf.*freq_comp_fram;

%disp(['uplink freq offset: ' num2str(angle_os1) '  ' num2str(angle_os2)]);
angle_os_ul = angle_os1;
angle_os_dl = -angle_os_ul;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% channel estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chan_coeff_matrix = zeros(num_ant_ap, num_valid_sc, num_frame);    
for frm_idx = 1:num_frame
    % reoganize the signal
    signal_frame = signal_stream_buf(:, 1600*(frm_idx-1)+1:1600*frm_idx);
    signal_matrix = zeros(num_ant_ap, num_valid_sc, num_symbol_per_frame);
    for rx_idx = 1:num_ant_ap
        data_t1 = signal_frame(rx_idx,:);
        data_t2 = reshape(data_t1, num_subcarrier+16, num_symbol_per_frame);
        data_t3 = data_t2(17:end,:);
        data_f1 = fft(data_t3)/sqrt(num_subcarrier);
        data_f2 = data_f1(POS_VALID_SC,:);
        signal_matrix(rx_idx,:,:) = data_f2;
    end
    % estimate channel
    for rx_idx = 1:num_ant_ap
        signal_ltf_avg = (signal_matrix(rx_idx, :, 3) + signal_matrix(rx_idx, :, 4))/2;
        chan_coeff_matrix(rx_idx,:,frm_idx) = signal_ltf_avg./LTF_DAT_GRP(1, :);
    end
end

% for ant_idx = 1:num_ant_ap
%     figure; 
%     subplot(2, 1, 1);
%     title(['UL channel Amp - Ant ' num2str(173+ant_idx)]);    
%     hold on;
%     grid on;    
%     for sym_idx = 1:num_symbol_per_frame
%         plot(abs(chan_coeff_matrix(ant_idx,:,sym_idx)),'*-');
%     end
%     subplot(2, 1, 2);
%     title(['UL channel Pha - Ant ' num2str(173+ant_idx)]);    
%     hold on;
%     grid on;    
%     for sym_idx = 1:num_symbol_per_frame
%         plot(unwrap(angle((chan_coeff_matrix(ant_idx,:,sym_idx)))),'*-');
%     end    
% end
ul_chan_coeff_matrix = chan_coeff_matrix;
ul_chan=ul_chan_coeff_matrix;
%save ul_chan.mat ul_chan

%==================================================
%================= DL Channel =====================
%==================================================
signal_stream_original = zeros(num_ant_sta, num_sample_process);
for ii = 1:num_ant_sta
    fileID = fopen(['signals/cali_user_signal_rx' num2str(ii-1) '.dat'], 'rb');
    if (fileID < 0)
        disp('Error: fail to open files!');
        pause;
    end
    frewind(fileID);
    fseek(fileID, num_sample_shift, 'bof');
    data_ant_f = fread(fileID, 2*num_sample_process, 'float');
    data_ant_c = transpose(data_ant_f(1:2:end) + 1i*data_ant_f(2:2:end));
    signal_stream_original(ii,:) = data_ant_c(1:num_sample_process);
    fclose(fileID);
end
signal_stream_original =  signal_stream_original*1000;
signal_stream_combined = signal_stream_original;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time synchronization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% Cross correlation
num_preamble = length(PREAMABLE);
local_signal = PREAMABLE;
cross_correlation_arr = zeros(1,num_sample_process);
for ii = 100:3*num_ant_ap*80*num_symbol_per_frame+20000
    signal_segment = signal_stream_combined(ii:ii+num_preamble-1);
    cross_correlation_arr(ii) =  local_signal*signal_segment'/(norm(local_signal)*norm(signal_segment));
end
[max_cc_value, max_cc_pos] = max(abs(cross_correlation_arr(1:10000)));

% figure; hold on; 
% plot(abs(cross_correlation_arr),'k*-');
% xlabel('sample index');
% ylabel('correlation');
% title('Cross correlation (DL)');
% box on;grid on; 
% axis([0 length(signal_stream_combined) 0 1]);

 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove redundant samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
begofframe = max_cc_pos + num_preamble - num_move_sample;
preamable_signal = signal_stream_original(begofframe-num_preamble:begofframe-1);
signal_stream_buf = signal_stream_original(:,begofframe:end);
signal_stream_combined = signal_stream_combined(:,begofframe:end);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% frequency synchronization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
preamable_signal1 = preamable_signal(1:num_preamble/2);
preamable_signal2 = preamable_signal(num_preamble/2+1:num_preamble);
preamable_signal1 = preamable_signal1(2:end-1);
preamable_signal2 = preamable_signal2(2:end-1);
phase_offset = preamable_signal1*preamable_signal2';
angle_delta = angle(phase_offset)/length(preamable_signal1);

% angle_delta = -0.00071893;

freq_comp = exp(1i*angle_delta*(1:length(signal_stream_combined)));
freq_comp_fram = repmat(freq_comp, size(signal_stream_buf,1), 1);
signal_stream_buf = signal_stream_buf.*freq_comp_fram;

%disp(['downlink freq offset: ' num2str(angle_delta)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% channel estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_pkt = 20;
pkt_len = 1800;
chan_coeff_matrix = zeros(num_ant_ap, num_valid_sc, num_pkt);
for pkt_idx = 1:num_pkt
    % re-orgainze the signal 
    time_signal_arr = signal_stream_buf((pkt_idx-1)*pkt_len+1:(pkt_idx-1)*pkt_len+num_ant_ap*80*3);
    time_signal_frame = reshape(time_signal_arr, 240,[]);
    time_signal_frame = time_signal_frame(2*80+16+1:2*80+16+64, :)+time_signal_frame(1*80+16+1:1*80+16+64, :);
    time_signal_frame = time_signal_frame/2;
    freq_signal_frame = fft(time_signal_frame)/sqrt(num_subcarrier);
    freq_signal_frame = freq_signal_frame(POS_VALID_SC,:);
    % channel estimation
    for ant_idx = 1:num_ant_ap
        x_arr = transpose(LTF_DAT_GRP(ant_idx,:));   
        y_arr = freq_signal_frame(:, ant_idx);
        h_arr = y_arr./x_arr;
        chan_coeff_matrix(ant_idx,:,pkt_idx) = h_arr;
    end
end
 

% for ant_idx = 1:num_ant_ap
%     figure; 
%     subplot(2, 1, 1);
%     title(['DL channel Amp - Ant ' num2str(173+ant_idx)]);    
%     hold on;
%     grid on;    
%     for sym_idx = 1:num_pkt
%         plot(abs(chan_coeff_matrix(ant_idx,:,sym_idx)),'*-');
%     end
%     subplot(2, 1, 2);
%     title(['DL channel Pha - Ant ' num2str(173+ant_idx)]);    
%     hold on;
%     grid on;    
%     for sym_idx = 1:num_pkt
%         plot(unwrap(angle((chan_coeff_matrix(ant_idx,:,sym_idx)))),'*-');
%     end    
% end
dl_chan_coeff_matrix = chan_coeff_matrix;
dl_chan=dl_chan_coeff_matrix;
%save dl_chan.mat dl_chan
 

%============================================================================================
% calculate the compensation coeff
%============================================================================================
dl_chan = mean(dl_chan_coeff_matrix,3);
ul_chan = mean(ul_chan_coeff_matrix,3);
% dl_chan = sum(dl_chan_coeff_matrix(:,:,1:2),3);
% ul_chan = sum(ul_chan_coeff_matrix(:,:,1:2),3);

dl_chan_ratio = dl_chan./repmat(dl_chan(1,:), num_ant_ap, 1);
ul_chan_ratio = ul_chan./repmat(ul_chan(1,:), num_ant_ap, 1);

calib_coeff = ul_chan_ratio./dl_chan_ratio;

for rx_idx = 1:num_ant_ap
    figure; 
    subplot(2,1,1);
    hold on; box on;grid on; 
    title('channel calibration: amplitude');
    axis([0 52 -2 2]);
    plot(abs(calib_coeff(rx_idx,:)),'b*-');
    subplot(2,1,2);
    hold on; box on;grid on; 
    title('channel calibration: phase');
    axis([0 52 -pi pi]);
    plot(unwrap(angle(calib_coeff(rx_idx,:))),'r*-');
    calib_coeff(rx_idx) = mean(mean(calib_coeff(rx_idx,:,:)));
    disp(['antenna ' num2str(rx_idx) ': ' num2str(calib_coeff(rx_idx)) '; ampl=' num2str(abs(calib_coeff(rx_idx))) '; angle=' num2str(unwrap(angle(calib_coeff(rx_idx))))]);
end

chan_calib_coeff = mean(calib_coeff, 2);
c_calib=chan_calib_coeff(1:4);
save('c_calib.mat','c_calib')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%        D2D %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_ant_ap = 3;
num_valid_sc = 52;
num_subcarrier = 64;
num_symbol_per_frame = 20;
num_frame = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LTF
LTF_DAT_GRP = zeros(10, 52);
LTF_DAT_GRP(1,:) = [-1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1 1 1 -1 -1 -1];
LTF_DAT_GRP(2,:) = [-1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1 1 1 -1 -1 -1];
LTF_DAT_GRP(3,:) = [-1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1 1 1 -1 -1 -1];
LTF_DAT_GRP(4,:) = [-1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1 1 1 -1 -1 -1];
LTF_DAT_GRP(5,:) = [-1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1 1 1 -1 -1 -1];
LTF_DAT_GRP(6,:) = [-1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1 1 1 -1 -1 -1];
LTF_DAT_GRP(7,:) = [-1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1 1 1 -1 -1 -1];
LTF_DAT_GRP(8,:) = [-1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1 1 1 -1 -1 -1];
LTF_DAT_GRP(9,:) = [-1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1 1 1 -1 -1 -1];


%%%% pilot
POS_VALID_SC = [39:64 2:27];
POS_PAYLOAD_SC = [1:5 7:19 21:32 34:46 48:52];
POS_PILOT_SC = [6 20 33 47];
PILOT_DAT_GRP = zeros(4,num_symbol_per_frame);
PILOT_DAT_GRP(1,:) = [-1 -1 -1 -1 1 1 1 -1 1 -1 -1 1 1 1 1 -1 -1 1 -1 1];
PILOT_DAT_GRP(2,:) = [1 -1 -1 1 -1 1 1 1 1 1 1 -1 1 1 1 -1 -1 1 1 1];
PILOT_DAT_GRP(3,:) = [-1 -1 1 1 -1 -1 -1 1 1 1 -1 -1 1 -1 1 1 1 1 1 1];
PILOT_DAT_GRP(4,:) = [1 1 1 -1 1 1 -1 -1 -1 -1 1 -1 -1 -1 -1 -1 -1 1 1 -1];

PREAMABLE = [1,1,-1,1,-1,-1,-1,-1,1,1,1,1,-1,1,-1,1,-1,1,-1,1,1,1,1,1,-1,1,-1,1,1,1,-1,-1,1,-1,1,1,1,1,1,1,1,-1,-1,-1,1,-1,1,-1,-1,-1,-1,1,-1,-1,-1,1,1,1,1,-1,-1,1,1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,-1,-1,-1,1,-1,1,1,1,-1,1,-1,-1,-1,-1,-1,-1,1,1,-1,1,-1,-1,1,-1,1];
PREAMABLE = [PREAMABLE PREAMABLE];
PREAMABLE = PREAMABLE/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate freq data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
time_data_matrix_all_ant = zeros(num_ant_ap, 80, num_symbol_per_frame);
for ant_idx = 1:num_ant_ap
    % generate freq data
    freq_data_matrix = zeros(num_valid_sc, num_symbol_per_frame);
    for sym_idx = 3*num_ant_ap+1:num_symbol_per_frame %2*num_ant_ap+1
        freq_data_matrix(:,sym_idx) = (1/2)*(sign(randn(num_valid_sc,1))+1i*sign(randn(num_valid_sc,1)));
    end
    freq_data_matrix(:,3*ant_idx-2) = transpose(LTF_DAT_GRP(ant_idx,:));
    freq_data_matrix(:,3*ant_idx-1) = transpose(LTF_DAT_GRP(ant_idx,:));
    freq_data_matrix(:,3*ant_idx-0) = transpose(LTF_DAT_GRP(ant_idx,:));
    % Convert to time domain
    freq_data_matrix_tmp = zeros(num_subcarrier, num_symbol_per_frame);
    freq_data_matrix_tmp(POS_VALID_SC,:) = freq_data_matrix;
    time_data_matrix = ifft(freq_data_matrix_tmp)*sqrt(num_subcarrier);
    time_data_matrix = [time_data_matrix(end-15:end,:); time_data_matrix];
    % % add preamble
    time_data_matrix_all_ant(ant_idx,:,:) = time_data_matrix;
end
 



    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build orthogonal tx signal for channel estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_data_matrix_grp = zeros(num_ant_ap, 80*num_symbol_per_frame);
for ant_idx = 1:num_ant_ap
    time_data_matrix_grp(ant_idx,:) = reshape(time_data_matrix_all_ant(ant_idx,:,:), 1, []);
end
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add preamble
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ampl_avg = mean(mean(mean(abs(time_data_matrix_grp))));
preamable_matrix = repmat(num_ant_ap*ampl_avg*PREAMABLE, num_ant_ap, 1);
tx_signal_in_time = [preamable_matrix time_data_matrix_grp];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% how much data to write
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_data_matrix = repmat(tx_signal_in_time, 1, 20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write to the file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:num_ant_ap
    %fileID= fopen(['/tmp/data/cali_ap_signal_tx' num2str(ii-1) '.dat'], 'wb');
    fileID= fopen(['signals/d2d_ap_signal_tx' num2str(ii-1) '.dat'], 'wb');
    if (fileID < 0)
        disp('Error: fail to open files!');
        pause;
    end
    data_tmp_c = time_data_matrix(ii,:);
    data_tmp_f = [real(data_tmp_c); imag(data_tmp_c)];
    data_tmp_f = data_tmp_f(:);
    fwrite(fileID,data_tmp_f,'float');
    fclose(fileID);
end

disp('D2D communication')

%% Parameters
num_ant_sta = 1;
num_valid_sc = 52;
num_subcarrier = 64;
num_symbol_per_frame = 20;
num_frame = 20;
num_sample_process = 50000;
num_sample_shift = 2e6;
num_move_sample = -1;

%% STF 
STF_POS = [3  7  11  15  19  23  30  34  38  42  46  50];
STF_DAT_GRP = zeros(10, 52);    
STF_DAT_GRP(1,STF_POS) = sqrt(2)*[-1-1i  -1+1i  -1-1i  -1-1i  -1-1i  1-1i  -1-1i  -1-1i  1+1i  1-1i  1+1i  -1+1i];
STF_DAT_GRP(2,STF_POS) = sqrt(2)*[-1+1i  -1+1i  1-1i  1-1i  -1-1i  1-1i  -1+1i  1+1i  1-1i  -1+1i  1-1i  -1+1i];
STF_DAT_GRP(3,STF_POS) = sqrt(2)*[-1+1i  1+1i  -1+1i  1+1i  1-1i  1+1i  -1+1i  -1+1i  1-1i  1-1i  -1+1i  -1+1i];
STF_DAT_GRP(4,STF_POS) = sqrt(2)*[1+1i  -1+1i  -1+1i  -1+1i  -1-1i  1+1i  1+1i  1-1i  1-1i  -1+1i  -1+1i  -1-1i];
STF_DAT_GRP(5,STF_POS) = sqrt(2)*[-1+1i  -1+1i  1+1i  -1+1i  -1-1i  1+1i  1-1i  -1-1i  -1-1i  -1+1i  -1+1i  1-1i];
STF_DAT_GRP(6,STF_POS) = sqrt(2)*[1+1i  1-1i  -1-1i  1+1i  -1+1i  1+1i  1-1i  -1-1i  1+1i  -1-1i  -1+1i  -1+1i];
STF_DAT_GRP(7,STF_POS) = sqrt(2)*[1+1i  1+1i  1+1i  1-1i  -1-1i  1-1i  -1-1i  1-1i  1+1i  1+1i  -1-1i  1-1i];
STF_DAT_GRP(8,STF_POS) = sqrt(2)*[-1+1i  1-1i  1+1i  -1+1i  1+1i  1-1i  -1-1i  -1-1i  -1+1i  1-1i  -1+1i  -1-1i];
STF_DAT_GRP(9,STF_POS) = sqrt(2)*[1-1i  -1+1i  1-1i  -1+1i  1-1i  1+1i  1+1i  -1-1i  -1-1i  -1+1i  1+1i  -1-1i];
STF_DAT_GRP(10,STF_POS) = sqrt(2)*[1-1i  1-1i  -1-1i  1+1i  1+1i  1+1i  1+1i  1-1i  -1-1i  -1+1i  -1+1i  -1-1i];  
%% LTF
LTF_DAT_GRP = zeros(10, 52);
LTF_DAT_GRP(1,:) = [-1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1 1 1 -1 -1 -1];
LTF_DAT_GRP(2,:) = [-1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1 1 1 -1 -1 -1];
LTF_DAT_GRP(3,:) = [-1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1 1 1 -1 -1 -1];
LTF_DAT_GRP(4,:) = [-1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1 1 1 -1 -1 -1];
LTF_DAT_GRP(5,:) = [-1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1 1 1 -1 -1 -1];
LTF_DAT_GRP(6,:) = [-1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1 1 1 -1 -1 -1];
LTF_DAT_GRP(7,:) = [-1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1 1 1 -1 -1 -1];
LTF_DAT_GRP(8,:) = [-1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1 1 1 -1 -1 -1];
LTF_DAT_GRP(9,:) = [-1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1 1 1 -1 -1 -1];
%% Pilot
POS_VALID_SC = [39:64 2:27];
POS_PAYLOAD_SC = [1:5 7:19 21:32 34:46 48:52];
POS_PILOT_SC = [6 20 33 47];
PILOT_DAT_GRP = zeros(4,num_symbol_per_frame);
PILOT_DAT_GRP(1,:) = [-1 -1 -1 -1 1 1 1 -1 1 -1 -1 1 1 1 1 -1 -1 1 -1 1];
PILOT_DAT_GRP(2,:) = [1 -1 -1 1 -1 1 1 1 1 1 1 -1 1 1 1 -1 -1 1 1 1];
PILOT_DAT_GRP(3,:) = [-1 -1 1 1 -1 -1 -1 1 1 1 -1 -1 1 -1 1 1 1 1 1 1];
PILOT_DAT_GRP(4,:) = [1 1 1 -1 1 1 -1 -1 -1 -1 1 -1 -1 -1 -1 -1 -1 1 1 -1];
%% Preamble
PREAMABLE = [1,1,-1,1,-1,-1,-1,-1,1,1,1,1,-1,1,-1,1,-1,1,-1,1,1,1,1,1,-1,1,-1,1,1,1,-1,-1,1,-1,1,1,1,1,1,1,1,-1,-1,-1,1,-1,1,-1,-1,-1,-1,1,-1,-1,-1,1,1,1,1,-1,-1,1,1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,-1,-1,-1,1,-1,1,1,1,-1,1,-1,-1,-1,-1,-1,-1,1,1,-1,1,-1,-1,1,-1,1];
PREAMABLE = [PREAMABLE PREAMABLE];
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ==================== UL Channels ==================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
signal_stream_original = zeros(num_ant_ap, num_sample_process);
for ii = 1:num_ant_ap
    fileID = fopen(['signals/d2d_ap_signal_rx' num2str(ii-1) '.dat'], 'rb');
    if (fileID < 0)
        error('Error: fail to open files!');
    end
    frewind(fileID);
    fseek(fileID, num_sample_shift, 'bof');
    data_ant_f = fread(fileID, 2*num_sample_process, 'float');
    data_ant_c = transpose(data_ant_f(1:2:end) + 1i*data_ant_f(2:2:end));
    signal_stream_original(ii,:) = data_ant_c(1:num_sample_process);
    fclose(fileID);
end
signal_stream_original =  signal_stream_original*1000;
signal_stream_combined = sum(signal_stream_original,1);

%% time synchronization
% Cross correlation
data_freq = zeros(1, num_subcarrier);
data_freq(POS_VALID_SC) = LTF_DAT_GRP(1,:);
data_time = ifft(data_freq)*sqrt(num_subcarrier);
ltf_signal = [data_time(end-15:end) data_time];
local_signal = [ltf_signal  ltf_signal];
cross_correlation_arr = zeros(1,num_sample_process);
for ii = 200:10000
    signal_segment = signal_stream_combined(ii:ii+160-1);
    cross_correlation_arr(ii) =  local_signal*signal_segment'/(norm(local_signal)*norm(signal_segment));
end
% figure; hold on; grid on; 
% plot(abs(cross_correlation_arr),'k*-');
% xlabel('sample index');
% ylabel('correlation');
% title('Cross correlation (UL)');
% box on;
% axis([0 length(signal_stream_combined) 0 1]);

%% Capture frames
[~, begofframe] = max(abs(cross_correlation_arr(1:1600)));
begofframe = begofframe(1) + 1600 - 160 - num_move_sample;
signal_stream_buf = signal_stream_original(:, begofframe:end);
signal_stream_combined = signal_stream_combined(:, begofframe:end);


%% frequency synchronization
freq_offset1 = 0;
freq_offset2 = 0;
for ii = 0:1600:length(signal_stream_combined)-1600
    % CP
    for jj = 1:20
        freq_offset1 = freq_offset1 + signal_stream_combined(ii+80*(jj-1)+1:ii+80*(jj-1)+16)*signal_stream_combined(ii+80*(jj-1)+65:ii+80*(jj-1)+80)';
    end
    % LTF
    freq_offset2 = freq_offset2 + signal_stream_combined(ii+177:ii+240)*signal_stream_combined(ii+257:ii+320)';
end
angle_os1 = angle(freq_offset1)/64;
angle_os2 = angle(freq_offset2)/80;
freq_comp = exp(1i*angle_os1*(1:length(signal_stream_combined)));
freq_comp_fram = repmat(freq_comp, size(signal_stream_buf,1), 1);
signal_stream_buf = signal_stream_buf.*freq_comp_fram;

%disp(['uplink freq offset: ' num2str(angle_os1) '  ' num2str(angle_os2)]);
angle_os_ul = angle_os1;
angle_os_dl = -angle_os_ul;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% channel estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chan_coeff_matrix = zeros(num_ant_ap, num_valid_sc, num_frame);    
for frm_idx = 1:num_frame
    % reoganize the signal
    signal_frame = signal_stream_buf(:, 1600*(frm_idx-1)+1:1600*frm_idx);
    signal_matrix = zeros(num_ant_ap, num_valid_sc, num_symbol_per_frame);
    for rx_idx = 1:num_ant_ap
        data_t1 = signal_frame(rx_idx,:);
        data_t2 = reshape(data_t1, num_subcarrier+16, num_symbol_per_frame);
        data_t3 = data_t2(17:end,:);
        data_f1 = fft(data_t3)/sqrt(num_subcarrier);
        data_f2 = data_f1(POS_VALID_SC,:);
        signal_matrix(rx_idx,:,:) = data_f2;
    end
    % estimate channel
    for rx_idx = 1:num_ant_ap
        signal_ltf_avg = (signal_matrix(rx_idx, :, 3) + signal_matrix(rx_idx, :, 4))/2;
        chan_coeff_matrix(rx_idx,:,frm_idx) = signal_ltf_avg./LTF_DAT_GRP(1, :);
    end
end

% for ant_idx = 1:num_ant_ap
%     figure; 
%     subplot(2, 1, 1);
%     title(['UL channel Amp - Ant ' num2str(173+ant_idx)]);    
%     hold on;
%     grid on;    
%     for sym_idx = 1:num_symbol_per_frame
%         plot(abs(chan_coeff_matrix(ant_idx,:,sym_idx)),'*-');
%     end
%     subplot(2, 1, 2);
%     title(['UL channel Pha - Ant ' num2str(173+ant_idx)]);    
%     hold on;
%     grid on;    
%     for sym_idx = 1:num_symbol_per_frame
%         plot(unwrap(angle((chan_coeff_matrix(ant_idx,:,sym_idx)))),'*-');
%     end    
% end
ul_chan_coeff_matrix = chan_coeff_matrix;
ul_chan=ul_chan_coeff_matrix;
%save ul_chan.mat ul_chan

%==================================================
%================= DL Channel =====================
%==================================================
signal_stream_original = zeros(num_ant_sta, num_sample_process);
for ii = 1:num_ant_sta
    fileID = fopen(['signals/d2d_user_signal_rx' num2str(ii-1) '.dat'], 'rb');
    if (fileID < 0)
        disp('Error: fail to open files!');
        pause;
    end
    frewind(fileID);
    fseek(fileID, num_sample_shift, 'bof');
    data_ant_f = fread(fileID, 2*num_sample_process, 'float');
    data_ant_c = transpose(data_ant_f(1:2:end) + 1i*data_ant_f(2:2:end));
    signal_stream_original(ii,:) = data_ant_c(1:num_sample_process);
    fclose(fileID);
end
signal_stream_original =  signal_stream_original*1000;
signal_stream_combined = signal_stream_original;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time synchronization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% Cross correlation
num_preamble = length(PREAMABLE);
local_signal = PREAMABLE;
cross_correlation_arr = zeros(1,num_sample_process);
for ii = 100:3*num_ant_ap*80*num_symbol_per_frame+20000
    signal_segment = signal_stream_combined(ii:ii+num_preamble-1);
    cross_correlation_arr(ii) =  local_signal*signal_segment'/(norm(local_signal)*norm(signal_segment));
end
[max_cc_value, max_cc_pos] = max(abs(cross_correlation_arr(1:10000)));

% figure; hold on; 
% plot(abs(cross_correlation_arr),'k*-');
% xlabel('sample index');
% ylabel('correlation');
% title('Cross correlation (DL)');
% box on;grid on; 
% axis([0 length(signal_stream_combined) 0 1]);

 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove redundant samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
begofframe = max_cc_pos + num_preamble - num_move_sample;
preamable_signal = signal_stream_original(begofframe-num_preamble:begofframe-1);
signal_stream_buf = signal_stream_original(:,begofframe:end);
signal_stream_combined = signal_stream_combined(:,begofframe:end);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% frequency synchronization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
preamable_signal1 = preamable_signal(1:num_preamble/2);
preamable_signal2 = preamable_signal(num_preamble/2+1:num_preamble);
preamable_signal1 = preamable_signal1(2:end-1);
preamable_signal2 = preamable_signal2(2:end-1);
phase_offset = preamable_signal1*preamable_signal2';
angle_delta = angle(phase_offset)/length(preamable_signal1);

% angle_delta = -0.00071893;

freq_comp = exp(1i*angle_delta*(1:length(signal_stream_combined)));
freq_comp_fram = repmat(freq_comp, size(signal_stream_buf,1), 1);
signal_stream_buf = signal_stream_buf.*freq_comp_fram;

%disp(['downlink freq offset: ' num2str(angle_delta)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% channel estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_pkt = 20;
pkt_len = 1800;
chan_coeff_matrix = zeros(num_ant_ap, num_valid_sc, num_pkt);
for pkt_idx = 1:num_pkt
    % re-orgainze the signal 
    time_signal_arr = signal_stream_buf((pkt_idx-1)*pkt_len+1:(pkt_idx-1)*pkt_len+num_ant_ap*80*3);
    time_signal_frame = reshape(time_signal_arr, 240,[]);
    time_signal_frame = time_signal_frame(2*80+16+1:2*80+16+64, :)+time_signal_frame(1*80+16+1:1*80+16+64, :);
    time_signal_frame = time_signal_frame/2;
    freq_signal_frame = fft(time_signal_frame)/sqrt(num_subcarrier);
    freq_signal_frame = freq_signal_frame(POS_VALID_SC,:);
    % channel estimation
    for ant_idx = 1:num_ant_ap
        x_arr = transpose(LTF_DAT_GRP(ant_idx,:));   
        y_arr = freq_signal_frame(:, ant_idx);
        h_arr = y_arr./x_arr;
        chan_coeff_matrix(ant_idx,:,pkt_idx) = h_arr;
    end
end
 

% for ant_idx = 1:num_ant_ap
%     figure; 
%     subplot(2, 1, 1);
%     title(['DL channel Amp - Ant ' num2str(173+ant_idx)]);    
%     hold on;
%     grid on;    
%     for sym_idx = 1:num_pkt
%         plot(abs(chan_coeff_matrix(ant_idx,:,sym_idx)),'*-');
%     end
%     subplot(2, 1, 2);
%     title(['DL channel Pha - Ant ' num2str(173+ant_idx)]);    
%     hold on;
%     grid on;    
%     for sym_idx = 1:num_pkt
%         plot(unwrap(angle((chan_coeff_matrix(ant_idx,:,sym_idx)))),'*-');
%     end    
% end
dl_chan_coeff_matrix = chan_coeff_matrix;
dl_chan=dl_chan_coeff_matrix;
%save dl_chan.mat dl_chan
 

%============================================================================================
% calculate the compensation coeff
%============================================================================================
dl_chan = mean(dl_chan_coeff_matrix,3);
ul_chan = mean(ul_chan_coeff_matrix,3);
% dl_chan = sum(dl_chan_coeff_matrix(:,:,1:2),3);
% ul_chan = sum(ul_chan_coeff_matrix(:,:,1:2),3);

dl_chan_ratio = dl_chan./repmat(dl_chan(1,:), num_ant_ap, 1);
ul_chan_ratio = ul_chan./repmat(ul_chan(1,:), num_ant_ap, 1);

calib_coeff = ul_chan_ratio./dl_chan_ratio;
for rx_idx = 1:num_ant_ap
    figure; 
    subplot(2,1,1);
    hold on; box on;grid on; 
    title('D2D-channel calibration: amplitude');
    axis([0 52 -2 2]);
    plot(abs(calib_coeff(rx_idx,:)),'b*-');
    subplot(2,1,2);
    hold on; box on;grid on; 
    title('D2D-channel calibration: phase');
    axis([0 52 -pi pi]);
    plot(unwrap(angle(calib_coeff(rx_idx,:))),'r*-');
    calib_coeff(rx_idx) = mean(mean(calib_coeff(rx_idx,:,:)));
    disp(['antenna ' num2str(rx_idx) ': ' num2str(calib_coeff(rx_idx)) '; ampl=' num2str(abs(calib_coeff(rx_idx))) '; angle=' num2str(unwrap(angle(calib_coeff(rx_idx))))]);
end
chan_calib_coeff = mean(calib_coeff, 2);
d_calib=chan_calib_coeff(1:3);
save('d_calib.mat','d_calib')