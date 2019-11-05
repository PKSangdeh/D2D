clear;
clc;
close all;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_data_stream = 1;
num_valid_sc = 52;
num_subcarrier = 64;
num_symbol_per_frame = 20;
num_frame = 50;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STF 
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
%%% LTF
LTF_DAT_GRP = zeros(10, 52);
LTF_DAT_GRP(1,:) = [-1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1 1 1 -1 -1 -1];
LTF_DAT_GRP(2,:) = [1 1 -1 -1 1 1 -1 -1 1 1 -1 -1 -1 1 -1 1 1 1 -1 -1 -1 -1 1 1 1 -1 1 1 1 -1 -1 -1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 -1 -1 1 -1 -1 1 1 -1];
LTF_DAT_GRP(3,:) = [1 -1 -1 -1 -1 1 -1 -1 -1 1 1 -1 1 1 1 1 1 -1 -1 1 -1 -1 -1 -1 -1 1 -1 1 -1 -1 1 1 -1 -1 -1 1 1 -1 1 -1 1 -1 1 1 -1 -1 1 1 1 -1 1 -1];
LTF_DAT_GRP(4,:) = [1 1 1 -1 -1 1 -1 1 -1 -1 -1 1 1 1 -1 -1 1 1 -1 1 1 -1 -1 1 -1 -1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 -1 1 -1 -1 1 -1 1 -1 1 1];
LTF_DAT_GRP(5,:) = [1 -1 -1 1 1 1 1 1 -1 -1 -1 1 1 1 1 1 -1 -1 1 1 1 1 -1 -1 -1 1 1 1 -1 1 1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 -1 -1 -1 -1 1 -1 1 -1 -1];
LTF_DAT_GRP(6,:) = [-1 1 -1 1 -1 1 1 1 1 1 -1 1 1 1 -1 -1 -1 1 -1 1 1 1 1 1 1 1 -1 1 1 1 1 1 -1 1 -1 1 -1 -1 1 1 -1 -1 1 1 -1 -1 1 -1 -1 1 1 1];
LTF_DAT_GRP(7,:) = [-1 1 1 -1 1 -1 -1 1 1 -1 -1 -1 -1 -1 -1 1 -1 1 1 -1 1 1 -1 -1 1 -1 1 1 -1 1 1 -1 -1 1 -1 -1 -1 1 -1 -1 1 1 1 1 -1 -1 -1 -1 -1 1 1 1];
LTF_DAT_GRP(8,:) = [1 -1 1 1 1 -1 -1 -1 1 -1 1 1 -1 -1 -1 1 1 -1 1 1 1 -1 1 -1 -1 -1 1 -1 1 1 1 1 -1 -1 -1 1 1 1 -1 -1 1 1 -1 1 -1 1 -1 1 -1 1 1 1];
LTF_DAT_GRP(9,:) = [-1 1 1 -1 1 -1 -1 1 1 1 1 1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 -1 -1 1 -1 1 1 1 1 1 1 -1 1 -1 1 1 -1 -1 -1 1 -1 -1 1 1 1 -1 -1 -1 1 -1];
LTF_DAT_GRP(10,:) = [1 -1 -1 -1 -1 1 -1 -1 1 1 1 1 -1 -1 1 -1 -1 -1 1 -1 -1 -1 -1 1 1 1 1 1 1 1 1 -1 -1 -1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 -1 1 -1 -1 -1 -1];
%%%% pilot
POS_VALID_SC = [39:64 2:27];
POS_PAYLOAD_SC = [1:5 7:19 21:32 34:46 48:52];
POS_PILOT_SC = [6 20 33 47];
PILOT_DAT_GRP = zeros(4,num_symbol_per_frame);
PILOT_DAT_GRP(1,:) = [-1 -1 -1 -1 1 1 1 -1 1 -1 -1 1 1 1 1 -1 -1 1 -1 1];
PILOT_DAT_GRP(2,:) = [1 -1 -1 1 -1 1 1 1 1 1 1 -1 1 1 1 -1 -1 1 1 1];
PILOT_DAT_GRP(3,:) = [-1 -1 1 1 -1 -1 -1 1 1 1 -1 -1 1 -1 1 1 1 1 1 1];
PILOT_DAT_GRP(4,:) = [1 1 1 -1 1 1 -1 -1 -1 -1 1 -1 -1 -1 -1 -1 -1 1 1 -1];

   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate freq data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
freq_data_grp = zeros(num_data_stream, num_valid_sc, num_symbol_per_frame);
for ii = 1:num_data_stream
    stf_dat = STF_DAT_GRP(ii,:);
    ltf_dat = LTF_DAT_GRP(ii,:);
    data_freq = zeros(num_valid_sc,num_symbol_per_frame);
    data_freq(:,1) = stf_dat.';
    data_freq(:,2) = stf_dat.';
    data_freq(:,3) = ltf_dat.';
    data_freq(:,4) = ltf_dat.';
    for jj = 5:num_symbol_per_frame
        data_freq(:,jj) = (1/sqrt(2))*(sign(randn(num_valid_sc,1))+1i*sign(randn(num_valid_sc,1)));
        data_freq(POS_PILOT_SC,jj) = PILOT_DAT_GRP(:,jj);
    end
    freq_data_grp(ii,:,:) = data_freq;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert to time domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_signal_grp = zeros(num_data_stream, 80*num_symbol_per_frame);
for ii = 1:num_data_stream
    data_ant = squeeze(freq_data_grp(ii,:,:));
    data_f_frm = zeros(num_subcarrier, num_symbol_per_frame);
    data_f_frm(POS_VALID_SC,:) = data_ant;
    data_t_frm = ifft(data_f_frm)*sqrt(num_subcarrier);
    data_s_frm = [data_t_frm(end-15:end,:); data_t_frm];
    data_s_frm = data_s_frm(:).';
    time_signal_grp(ii,:) = data_s_frm;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Test Only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H = randn(num_data_stream, num_data_stream)+i*randn(num_data_stream, num_data_stream);
% H = eye(num_data_stream);
% time_signal_grp = H*time_signal_grp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% how much data to write
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_signal_grp = repmat(time_signal_grp, 1, num_frame);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write to the file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%system('mkdir /tmp/data_implicit');
for ii = 1:num_data_stream
    fileID= fopen(['signals/sounding_user_tx_signal' num2str(ii-1) '.dat'], 'wb');
    if (fileID < 0)
        error('Error: fail to open files!');
    end
    data_tmp_c = time_signal_grp(ii,:);
    data_tmp_f = [real(data_tmp_c); imag(data_tmp_c)];
    data_tmp_f = data_tmp_f(:);
    fwrite(fileID,data_tmp_f,'float');
    fclose(fileID);
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% verify the written data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V_DATA_S_GRP = zeros(num_data_stream, 80*num_symbol_per_frame);
time_signal_grp = time_signal_grp(:,1:80*num_symbol_per_frame);
for ii = 1:num_data_stream
    fileID= fopen(['signals/sounding_user_tx_signal' num2str(ii-1) '.dat'], 'rb');
    if (fileID < 0)
        disp('Error: fail to open files!');
        pause;
    end
    data_tmp_f = fread(fileID, 'float');
    data_tmp_c = transpose(data_tmp_f(1:2:end) + 1i*data_tmp_f(2:2:end));
    V_DATA_S_GRP(ii,:) = data_tmp_c(1:80*num_symbol_per_frame);
    fclose(fileID);
end
figure; plot(abs(V_DATA_S_GRP(:) - time_signal_grp(:)));


% get frequency data
V_ENC_DATA_FREQ_GRP = zeros(num_data_stream, num_valid_sc, num_symbol_per_frame);
for ii = 1:num_data_stream
    data_tmp = V_DATA_S_GRP(ii,:);
    data_tmp = reshape(data_tmp, 80, num_symbol_per_frame);
    data_tmp = data_tmp(17:end,:);
    data_tmp = fft(data_tmp)/sqrt(num_subcarrier);
    data_tmp = data_tmp(POS_VALID_SC,:);
    V_ENC_DATA_FREQ_GRP(ii,:,:) = data_tmp;
end
figure; plot(abs(V_ENC_DATA_FREQ_GRP(:)-freq_data_grp(:)));

% revert precoding
V_DATA_FREQ_GRP = V_ENC_DATA_FREQ_GRP(1:num_data_stream, :, :);
figure; plot(abs(V_DATA_FREQ_GRP(:)-freq_data_grp(:)));


 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check pilot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D_SUM = zeros(num_data_stream,num_symbol_per_frame);
for ii = 1:num_data_stream
    data_frm = squeeze(V_DATA_FREQ_GRP(ii,:,:));
    % stf
    stf_dat = data_frm(:,1);
    D_SUM(ii,1) = sum(abs(stf_dat - STF_DAT_GRP(ii,:).'));
    stf_dat = data_frm(:,2);
    D_SUM(ii,2) = sum(abs(stf_dat - STF_DAT_GRP(ii,:).'));
    % ltf
    ltf_dat = data_frm(:,3);
    D_SUM(ii,3) = sum(abs(ltf_dat - LTF_DAT_GRP(ii,:).'));
    ltf_dat = data_frm(:,3);
    D_SUM(ii,4) = sum(abs(ltf_dat - LTF_DAT_GRP(ii,:).'));
    % others
    for jj = 5:num_symbol_per_frame
        D_SUM(ii,jj) = sum(abs(data_frm(POS_PILOT_SC,jj) - PILOT_DAT_GRP(:,jj)));
    end
end

if sum(sum(D_SUM)) < 1e-3
    msg_info = (['Diff: ' num2str(sum(sum(D_SUM))) ': AP TX signals have been successfully generated.']);
else
    msg_info = (['Diff: ' num2str(sum(sum(D_SUM))) ': AP TX signal generation failed.']);
end
% msgbox(msg_info);


