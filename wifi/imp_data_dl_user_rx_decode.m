clear;
clc;
close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% system parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
is_debug_enabled = 1;
num_ant_ap = 1;
num_stm_ap = 2;
num_sample_processing = 16000;
num_sc = 64;
num_symbol_per_frame = 20;
num_sample_shift = 1e6;
num_averaging_sc = 4;
num_valid_sc = 52;
const_correl_threhold = 0.5;

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


USER_DECODED_SIGNAL_EVM = zeros(1,num_stm_ap);
for  user_index = 1:num_stm_ap
    
    data_stream_index = user_index;
    file_index = user_index;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    file_id = fopen(['/tmp/data_implicit/data_user_rx_signal' num2str(file_index-1) '.dat'], 'rb');%
    if (file_id < 0)
        error('Error: fail to open files!');
    end
    frewind(file_id);
    fseek(file_id, num_sample_shift, 'bof');
    rx_signal_f = fread(file_id, 2*num_sample_processing, 'float');
    rx_signal_buf = transpose(rx_signal_f(1:2:end) + 1i*rx_signal_f(2:2:end));
    fclose(file_id);
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % time synchronization
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    % Cross correlation
    ltf_freq_signal = zeros(1, num_sc);
    ltf_freq_signal(POS_VALID_SC) = LTF_DAT_GRP(user_index,:);
    lte_time_signal = ifft(ltf_freq_signal)*sqrt(num_sc);
    ltf_waveform = [lte_time_signal(end-15:end) lte_time_signal lte_time_signal(end-15:end) lte_time_signal];
    cross_correl_value = zeros(1,num_sample_processing);
    for ii = 161:num_sample_processing-160
        signal_segment = rx_signal_buf(ii:ii+160-1);
        cross_correl_value(ii) =  ltf_waveform*signal_segment'/(norm(ltf_waveform)*norm(signal_segment));
    end
    [max_correl_value, max_correl_pos] = max(abs(cross_correl_value(1:end-1600)));
    begofframe = max_correl_pos - 160;

    if (is_debug_enabled == 1 && max_correl_value > const_correl_threhold)
        figure(10*data_stream_index+0); hold on; grid on; box on;
        plot(abs(cross_correl_value),'k*-');
        xlabel('sample index');
        ylabel('correlation');
        axis([0 length(rx_signal_buf) 0 1]);
        title(['Cross correlation: ' num2str(data_stream_index)]);
    end
     
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find the begining of a frame
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rx_signal_buf = rx_signal_buf(begofframe:end);
    %rx_signal_buf = rx_signal_buf(begofframe-1:end); %????
    rx_signal_frame = rx_signal_buf(1:1600);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % frequency synchronization
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    freq_offset_cp_based = 0;
    freq_offset_ltf_based = 0;
    % CP
    for jj = 1:20
        freq_offset_cp_based = freq_offset_cp_based + rx_signal_frame(80*(jj-1)+1:80*(jj-1)+16)*rx_signal_frame(80*(jj-1)+65:80*(jj-1)+80)';
    end
    % LTF
    freq_offset_ltf_based = freq_offset_ltf_based + rx_signal_frame(177:240)*rx_signal_frame(257:320)';
    delta_phase_per_sample_cp_based = angle(freq_offset_cp_based)/64;
    delta_phase_per_sample_ltf_based = angle(freq_offset_ltf_based)/80;
    delta_phase_per_sample = delta_phase_per_sample_ltf_based;
    rx_dds = exp(1i*delta_phase_per_sample*(1:size(rx_signal_frame,2)));
    rx_dds_frame = repmat(rx_dds, size(rx_signal_frame,1), 1);
    rx_signal_frame = rx_signal_frame.*rx_dds_frame;

    if is_debug_enabled == 1
        disp(['cp-based freq offset estimate: ' num2str(delta_phase_per_sample_cp_based)]);
        disp(['ltf-based freq offset estimate: ' num2str(delta_phase_per_sample_ltf_based)]);
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % received frequency frame
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    time_signal_one_frame = rx_signal_frame;
    freq_signal_one_frame = zeros(num_ant_ap, num_valid_sc, num_symbol_per_frame);
    for rx_idx = 1:num_ant_ap
        time_signal_tmp1 = time_signal_one_frame(rx_idx,:);
        time_signal_tmp2 = reshape(time_signal_tmp1, num_sc+16, num_symbol_per_frame);
        time_signal_tmp3 = time_signal_tmp2(17:end,:);
        freq_signal_tmp1 = fft(time_signal_tmp3)/sqrt(num_sc);
        freq_signal_tmp2 = freq_signal_tmp1(POS_VALID_SC,:);
        freq_signal_one_frame(rx_idx,:,:) = freq_signal_tmp2;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute G
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute EYY and EYX
    EYY_SUM = zeros(num_ant_ap, num_ant_ap, num_valid_sc);
    EYX_SUM = zeros(num_ant_ap, num_valid_sc);
    for sc_idx = 1:num_valid_sc
        % TSF (OFDM 1)
        Y = freq_signal_one_frame(:, sc_idx, 1);
        EYY_SUM(:,:,sc_idx) = squeeze(EYY_SUM(:,:,sc_idx)) + Y*Y';
        X = STF_DAT_GRP(data_stream_index, sc_idx);
        EYX_SUM(:,sc_idx) = squeeze(EYX_SUM(:,sc_idx)) + Y*X';
        % TSF (OFDM 2)
        Y = freq_signal_one_frame(:, sc_idx, 2);
        EYY_SUM(:,:,sc_idx) = squeeze(EYY_SUM(:,:,sc_idx)) + Y*Y';
        X = STF_DAT_GRP(data_stream_index, sc_idx);
        EYX_SUM(:,sc_idx) = squeeze(EYX_SUM(:,sc_idx)) + Y*X';
        % LTF (OFDM 3)
        Y = freq_signal_one_frame(:, sc_idx, 3);
        EYY_SUM(:,:,sc_idx) = squeeze(EYY_SUM(:,:,sc_idx)) + Y*Y';
        X = LTF_DAT_GRP(data_stream_index, sc_idx);
        EYX_SUM(:,sc_idx) = squeeze(EYX_SUM(:,sc_idx)) + Y*X';
        % LTF (OFDM 4)
        Y = freq_signal_one_frame(:, sc_idx, 4);
        EYY_SUM(:,:,sc_idx) = squeeze(EYY_SUM(:,:,sc_idx)) + Y*Y';
        X = LTF_DAT_GRP(data_stream_index, sc_idx);
        EYX_SUM(:,sc_idx) = squeeze(EYX_SUM(:,sc_idx)) + Y*X';
        % Scattered pilots
        % if ismember(sc_idx, pos_pilot_sc)
        %     plt_idx = find(pos_pilot_sc == sc_idx);
        %     for sb_idx = 5:20
        %         y = freq_signal_one_frame(:, sc_idx, sb_idx);
        %         eyy_sum(:,:,sc_idx) = squeeze(eyy_sum(:,:,sc_idx)) + y*y';
        %         x = pilot_dat(plt_idx, sb_idx);
        %         eyx_sum(:,sc_idx) = squeeze(eyx_sum(:,sc_idx)) + y*x';
        %     end
        % end
    end
    % compute G
    G_ARR_PER_FRAME = zeros(num_ant_ap, num_valid_sc);
    for sc_idx = 1:num_valid_sc
        EYY = zeros(num_ant_ap, num_ant_ap);
        EYX = zeros(num_ant_ap, 1);
        sc_lwbd = max(1, sc_idx-num_averaging_sc);
        sc_upbd = min(sc_idx+num_averaging_sc, num_valid_sc);
        for kk = sc_lwbd:sc_upbd
            EYY = EYY + squeeze(EYY_SUM(:,:,kk));
            EYX = EYX + squeeze(EYX_SUM(:,kk));
        end
        G = pinv(EYY)*(EYX);
        G_ARR_PER_FRAME(:,sc_idx) = G;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % decoded data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    decoded_data_frame = zeros(num_valid_sc,num_symbol_per_frame);
    for sc_idx = 1:num_valid_sc
        G = G_ARR_PER_FRAME(:, sc_idx);
        Y = squeeze(freq_signal_one_frame(:, sc_idx, :));
        decoded_data_frame(sc_idx,:) = G'*Y;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % display channel
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if is_debug_enabled == 1
        h = figure(user_index*10+1); 
        hold on; grid on;box on;
        title(['G for stream ' num2str(data_stream_index)]);
        plot(real(G_ARR_PER_FRAME(1,:)),'ro-');
        plot(imag(G_ARR_PER_FRAME(1,:)),'bo-');
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % decoded data without phase compension
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    decoded_payload_wo_phase_comp = decoded_data_frame(POS_PAYLOAD_SC,5:num_symbol_per_frame); 
    decoded_payload_wo_phase_comp_arr = decoded_payload_wo_phase_comp(:);
    if is_debug_enabled == 11
        h = figure(user_index*10+2); 
        hold on; grid on; box on;
        axis([-2, 2, -2, 2]);
        title(['constellation of stream: ' num2str(data_stream_index)  ' (w/o phase compension)']);
        scatter(real(decoded_payload_wo_phase_comp_arr),imag(decoded_payload_wo_phase_comp_arr),'ko');
        scatter(real(decoded_payload_wo_phase_comp_arr),imag(decoded_payload_wo_phase_comp_arr),'k*');    
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % phase calibration for decoded data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    phase_offset1 = decoded_data_frame(POS_PILOT_SC(1),:)./PILOT_DAT_GRP(1,:);
    phase_offset2 = decoded_data_frame(POS_PILOT_SC(2),:)./PILOT_DAT_GRP(2,:);
    phase_offset3 = decoded_data_frame(POS_PILOT_SC(3),:)./PILOT_DAT_GRP(3,:);
    phase_offset4 = decoded_data_frame(POS_PILOT_SC(4),:)./PILOT_DAT_GRP(4,:);
    decoded_data_frame(1:13,:) = decoded_data_frame(1:13,:)./repmat(phase_offset1,13,1);
    decoded_data_frame(14:27,:) = decoded_data_frame(14:27,:)./repmat(phase_offset2,14,1);
    decoded_data_frame(28:40,:) = decoded_data_frame(28:40,:)./repmat(phase_offset3,13,1);
    decoded_data_frame(41:52,:) = decoded_data_frame(41:52,:)./repmat(phase_offset4,12,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % decoded data with phase compension
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    decoded_payload_w_phase_comp = decoded_data_frame(POS_PAYLOAD_SC,5:num_symbol_per_frame);  
    decoded_payload_w_phase_comp_arr = decoded_payload_w_phase_comp(:);
    if is_debug_enabled == 1
        h = figure(user_index*10+3); 
        hold on; grid on;box on;
        title(['constellation of stream: ' num2str(data_stream_index) ' (w/ phase compension)']);            
        scatter(real(decoded_payload_w_phase_comp_arr),imag(decoded_payload_w_phase_comp_arr),'ko');
        scatter(real(decoded_payload_w_phase_comp_arr),imag(decoded_payload_w_phase_comp_arr),'k*');    
        axis([-2, 2, -2, 2])
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute SNR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    D_FRAME = abs(real(decoded_payload_w_phase_comp_arr)) + 1i*abs(imag(decoded_payload_w_phase_comp_arr));
    P_FRAME = abs(mean(D_FRAME))^2;
    E_FRAME = mean(abs(D_FRAME - mean(D_FRAME)).^2);
    decoded_signal_evm = 10*log10(P_FRAME/E_FRAME);

    USER_DECODED_SIGNAL_EVM(user_index) = decoded_signal_evm;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% store data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = clock;
s = system('mkdir implicit_evm');
fname = 'implicit_evm/evm_';
fname = strcat(fname, sprintf('%04d',a(1)));
fname = strcat(fname, sprintf('%02d',a(2)));
fname = strcat(fname, sprintf('%02d',a(3)));
fname = strcat(fname, '_');
fname = strcat(fname, sprintf('%02d',a(4)));
fname = strcat(fname, sprintf('%02d',a(5)));
fname = strcat(fname, sprintf('%02d',round(a(6))));
fname = strcat(fname, '.mat');
save(fname, 'USER_DECODED_SIGNAL_EVM');
disp(USER_DECODED_SIGNAL_EVM);


