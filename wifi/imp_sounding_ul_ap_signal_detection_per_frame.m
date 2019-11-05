function [output_evm, output_sync_filter, decode_filter_arr] = ...
    imp_sounding_ul_ap_signal_detection_per_frame(num_ant_ap, rx_signal_frame, ...
    data_stream_index, sync_filter, is_debug_enabled,M)
    global evm
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % system parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    const_num_sc = 64;
    num_symbol_per_frame = 20;
    num_valid_sc = 52;
    const_num_averaging_sc = 2;    
    num_sc_for_filter_average = 2;   
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % frequency synchronization
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rx_combined_signal = sync_filter'*rx_signal_frame;
    freq_offset_cp_based = 0;
    freq_offset_ltf_based = 0;
    % CP
    for jj = 1:20
        freq_offset_cp_based = freq_offset_cp_based + ...
                rx_combined_signal(80*(jj-1)+1:80*(jj-1)+16)*rx_combined_signal(80*(jj-1)+65:80*(jj-1)+80)';
    end
    % LTF
    freq_offset_ltf_based = freq_offset_ltf_based + rx_combined_signal(177:240)*rx_combined_signal(257:320)';
    delta_phase_per_sample_cp_based = angle(freq_offset_cp_based)/64;
    delta_phase_per_sample_ltf_based = angle(freq_offset_ltf_based)/80;
    delta_phase_per_sample = delta_phase_per_sample_cp_based;
    rx_dds = exp(1i*delta_phase_per_sample*([1:size(rx_signal_frame,2)]));
    rx_dds_frame = repmat(rx_dds,  size(rx_signal_frame,1), 1);
    rx_signal_frame = rx_signal_frame.*rx_dds_frame;

    if is_debug_enabled == 11
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
        time_signal_tmp2 = reshape(time_signal_tmp1, const_num_sc+16, num_symbol_per_frame);
        time_signal_tmp3 = time_signal_tmp2(17:end,:);
        freq_signal_tmp1 = fft(time_signal_tmp3)/sqrt(const_num_sc);
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
        % if ismember(sc_idx, POS_PILOT_SC)
        %     plt_idx = find(POS_PILOT_SC == sc_idx);
        %     for sb_idx = 5:20
        %         Y = freq_signal_one_frame(:, sc_idx, sb_idx);
        %         EYY_SUM(:,:,sc_idx) = squeeze(EYY_SUM(:,:,sc_idx)) + Y*Y';
        %         X = PILOT_DAT_GRP(plt_idx, sb_idx);
        %         EYX_SUM(:,sc_idx) = squeeze(EYX_SUM(:,sc_idx)) + Y*X';
        %     end
        % end
    end
    % compute G
    decode_filter_arr = zeros(num_ant_ap, num_valid_sc);
    for sc_idx = 1:num_valid_sc
        EYY = zeros(num_ant_ap, num_ant_ap);
        EYX = zeros(num_ant_ap, 1);
        sc_lwbd = max(1, sc_idx-const_num_averaging_sc);
        sc_upbd = min(sc_idx+const_num_averaging_sc, num_valid_sc);
        for kk = sc_lwbd:sc_upbd
            EYY = EYY + squeeze(EYY_SUM(:,:,kk));
            EYX = EYX + squeeze(EYX_SUM(:,kk));
        end
        G = pinv(EYY)*(EYX);
        decode_filter_arr(:,sc_idx) = G;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % average filter
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    decode_filter_arr_tmp = decode_filter_arr;
    decode_filter_arr = zeros(num_ant_ap, num_valid_sc);
    for sc_idx = 1:num_valid_sc
        sc_lwbd = max(1, sc_idx-num_sc_for_filter_average);
        sc_upbd = min(sc_idx+num_sc_for_filter_average, num_valid_sc);
        decode_filter_arr(:,sc_idx) = mean(decode_filter_arr_tmp(:,sc_lwbd:sc_upbd),2);
    end 
        
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % decoded data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    decoded_data_frame = zeros(num_valid_sc,num_symbol_per_frame);
    for sc_idx = 1:num_valid_sc
        G = decode_filter_arr(:, sc_idx);
        Y = squeeze(freq_signal_one_frame(:, sc_idx, :));
        decoded_data_frame(sc_idx,:) = G'*Y;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % display channel
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if is_debug_enabled == 11
        h = figure(10*data_stream_index+1);
        subplot(2,1,1);
        hold on; grid on;box on;
        title(['|G| for stream ' num2str(data_stream_index)]);
        plot(abs(decode_filter_arr(1,:)),'r-');
        subplot(2,1,2);
        hold on; grid on;box on;
        title(['|G| for stream ' num2str(data_stream_index)]);
        plot(abs(decode_filter_arr(2,:)),'r-');
%         subplot(2,2,3);
%         hold on; grid on;box on;
%         title(['|G| for stream ' num2str(data_stream_index)]);
%         plot(abs(decode_filter_arr(3,:)),'r-');
%         subplot(2,2,4);
%         hold on; grid on;box on;
%         title(['|G| for stream ' num2str(data_stream_index)]);
%         plot(abs(decode_filter_arr(4,:)),'r-');
    end
    if is_debug_enabled == 11
        h = figure(10*data_stream_index+2);
        subplot(2,1,1);
        hold on; grid on;box on;
        title(['angle(G) for stream ' num2str(data_stream_index)]);
        plot(unwrap(angle(decode_filter_arr(1,:))),'b-');
        subplot(2,1,2);
        hold on; grid on;box on;
        title(['angle(G) for stream ' num2str(data_stream_index)]);
        plot(unwrap(angle(decode_filter_arr(2,:))),'b-');
%         subplot(2,2,3);
%         hold on; grid on;box on;
%         title(['angle(G) for stream ' num2str(data_stream_index)]);
%         plot(unwrap(angle(decode_filter_arr(3,:))),'b-');
%         subplot(2,2,4);
%         hold on; grid on;box on;
%         title(['angle(G) for stream ' num2str(data_stream_index)]);
%         plot(unwrap(angle(decode_filter_arr(4,:))),'b-');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % decoded data without phase compension
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    decoded_payload_wo_phase_comp = decoded_data_frame(POS_PAYLOAD_SC,5:num_symbol_per_frame); 
    decoded_payload_wo_phase_comp_arr = decoded_payload_wo_phase_comp(:);
    if is_debug_enabled == 11
        h = figure(10*data_stream_index+3); 
        hold on; grid on; box on;
        axis([-2, 2, -2, 2]);
        %title(['constellation of stream: ' num2str(data_stream_index) ' (w/o phase compension)']);  
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
    phase_offset1(1:4) = 1;
    phase_offset2(1:4) = 1;
    phase_offset3(1:4) = 1;
    phase_offset4(1:4) = 1;
    decoded_data_frame(1:13,:) = decoded_data_frame(1:13,:)./repmat(phase_offset1,13,1);
    decoded_data_frame(14:27,:) = decoded_data_frame(14:27,:)./repmat(phase_offset2,14,1);
    decoded_data_frame(28:40,:) = decoded_data_frame(28:40,:)./repmat(phase_offset3,13,1);
    decoded_data_frame(41:52,:) = decoded_data_frame(41:52,:)./repmat(phase_offset4,12,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % decoded data with phase compension
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    decoded_payload_w_phase_comp = decoded_data_frame(POS_PAYLOAD_SC,5:num_symbol_per_frame);  
    % decoded_payload_w_phase_comp = decoded_data_frame(1:26,5:num_symbol_per_frame);  
    % decoded_payload_w_phase_comp = decoded_data_frame(27:52,5:num_symbol_per_frame);  
    decoded_payload_w_phase_comp_arr = decoded_payload_w_phase_comp(:);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute SNR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    output_evm = 0;
    if is_debug_enabled == 1        
        D_FRAME = abs(real(decoded_payload_w_phase_comp_arr)) + 1i*abs(imag(decoded_payload_w_phase_comp_arr));
        P_FRAME = abs(mean(D_FRAME))^2;
        E_FRAME = mean(abs(D_FRAME - mean(D_FRAME)).^2);
        output_evm = 10*log10(P_FRAME/E_FRAME);
    end
%fprintf('evm1=%f\n',output_evm)
    %% MY EVM
    rx_data=decoded_payload_w_phase_comp_arr;
    rx_decoded=qamdemod(rx_data,M,'UnitAveragePower',true);
    ref_points=qammod(rx_decoded,M,'UnitAveragePower',true);
    err_vec=abs(rx_data-ref_points);
    evm_vec=abs(ref_points).^2./err_vec.^2;
    evm_vec=10*log10(mean(evm_vec));
    output_evm1 =mean(evm_vec);
    %%
    D_FRAME1 = abs(real(decoded_payload_w_phase_comp_arr)) + 1i*abs(imag(decoded_payload_w_phase_comp_arr));
    P_FRAME1 = abs(mean(D_FRAME1))^2;
    E_FRAME1 = mean(abs(D_FRAME1 - mean(D_FRAME1)).^2);
    output_evm2 = 10*log10(P_FRAME1/E_FRAME1);
    %output_evm=min([output_evm1,output_evm2]);
    %output_evm=max([output_evm1,output_evm2]);
    output_evm=output_evm2;
    evm=[evm output_evm];
fprintf('EVM=-%2.2f dB\n',output_evm)
    if is_debug_enabled == 1
        h = figure(10*data_stream_index+4); 
        grid on;
        box on;
        %title(['constellation of stream: ' num2str(data_stream_index) ' (w/ phase compension)']);            
        scatter(real(decoded_payload_w_phase_comp_arr),imag(decoded_payload_w_phase_comp_arr),'o','MarkerFaceColor','[0.07 .21 0.14]','MarkerEdgeColor','k');    
        axis([-1.5, 1.5, -1.5, 1.5])
        title(['EVM=-',num2str(output_evm),' dB'])
        %pause(0.2)
    end

output_sync_filter = mean(decode_filter_arr(:, 20:30));

