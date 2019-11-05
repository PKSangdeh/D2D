function [is_frame_detected, output_evm, output_sync_filter] = ...
         imp_sounding_ul_ap_decode_one_stream(num_ant_ap, data_stream_index, num_frame_for_detection, input_sync_filter, is_debug_enabled,M);
    global evm
    output_evm = -100;
    output_sync_filter = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % system parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    const_num_sc = 64;
    num_valid_sc = 52;    
    const_num_sample_shift = 2e6;
    const_correl_threhold = 0.5;
    LTF_DAT_GRP = zeros(7, 52);
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
    POS_VALID_SC = [39:64 2:27];
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % fetch received signals
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    num_sample_processing = 2000;    
    rx_signal_buf = zeros(num_ant_ap, num_sample_processing);
    for ii = 1:num_ant_ap
%         file_id = fopen(['/tmp/data_implicit/sounding_ap_rx_signal' num2str(ii-1) '.dat'], 'rb');
          file_id = fopen(['signals/data_rx_signal' num2str(ii-1) '.dat'], 'rb');
        if (file_id < 0)
            error('Error: fail to open files!');
        end
        frewind(file_id);
        fseek(file_id, const_num_sample_shift, 'bof');
        rx_signal_f = fread(file_id, 2*num_sample_processing, 'float');
        rx_signal_c = transpose(rx_signal_f(1:2:end) + 1i*rx_signal_f(2:2:end));
        rx_signal_buf(ii,:) = rx_signal_c;
        fclose(file_id);
    end
    rx_signal_buf = rx_signal_buf*1000;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % frame detection of received signals
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [U, ~] = eig(rx_signal_buf*rx_signal_buf');
    cross_correl_result = zeros(2, num_ant_ap);
    for eigvec_idx = 1:num_ant_ap
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % calculate eigen vector
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isempty(input_sync_filter) == true
            sync_filter = U(:, eigvec_idx);
        else
            sync_filter = input_sync_filter;
        end
        rx_combined_signal = sync_filter'*rx_signal_buf;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % time synchronization
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        % Auto correlation
        if is_debug_enabled == 11
            auto_correl_value = zeros(1,num_sample_processing);
            for ii = 1:num_sample_processing-16
                sync_sig = rx_combined_signal(ii:ii+15);
                signal_segment = rx_combined_signal(ii+16:ii+16+15);
                auto_correl_value(ii) =  sync_sig*signal_segment'/(norm(sync_sig)*norm(signal_segment));
            end
            figure; hold on; grid on; 
            plot(abs(auto_correl_value),'k*-');
            xlabel('sample index');
            ylabel('correlation');
            title('Auto correlation after IC');
            box on;
            axis([0 length(rx_combined_signal) 0 1]);
        end

        % Cross correlation
        ltf_freq_signal = zeros(1, const_num_sc);
        ltf_freq_signal(POS_VALID_SC) = LTF_DAT_GRP(data_stream_index,:);
        lte_time_signal = ifft(ltf_freq_signal)*sqrt(const_num_sc);
        ltf_waveform = [lte_time_signal(end-15:end) lte_time_signal lte_time_signal(end-15:end) lte_time_signal];
        cross_correl_value1 = zeros(1,num_sample_processing);
        for ii = 3:num_sample_processing-160
            signal_segment = rx_combined_signal(ii:ii+160-1);
            cross_correl_value1(ii) = ltf_waveform*signal_segment'/(norm(ltf_waveform)*norm(signal_segment));
        end
        cross_correl_value1 = abs(cross_correl_value1);
        cross_correl_value2 = [0 cross_correl_value1(1:end-1)];
        cross_correl_value = cross_correl_value1 + cross_correl_value2;
        %cross_correl_value = cross_correl_value1;
        
        [max_correl_value, max_correl_pos] = max(abs(cross_correl_value));
        cross_correl_result(1, eigvec_idx) = max_correl_value;
        if (cross_correl_value1(max_correl_pos - 1) > 0.25)
            cross_correl_result(2, eigvec_idx) = max_correl_pos - 160 - 2;
        else
            cross_correl_result(2, eigvec_idx) = max_correl_pos - 160 - 1;
        end
        
        if (is_debug_enabled == 11 && max_correl_value > const_correl_threhold)
            figure(10*data_stream_index+0); hold on; grid on; box on;
            plot(abs(cross_correl_value),'k*-');
            plot(abs(cross_correl_value1),'ro-');
            xlabel('sample index');
            ylabel('correlation');
            axis([0 length(rx_combined_signal) 0 1]);
            title(['Cross correlation: ' num2str(data_stream_index)]);
        end
        if isempty(input_sync_filter) == false
            break;
        end        
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % is frame detected?
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [max_correl_value, pos_tmp] = max(cross_correl_result(1,:));
    begin_of_frame = cross_correl_result(2,pos_tmp);
    if max_correl_value > const_correl_threhold
        is_frame_detected = true;
        if isempty(input_sync_filter) == true
            sync_filter = U(:, pos_tmp);
        else
            sync_filter = input_sync_filter;
        end
    else
        is_frame_detected = false;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % detect signals in frame
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if max_correl_value > const_correl_threhold
        decode_filter_mtx = zeros(num_ant_ap, num_valid_sc);
        for frm_idx = 1:num_frame_for_detection
            % fetch one frame of data from file
            num_sample_processing = 1600; 
            num_sample_shift = const_num_sample_shift+8*(begin_of_frame+1600*(frm_idx-1));
            rx_signal_buf = zeros(num_ant_ap, num_sample_processing);
            for ii = 1:num_ant_ap
                file_id = fopen(['signals/data_rx_signal' num2str(ii-1) '.dat'], 'rb');
                if (file_id < 0)
                    error('Error: fail to open files!');
                end
                frewind(file_id);
                fseek(file_id, num_sample_shift, 'bof');
                rx_signal_f = fread(file_id, 2*num_sample_processing, 'float');
                rx_signal_buf(ii,:) = transpose(rx_signal_f(1:2:end) + 1i*rx_signal_f(2:2:end));
                fclose(file_id);
            end
            rx_signal_frame =  rx_signal_buf*1000;

            % decoding one frame of signal
            [output_evm, output_sync_filter, decode_filter_arr] = ...
                imp_sounding_ul_ap_signal_detection_per_frame(num_ant_ap, rx_signal_frame, ...
                data_stream_index, sync_filter, is_debug_enabled,M);
            
            % process results
            if (frm_idx == 1)
                [~, ref_pos] = max(mean(abs(decode_filter_arr).^2,2));
                ref_pos = ref_pos(1);
            end
            decode_filter_arr_normized = decode_filter_arr./repmat(decode_filter_arr(ref_pos,:),size(decode_filter_arr,1),1);
            decode_filter_mtx = decode_filter_mtx + decode_filter_arr_normized;
            if is_debug_enabled == 11
                h = figure(10*data_stream_index+2); 
                title(['channel ratio for stream: ' num2str(data_stream_index)]);
                hold on; grid on;box on;
                plot(real(decode_filter_arr_normized(1,:)),'ro-');
                plot(imag(decode_filter_arr_normized(1,:)),'r*-');
                plot(real(decode_filter_arr_normized(2,:)),'bo-');
                plot(imag(decode_filter_arr_normized(2,:)),'b*-');
                plot(real(decode_filter_arr_normized(3,:)),'ko-');
                plot(imag(decode_filter_arr_normized(3,:)),'k*-');
                plot(real(decode_filter_arr_normized(4,:)),'mo-');
                plot(imag(decode_filter_arr_normized(4,:)),'m*-');
            end

        end
        
        % write precoding coefficient
        file_id = fopen(['decode_filter_' num2str(data_stream_index-1) '.dat'], 'wb');
        if (file_id < 0)
            disp('Error: fail to open files!');
            pause;
        end
        decode_filter = decode_filter_arr(:).';    
        decode_filter_tmp = [real(decode_filter); imag(decode_filter)];
        decode_filter_f = decode_filter_tmp(:);
        fwrite(file_id, decode_filter_f, 'float');
        fclose(file_id);
    end
    
