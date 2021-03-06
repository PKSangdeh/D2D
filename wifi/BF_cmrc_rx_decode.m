clc;
close all
num_ant_ap = 3;
num_data_stream = 1;
num_frame_for_detection = 50;
M=4;
global evm
evm=[];
for data_stream_index = 1:num_data_stream
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % try singular vector filter
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    input_sync_filter = [];
    is_debug_enabled = 1;
    [is_frame_detectusered, output_evm, output_sync_filter] = ...
                imp_sounding_ul_ap_decode_one_stream(num_ant_ap, data_stream_index, ...
                num_frame_for_detection, input_sync_filter, is_debug_enabled,M);
    %disp(['output_evm[' num2str(data_stream_index) '] = ' num2str(output_evm)]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % try G
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    % input_sync_filter = output_sync_filter;
    % is_debug_enabled = 0;
    % num_frame_for_detection = 2;
    % [is_frame_detected, output_evm, output_sync_filter] = sounding_ul_ap_decode_one_stream(...
    %                             num_ant_ap, data_stream_index, num_frame_for_detection, ...
    %                             input_sync_filter, is_debug_enabled);

end

% is_debug_enabled = 0;
% imp_data_dl_ap_tx_encode(num_ant_ap, num_data_stream, is_debug_enabled);

a = 1;
%% Top five frames
evm=sort(evm);
if num_frame_for_detection>=5
    disp('====================================================')
fprintf('Top five evm: -%f -%f -%f -%f -%f ',evm(end), evm(end-1), evm(end-2),evm(end-3), evm(end-4))
end
