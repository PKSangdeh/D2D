%% Simplified NR Rx code
clc
clear all
close all
%% Control Parameters
is_debug_enabled=1;
cfo_enabled=1;
phase_track_enabled=1;
phase_track_type=2;   % 1 for regular and 2 for polynomial fit
samp_shift=-0;
%% Frame Parameters
cp_length=38;
num_str=3;
Nz=300;
n_ant=3;                        % Number of BS's antenna
nrb=25;                         % Number of resource block
n_sf_per_f=2;                   % Number of slot per frame
n_sym_per_sf=14;                % Number of subcarrier per frame
n_sc_per_rb=12;                 % Number of subcarrier per resource block
n_valid_sc=300;                 % Number of valid subcarrier
nfft=512;                       % FFT and IFFT points
n_sym=28;                       % Number of OFDM symbols in one frame
valid_sc_indx=[363:512,2:151];  % Index of valid subcarriers
n_frame=10;                     % Number of NR frames
len_frame=15950;
G_stack=zeros(num_str,n_ant,n_valid_sc);
evm_cel=zeros(num_str,n_frame);
%% Proccessing Parameters
load dmrs_patterns
%dmrs_val=[1,-1/sqrt(2)+1j/sqrt(2),1/sqrt(2)-1j/sqrt(2),-1/sqrt(2)-1j/sqrt(2)];
%% Preamble
load pn
%% PTRS
ptrs_cliche=zeros(n_sc_per_rb,n_sym_per_sf);
ptrs_val=1;
ptrs_cliche(end,:)=ptrs_val;
ptrs_pad=repmat(ptrs_cliche,nrb,n_sf_per_f);              % pattern of dmrs during whole resource grid
dmrs_pat_indx=1:num_str;
begf=zeros(1,num_str);
for st=3:3
    num_sample_processing = 150000;
    const_num_sample_shift = 10^5;
    const_correl_threhold = 0.25;
    num_frame_for_detection=1;
    input_sync_filter = [];
    rx_signal_buf = zeros(n_ant, num_sample_processing);
    for ii = 1:n_ant
        file_id = fopen(['signals/d2d_rx_' num2str(ii-1) '.dat'], 'rb');
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
    rx_signal_buf = rx_signal_buf*1;
    [U, ~] = eig(rx_signal_buf*rx_signal_buf');
    cross_correl_result = zeros(2, n_ant);
    for eigvec_idx = 1:n_ant
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % calculate eigen vector
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sync_filter = U(:, eigvec_idx);
        rx_combined_signal = sync_filter'*rx_signal_buf;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % time synchronization
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cross_correl_value1 = zeros(1,num_sample_processing);
        for ii = 1:num_sample_processing-length(pn(st,:))
            signal_segment = rx_combined_signal(ii:ii+length(pn(st,:))-1);
            cross_correl_value1(ii) =pn(st,:)*signal_segment'/(norm(pn(st,:))*norm(signal_segment));
        end
        cross_correl_value1 = abs(cross_correl_value1);
        cross_correl_value = cross_correl_value1;
        [max_correl_value, max_correl_pos] = max(abs(cross_correl_value));
        cross_correl_result(1, eigvec_idx) = max_correl_value;
        cross_correl_result(2, eigvec_idx) = max_correl_pos - samp_shift;
        %     if ((is_debug_enabled == 0||is_debug_enabled == 1||is_debug_enabled == 11) && cross_correl_result(1, eigvec_idx) > const_correl_threhold)
        %         %figure(10*data_stream_index+0); hold on; grid on; box on;
        %         figure
        %         stem(abs(cross_correl_value1),'b-');
        %         xlabel('sample index');
        %         ylabel('correlation');
        %         axis([0 length(rx_combined_signal) 0 1]);
        %         title(['Cross correlation-Eigenvector No.',num2str(eigvec_idx)]);
        %     end
    end
    [~,loc]=max(cross_correl_result(1,:));
    beg_of_frame=cross_correl_result(2,loc);
    synch_filter=U(:,loc);
    begf(st)=beg_of_frame;
    %% fecth one frame
    num_sample_processing=2.5*num_sample_processing;
    rx_signal_buf = zeros(n_ant, num_sample_processing);
    for ii = 1:n_ant
        file_id = fopen(['signals/d2d_rx_' num2str(ii-1) '.dat'], 'rb');
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
    for frm=1:n_frame
        % num_sample_shift = const_num_sample_shift+8*(begin_of_frame+len_frame*(frm_idx-1));
        rx_signal_frame = rx_signal_buf(:,(frm-1)*len_frame+beg_of_frame:beg_of_frame+(frm-1)*len_frame+15950-1);
        rx_combined_signal = sync_filter'*rx_signal_frame;
        if cfo_enabled==1
            freq_offset_cp_based = 0;
            % CP
            for jj = 1:29
                freq_offset_cp_based = freq_offset_cp_based + ...
                    rx_combined_signal(550*(jj-1)+1:550*(jj-1)+38)*rx_combined_signal(550*(jj-1)+513:550*jj)';
            end
            delta_phase_per_sample_cp_based = angle(freq_offset_cp_based)/nfft;
            delta_phase_per_sample = delta_phase_per_sample_cp_based;
            
            %delta_phase_per_sample = delta_phase_per_sample*0.85;
            
            rx_dds = exp(1i*delta_phase_per_sample*(1:size(rx_signal_frame,2)));
            
            
            
            
            rx_dds_frame = repmat(rx_dds,  size(rx_signal_frame,1), 1);
            rx_signal_frame = rx_signal_frame.*rx_dds_frame;
            disp(['cp-based freq offset estimate: ' num2str(delta_phase_per_sample_cp_based)]);
        end
        
        % fetch freq-domain frame
        time_signal_one_frame = rx_signal_frame;
        freq_signal_one_frame = zeros(n_ant, 300, n_sf_per_f*n_sym_per_sf);
        for rx_idx = 1:n_ant
            time_signal_tmp1 = time_signal_one_frame(rx_idx,:);
            time_signal_tmp2 = reshape(time_signal_tmp1,(nfft+cp_length), []);
            time_signal_tmp3 = time_signal_tmp2((cp_length+1):end,:);
            freq_signal_tmp1 = fft(time_signal_tmp3)/sqrt(nfft);
            freq_signal_tmp2 = freq_signal_tmp1(valid_sc_indx,2:end);
            freq_signal_one_frame(rx_idx,:,:) = freq_signal_tmp2;
        end
        rx_freq=freq_signal_one_frame;
        %
        % G FILTER Design
        load dmrs_patterns
        %ref_sig_val=dmrs_val(st);
        indx=3;
        cliche=squeeze(dmrs_patterns(indx,:,:));% dmrs pattern in one time slot
        pad=repmat(cliche,nrb,n_sf_per_f);              % pattern of dmrs during whole resource grid
        data_pad=ones(size(pad))-pad;
        f_wind=1; % other options: 5, 25
        t_wind=1; % other options: 2,5,10;
        n_f_wind=nrb/f_wind;
        n_t_wind=n_sf_per_f/t_wind;
        len_f=f_wind*n_sc_per_rb;
        len_t=t_wind*n_sym_per_sf;
        G_temp=zeros(n_ant,n_f_wind);
        decoded_frame=zeros(size(pad));
        aaa=[1 -1 1 -1 1 -1 -1 1 -1 -1 1]/2;
        for tt=1:n_t_wind
            for ff=1:n_f_wind
                % creating window pad
                win_pad=zeros(size(pad));
                tmp_pad=ones(len_f,len_t);
                win_pad(((ff-1)*len_f+1):ff*len_f,((tt-1)*len_t+1):tt*len_t)=tmp_pad;
                d_pad=data_pad.*win_pad;
                rs_pad=pad.*win_pad;
                eq_pad=d_pad+rs_pad;
                all_indx=find(eq_pad);
                rs_indx=find(rs_pad);
                d_indx=find(d_pad);
                X=pad(rs_indx);
                Y=zeros(n_ant,length(rs_indx));
                window=zeros(n_ant,length(all_indx));
                recovered=zeros(1,length(all_indx));
                for i=1:n_ant
                    tmp=squeeze(rx_freq(i,:,:));
                    Y(i,:)=tmp(rs_indx);
                    window(i,:)=tmp(all_indx);
                end
                EYY=0;
                EYX=0;
                X=aaa;
                for i=1:length(rs_indx)
                    EYY=EYY+Y(:,i)*Y(:,i)';
                    EYX=EYX+Y(:,i)*X(i);
                end
                G=pinv(EYY)*EYX;
                %tmp1=G/norm(G);
                tmp1=G;
                for ic=1:n_ant
                    G_stack(st,ic,(ff-1)*len_f+1:ff*len_f)=tmp1(ic);
                end
                for i=1:length(all_indx)
                    decoded_frame(all_indx(i))=G'*window(:,i);
                end
                
            end
        end
%         figure
%         subplot(2,1,1)
%         plot(abs(squeeze(G_stack(st,1,:))),'-r','linewidth',1)
%         hold on
%         plot(abs(squeeze(G_stack(st,2,:))),'-g','linewidth',1)
%         plot(abs(squeeze(G_stack(st,3,:))),'-m','linewidth',1)
%         plot(abs(squeeze(G_stack(st,4,:))),'-y','linewidth',1)
%         title('magnitude')
%         subplot(2,1,2)
%         plot(unwrap(angle(squeeze(G_stack(st,1,:)))),'-r');
%         hold on
%         plot(unwrap(angle(squeeze(G_stack(st,2,:)))),'-g');
%         plot(unwrap(angle(squeeze(G_stack(st,3,:)))),'-m');
%         plot(unwrap(angle(squeeze(G_stack(st,4,:)))),'-y');
%         title('phase')
        %% Phase offset correction
        len=n_sf_per_f*n_sym_per_sf;
        ang=zeros(nrb*n_sc_per_rb,len);
        if phase_track_enabled==1
            if phase_track_type==1
                for irb=1:nrb
                    tmp3=decoded_frame(irb*n_sc_per_rb,:)/ptrs_val;
                    ang((irb-1)*n_sc_per_rb+1:irb*n_sc_per_rb,:)=repmat(tmp3,n_sc_per_rb,1);
                end
                decoded_frame=decoded_frame./ang;
            elseif phase_track_type==2
                ang=zeros(nrb*n_sc_per_rb,len);
                phase_track=angle(decoded_frame(12:12:300,:));
                for irb=1:nrb
                    [a,~]=polyfit(1:len,phase_track(irb,:),1);
                    tmp=repmat([1:len],n_sc_per_rb,1);
                    tmp=exp(-1j*a(1)*tmp);
                    ang((irb-1)*n_sc_per_rb+1:irb*n_sc_per_rb,:)=tmp;
                end
                decoded_frame=decoded_frame.*ang;
                for irb=1:25
                    tmp3=decoded_frame(irb*n_sc_per_rb,:)/ptrs_val;
                    tmp4=mean(angle(tmp3));
                    tmp5=mean(abs(tmp3));
                    decoded_frame((irb-1)*n_sc_per_rb+1:irb*n_sc_per_rb,:)=decoded_frame((irb-1)*n_sc_per_rb+1:irb*n_sc_per_rb,:)*exp(-1j*tmp4)/tmp5;
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        data_pad=ones(size(pad))-pad-ptrs_pad;
        data=decoded_frame(find(data_pad));
        tmp=abs(real(data))+1i*abs(imag(data));
        tmp2=abs(mean(tmp))^2;
        tmp3=mean(abs(tmp - mean(tmp)).^2);
        output_evm = 10*log10(tmp2/tmp3);
        figure
        plot(real(data),imag(data),'o','MarkerFaceColor',[0.11,0.31,0.21],'MarkerEdgeColor','k')
        axis([-1.5 1.5 -1.5 1.5])
        title(['User ',num2str(st),':  EVM=-',num2str(output_evm),'dB'])
        evm_cel(st,frm)=output_evm;
    end
end
%%%%%%%%%%%%%%%%%%%%%% Interference resilience examination %%%%%%%%%%%%%%%
save('G_stack.mat','G_stack')
% diff=mod(abs(begf(1)-begf(2)),512)
evm=evm_cel;
%% EVM
for i=3:3
    disp(['EVM of D2D USER ',num2str(i)])
    evm_cel(i,:)
end
