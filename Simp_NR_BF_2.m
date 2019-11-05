%% Simplified NR_TX for two cellular users in the UL
clear all
n_ant_c=4;
n_ant_d=3;
num_str=2;
wifi_enabled=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%% Cellular  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
cp_length=38;
Nz=300;
nrb=25;                         % Number of resource block
n_sf_per_f=2;                   % Number of slot per frame
n_sym_per_sf=14;                % Number of subcarrier per frame
n_sc_per_rb=12;                 % Number of subcarrier per resource block
n_valid_sc=300;                 % Number of valid subcarrier
nfft=512;                       % FFT and IFFT points
n_sym=28;                       % Number of OFDM symbols in one frame
valid_sc_indx=[363:512,2:151];  % Index of valid subcarriers
n_frame=100;                     % Number of NR frames
load c_calib
%c_calib=ones(1,n_ant_c);
%% Preamble
load pn
preamble=[pn(1,:).',pn(2,:).',pn(3,:).',pn(4,:).'];
%% DMRS patterns
load dmrs_patterns
load G_stack
%% PTRS pattern
ptrs_cliche=zeros(n_sc_per_rb,n_sym_per_sf);
ptrs_val=1;
ptrs_cliche(end,:)=ptrs_val;
ptrs_pad=repmat(ptrs_cliche,nrb,n_sf_per_f);              % pattern of dmrs during whole resource grid
freq_signal=zeros(n_ant_c,300,28);

for ii=1:num_str
    cliche=squeeze(dmrs_patterns(ii,:,:));
    dmrs_pad=repmat(cliche,nrb,n_sf_per_f);
    data_pad=ones(size(dmrs_pad))-dmrs_pad-ptrs_pad;
    cliche(find(cliche))=[1 -1 1 -1 1 -1 -1 1 -1 -1 1]/2;
    dmrs_pad=repmat(cliche,nrb,n_sf_per_f);
    data_freq=1/sqrt(2)*(2*(randi(2,size(dmrs_pad))-1)-1)+1j*1/sqrt(2)*(2*(randi(2,size(dmrs_pad))-1)-1);
    data_freq=data_freq.*data_pad;
    data_freq(find(dmrs_pad))=dmrs_pad(find(dmrs_pad));
    data_freq(find(ptrs_pad))=ptrs_val;
    for jj=1:n_ant_c
        v=squeeze(G_stack(ii,jj,:));
        v=conj(v);
        v=repmat(v,1,28);
        v=ones(size(v));
        freq_signal(jj,:,:)=squeeze(freq_signal(jj,:,:))+data_freq.*v;       
    end
end

%% write the beamformed signal
for ii=1:n_ant_c
    frame_f=zeros(512,28);
    frame_f(valid_sc_indx,:)=squeeze(freq_signal(ii,:,:));
    frame_t=ifft(frame_f)*sqrt(nfft);
    frame_t_c = [frame_t(end-(cp_length-1):end,:); frame_t];
    frame_t_c_2=[preamble(:,1), frame_t_c];
    sig_t=frame_t_c_2(:).';
    sig_t1=repmat(sig_t,1,100);
    sig_t1=sig_t1*c_calib(ii);
    % write file
    fileID= fopen(['signals/nr_bf_tx_' num2str(ii-1) '.dat'], 'wb');
    if (fileID < 0)
        error('Error: fail to open files!');
    end
    data_tmp_c = sig_t1;
    data_tmp_f = [real(data_tmp_c); imag(data_tmp_c)];
    data_tmp_f = data_tmp_f(:);
    fwrite(fileID,data_tmp_f,'float');
    fclose(fileID);

end
disp('Beamforming vectors are applied to cellular')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%% Device to Device %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if wifi_enabled==1;
%     run('wifi_precode.m')
% else
     load d_calib
%     num_sample_shift=10^4;
%     %% PTRS pattern
%     ptrs_cliche=zeros(n_sc_per_rb,n_sym_per_sf);
%     ptrs_val=1;
%     ptrs_cliche(end,:)=ptrs_val;
%     ptrs_pad=repmat(ptrs_cliche,nrb,n_sf_per_f);              % pattern of dmrs during whole resource grid
%     %% USER 1 with DMRS pattern 1
%     cliche1=squeeze(dmrs_patterns(15,:,:));% dmrs pattern in one time slot
%     dmrs_pad1=repmat(cliche1,nrb,n_sf_per_f);              % pattern of dmrs during whole resource grid
%     %ref_sig_val=1;
%     data_pad1=ones(size(dmrs_pad1))-dmrs_pad1-ptrs_pad;
%     cliche1(find(cliche1))=[1 -1 1 -1 1 -1 -1 1 -1 -1 1];
%     dmrs_pad1=repmat(cliche1,nrb,n_sf_per_f);              % pattern of dmrs during whole resource grid
%     % freq domain data
%     data_freq=1/sqrt(2)*(2*(randi(2,size(dmrs_pad1))-1)-1)+1j*1/sqrt(2)*(2*(randi(2,size(dmrs_pad1))-1)-1);
%     data_freq=data_freq.*data_pad1;
%     data_freq(find(dmrs_pad1))=dmrs_pad1(find(dmrs_pad1));
%     data_freq(find(ptrs_pad))=ptrs_val;
%     %data_freq=[p(2,:).' data_freq];
%     %% Precoding
%     num_sample_process = nfft*150;
%     for ii = 1:n_ant_d
%         fileID = fopen(['signals/intf_' num2str(ii-1) '.dat'], 'rb');
%         if (fileID < 0)
%             error('Error: fail to open files!');
%         end
%         frewind(fileID);
%         fseek(fileID, num_sample_shift, 'bof');
%         data_ant_f = fread(fileID, 2*num_sample_process, 'float');
%         data_ant_c = transpose(data_ant_f(1:2:end) + 1i*data_ant_f(2:2:end));
%         rx(ii,:)= data_ant_c(1:num_sample_process);
%         fclose(fileID);
%     end
%     rx=rx*1000;
%     frx1=reshape(rx(1,:),nfft,[]);
%     frx2=reshape(rx(2,:),nfft,[]);
%     frx1=fft(frx1);
%     frx2=fft(frx2);
%     pv=zeros(2,nfft,28);
%     grp=32;
%     rr=512/grp;
%     for jj=1:rr
%         y=[reshape(frx1((jj-1)*grp+1:jj*grp,:),1,[]);reshape(frx1((jj-1)*grp+1:jj*grp,:),1,[])];
%         YYH=y*y';
%         [U S V]=svd(YYH);
%         for ii=(jj-1)*grp+1:jj*grp
%             pv(:,ii,:)=repmat(conj(V(:,2)),1,28);
%         end
%     end
%     %%
%     pv=1/sqrt(2)*ones(size(pv));
%     for ii=1:n_ant_d;
%         v=squeeze(pv(ii,:,:));
%         frame_f=zeros(512,28);
%         frame_f(valid_sc_indx,:)=data_freq;
%         frame_f=frame_f.*v;
%         frame_t=ifft(frame_f)*sqrt(nfft);
%         frame_t_c = [frame_t(end-(cp_length-1):end,:); frame_t];
%         frame_t_c=[pn(2,:).',frame_t_c];
%         %frame_t_c=[preamble(:,1), frame_t_c];
%         sig_t=frame_t_c(:).';
%         sig_t1=repmat(sig_t,1,n_frame);
%         sig_t1=[zeros(1,140) sig_t1];
%         % intentional CFO
%         delta_phase_per_sample = 0.00;
%         rx_dds = exp(-1i*delta_phase_per_sample*(1:size(sig_t1,2)));
%         rx_dds_frame = repmat(rx_dds,  size(sig_t1,1), 1);
%         sig_t1=sig_t1.*rx_dds_frame;
%         sig_t1=sig_t1*d_calib(ii);
%         % write file
%         fileID= fopen(['signals/d2d_bf_tx_' num2str(ii-1) '.dat'], 'wb');
%         if (fileID < 0)
%             error('Error: fail to open files!');
%         end
%         data_tmp_c = sig_t1;
%         data_tmp_f = [real(data_tmp_c); imag(data_tmp_c)];
%         data_tmp_f = data_tmp_f(:);
%         fwrite(fileID,data_tmp_f,'float');
%         fclose(fileID);
%     end
%     disp('Beamforming vectors are applied to D2D communications')
% end