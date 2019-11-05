%% Simplified NR_TX for two cellular users in the UL
clc
clear all
close all
%% Parameters
num_user=3;
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
%% Preamble
load pn
preamble=[pn(1,:).',pn(2,:).',pn(3,:).'];
%% DMRS patterns
load dmrs_patterns
%% PTRS pattern
ptrs_cliche=zeros(n_sc_per_rb,n_sym_per_sf);
ptrs_val=1;
ptrs_cliche(end,:)=ptrs_val;
ptrs_pad=repmat(ptrs_cliche,nrb,n_sf_per_f);              % pattern of dmrs during whole resource grid
for ix=1:num_user
cliche=squeeze(dmrs_patterns(ix,:,:));                    % dmrs pattern in one time slot    
dmrs_pad=repmat(cliche,nrb,n_sf_per_f);                   % pattern of dmrs during whole resource grid  
data_pad=ones(size(dmrs_pad))-dmrs_pad-ptrs_pad;
cliche(find(cliche))=[1 -1 1 -1 1 -1 -1 1 -1 -1 1]/2;
dmrs_pad=repmat(cliche,nrb,n_sf_per_f);                   % pattern of dmrs during whole resource grid
% freq domain data
data_freq=1/sqrt(2)*(2*(randi(2,size(dmrs_pad))-1)-1)+1j*1/sqrt(2)*(2*(randi(2,size(dmrs_pad))-1)-1);
data_freq=data_freq.*data_pad;
data_freq(find(dmrs_pad))=dmrs_pad(find(dmrs_pad));
data_freq(find(ptrs_pad))=ptrs_val;
frame_f=zeros(512,28);
frame_f(valid_sc_indx,:)=data_freq;
frame_t=ifft(frame_f)*sqrt(nfft);
frame_t_c = [frame_t(end-(cp_length-1):end,:); frame_t];
frame_t_c_2=[preamble(:,ix), frame_t_c];
sig_t=frame_t_c_2(:).';
figure
plot(abs(sig_t))
sig_t=repmat(sig_t,1,100);
sig_t_1=[zeros(1,(ix-1)*40) sig_t];
% intentional CFO
% delta_phase_per_sample = 0.00;
% rx_dds = exp(-1i*delta_phase_per_sample*(1:size(sig_t1,2)));
% rx_dds_frame = repmat(rx_dds,  size(sig_t1,1), 1);
% sig_t1=sig_t1.*rx_dds_frame;
% write file
fileID= fopen(['signals/nr_tx_',num2str(ix-1),'.dat'], 'wb');
if (fileID < 0)
    error('Error: fail to open files!');
end
data_tmp_c = sig_t_1;
data_tmp_f = [real(data_tmp_c); imag(data_tmp_c)];
data_tmp_f = data_tmp_f(:);
fwrite(fileID,data_tmp_f,'float');
fclose(fileID);
end