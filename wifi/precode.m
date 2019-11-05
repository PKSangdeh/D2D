%% calculate precoders
clc
close all
clear all
%%
ap_antenna=3;
num_sample_process = 64*1000;
num_sample_shift = 2e6;
num_move_sample=1;
rx=zeros(ap_antenna,num_sample_process);

apr=2;



POS_VALID_SC = [39:64 2:27];
% load dl_chan.mat
% load ul_chan.mat
% average over 20 samples of channel estimation
% dl_avg_ch=zeros(2,length(POS_VALID_SC));
% for i=1:20;
%     dl_avg_ch=dl_avg_ch+squeeze(dl_chan(:,:,i));
% end
% dl_avg_ch=dl_avg_ch/length(POS_VALID_SC);
% ul_avg_ch=zeros(2,length(POS_VALID_SC));
%% Capture 1000 samples of received signal
for ii = 1:ap_antenna
    fileID = fopen(['signals/intf_signal_rx' num2str(ii-1) '.dat'], 'rb');
    if (fileID < 0)
        error('Error: fail to open files!');
    end
    frewind(fileID);
    fseek(fileID, num_sample_shift, 'bof');
    data_ant_f = fread(fileID, 2*num_sample_process, 'float');
    data_ant_c = transpose(data_ant_f(1:2:end) + 1i*data_ant_f(2:2:end));
    rx(ii,:)= data_ant_c(1:num_sample_process);
    fclose(fileID);
end
rx=rx*1000;
frx1=reshape(rx(1,:),64,[]);
frx2=reshape(rx(2,:),64,[]);
frx3=reshape(rx(3,:),64,[]);
r=[rx(1,:);rx(2,:);rx(3,:)];
frx1=fft(frx1);
frx2=fft(frx2);
frx3=fft(frx3);
sing_val1=zeros(ap_antenna,64);
sing_val2=sing_val1;

figure; plot(abs(rx(:,1:100).'))

%% frequency-domain data
%% 
% calib=chan_calib_coeff;
% calib_amp=abs(calib);
% calib_ang=angle(calib);
P1=zeros(ap_antenna,64);
% approach 1
P2=P1;
%% Approach 2
for jj=1:64
    YYH=zeros(ap_antenna,ap_antenna);
    y=[frx1(jj,:);frx2(jj,:);frx3(jj,:)];
    for ii=1:length(y(1,ii))
        YYH=YYH+y(:,ii)*y(:,ii)';
    end
    YYH=YYH/ii;
    [U S V]=svd(YYH);
    P2(:,jj)=V(:,3);
    sing_val2(:,jj)=[S(1,1);S(2,2);S(3,3)];
end
% %% approach 3- ZF as a benchmark - uses uplink channel coefficient+calibration
% ul_avg_ch=zeros(ap_antenna,length(POS_VALID_SC));
% for i=1:1;
%     ul_avg_ch=ul_avg_ch+squeeze(ul_chan(:,:,i));
% end
% ul_avg_ch=ul_avg_ch/i;
% P3=zeros(ap_antenna,64);
% for i=1:52
%     pr1=null(ul_avg_ch(:,i).');
%     pr=pr1(:,1);
%     pr=pr/norm(pr);
%     P3(:,POS_VALID_SC(i))=pr;
% end
% %% appr 4
% dl_avg_ch=zeros(ap_antenna,length(POS_VALID_SC));
% for i=1:1;
%     dl_avg_ch=dl_avg_ch+squeeze(dl_chan(:,:,i));
% end
% dl_avg_ch=dl_avg_ch/i;
% P4=zeros(ap_antenna,64);
% for i=1:52
%     pr1=null(dl_avg_ch(:,i).');
%     pr=pr1(:,1);
%     pr=pr/norm(pr);
%     P4(:,POS_VALID_SC(i))=pr;
% end
% %% Approach 3
% nn=length(frx1(1,:));
% for i1=1:64
%     tmp_p=[0;0];
%     tmp_2=[];
% for i2=1:nn
%     tmp1=[frx1(i1,:);frx2(i1,:)];
%     [U S V]=svd(tmp1(:,i2)*tmp1(:,i2)');
%     tmp_p=tmp_p+V(:,2);
%     tmp_2=[tmp_2 V(:,2)];
% end
% save V.mat tmp_2
% tmp_p=tmp_p/nn;
% P(:,i1)=tmp_p/norm(tmp_p);
% end
    





% %% Generating Precoded Data;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Parameters
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% num_ant_ap = ap_antenna;
% num_valid_sc = 52;
% const_num_sc = 64;
% num_symbol_per_frame = 20;
% const_num_frame = 50;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % data parameters
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% STF 
% STF_POS = [3  7  11  15  19  23  30  34  38  42  46  50];
% STF_DAT_GRP = zeros(1, 52);    
% STF_DAT_GRP(1,STF_POS) = sqrt(2)*[-1-1i  -1+1i  -1-1i  -1-1i  -1-1i  1-1i  -1-1i  -1-1i  1+1i  1-1i  1+1i  -1+1i];
% %%% LTF
% LTF_DAT_GRP = zeros(1, 52);
% LTF_DAT_GRP(1,:) = [-1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1 1 1 -1 -1 -1];
% %%%% pilot
% POS_VALID_SC = [39:64 2:27];
% POS_PAYLOAD_SC = [1:5 7:19 21:32 34:46 48:52];
% POS_PILOT_SC = [6 20 33 47];
% PILOT_DAT_GRP = zeros(4,num_symbol_per_frame);
% PILOT_DAT_GRP(1,:) = [-1 -1 -1 -1 1 1 1 -1 1 -1 -1 1 1 1 1 -1 -1 1 -1 1];
% PILOT_DAT_GRP(2,:) = [1 -1 -1 1 -1 1 1 1 1 1 1 -1 1 1 1 -1 -1 1 1 1];
% PILOT_DAT_GRP(3,:) = [-1 -1 1 1 -1 -1 -1 1 1 1 -1 -1 1 -1 1 1 1 1 1 1];
% PILOT_DAT_GRP(4,:) = [1 1 1 -1 1 1 -1 -1 -1 -1 1 -1 -1 -1 -1 -1 -1 1 1 -1]; 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % generate freq data
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% freq_data_grp = zeros(num_valid_sc, num_symbol_per_frame);
%     stf_dat = STF_DAT_GRP(1,:);
%     ltf_dat = LTF_DAT_GRP(1,:);
%     data_freq = zeros(num_valid_sc,num_symbol_per_frame);
%     data_freq(:,1) = stf_dat.';
%     data_freq(:,2) = stf_dat.';
%     data_freq(:,3) = ltf_dat.';
%     data_freq(:,4) = ltf_dat.';
% for jj = 5:num_symbol_per_frame
%     data_freq(:,jj) = (1/sqrt(2))*(sign(randn(num_valid_sc,1))+1i*sign(randn(num_valid_sc,1)));
%     data_freq(POS_PILOT_SC,jj) = PILOT_DAT_GRP(:,jj);
% end
% data_f_frm = zeros(const_num_sc, num_symbol_per_frame);
% data_f_frm(POS_VALID_SC,:) = data_freq;
% f_tx1=zeros(const_num_sc, num_symbol_per_frame);
% f_tx2=zeros(const_num_sc, num_symbol_per_frame);
% f_tx3=zeros(const_num_sc, num_symbol_per_frame);
% if apr==1
%    for ii=1:64
%        f_tx1(ii,:)=data_f_frm(ii,:)*conj(P1(1,ii));
%        f_tx2(ii,:)=data_f_frm(ii,:)*conj(P1(2,ii));
%        f_tx3(ii,:)=data_f_frm(ii,:)*conj(P1(3,ii));
%    end
% elseif apr==2   
%    for ii=1:64
%        f_tx1(ii,:)=data_f_frm(ii,:)*(P2(1,ii));
%        f_tx2(ii,:)=data_f_frm(ii,:)*conj(P2(2,ii));
%        f_tx3(ii,:)=data_f_frm(ii,:)*conj(P2(3,ii));
%    end
% elseif apr==3
%    for ii=1:64
%        f_tx1(ii,:)=data_f_frm(ii,:)*P3(1,ii);
%        f_tx2(ii,:)=data_f_frm(ii,:)*P3(2,ii);
%        f_tx3(ii,:)=data_f_frm(ii,:)*P3(3,ii);
%    end   
%    elseif apr==4
%    for ii=1:64
%        f_tx1(ii,:)=data_f_frm(ii,:)*P4(1,ii);
%        f_tx2(ii,:)=data_f_frm(ii,:)*P4(2,ii);
%        f_tx3(ii,:)=data_f_frm(ii,:)*P4(3,ii);
%    end 
% end
% time_signal_tx1=ifft(f_tx1)*sqrt(const_num_sc);
% time_signal_tx2=ifft(f_tx2)*sqrt(const_num_sc);
% time_signal_tx3=ifft(f_tx3)*sqrt(const_num_sc);
% tx1=reshape(time_signal_tx1,1,[]);
% tx2=reshape(time_signal_tx2,1,[]);
% tx3=reshape(time_signal_tx3,1,[]);
% tx1=repmat(tx1,1,1000);
% tx2=repmat(tx2,1,1000);
% tx3=repmat(tx3,1,1000);
% %% channel resiprosity consideration
% load  chan_calib_coeff.mat
% if apr==4
% tx1=1*tx1;
% tx2=1*tx2;
% tx3=1*tx3;
% else
% tx1=chan_calib_coeff(1)*tx1;
% tx2=chan_calib_coeff(2)*tx2;
% tx3=chan_calib_coeff(3)*tx3;
% end
% time_signal_grp=[tx1;tx2;tx3];
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % write to the file
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for ii = 1:ap_antenna
%     fileID= fopen(['precoded_signal' num2str(ii-1) '.dat'], 'wb');
%     if (fileID < 0)
%         disp('Error: fail to open files!');
%         pause;
%     end
%     data_tmp_c = time_signal_grp(ii,:);
%     data_tmp_f = [real(data_tmp_c); imag(data_tmp_c)];
%     data_tmp_f = data_tmp_f(:);
%     fwrite(fileID,data_tmp_f,'float');
%     fclose(fileID);
% end

