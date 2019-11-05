%% generate PN sequence
clc
clear all
close all
%%
cp_length=38;
Nz=300;
n_valid_sc=300;                 % Number of valid subcarrier
nfft=512;                       % FFT and IFFT points
valid_sc=[363:512,2:151];  % Index of valid subcarriers
n_frame=1;                     % Number of NR frames
%% pn 1
p1=2*(randi(2,1,300)-1)-1;
pn1_f=zeros(1,512);
pn1_f(valid_sc)=p1;
pn1_f=ifft(pn1_f,nfft)*sqrt(nfft);
pn1=[pn1_f(end-(cp_length-1):end),pn1_f];
%% pn 3
p2=2*(randi(2,1,300)-1)-1;
pn2_f=zeros(1,512);
pn2_f(valid_sc)=p2;
pn2_f=ifft(pn2_f,nfft)*sqrt(nfft);
pn2=[pn2_f(end-(cp_length-1):end),pn2_f];
%% pn 3
p3=2*(randi(2,1,300)-1)-1;
pn3_f=zeros(1,512);
pn3_f(valid_sc)=p3;
pn3_f=ifft(pn3_f,nfft)*sqrt(nfft);
pn3=[pn3_f(end-(cp_length-1):end),pn3_f];
%% pn 4
p4=2*(randi(2,1,300)-1)-1;
pn4_f=zeros(1,512);
pn4_f(valid_sc)=p4;
pn4_f=ifft(pn4_f,nfft)*sqrt(nfft);
pn4=[pn4_f(end-(cp_length-1):end),pn4_f];
%% 
y=xcorr(pn1,pn2)/(norm(pn1)*norm(pn2));
figure
plot(abs(y))
title('cross-12')
%% 
y=xcorr(pn1,pn3)/(norm(pn1)*norm(pn3));
figure
plot(abs(y))
title('cross-13')
%% 
y=xcorr(pn2,pn3)/(norm(pn3)*norm(pn2));
figure
plot(abs(y))
title('cross-23')

%%
%% 
y=xcorr(pn1,pn4)/(norm(pn1)*norm(pn4));
figure
plot(abs(y))
title('cross-14')
%% 
y=xcorr(pn2,pn4)/(norm(pn4)*norm(pn2));
figure
plot(abs(y))
title('cross-24')
%% 
y=xcorr(pn4,pn3)/(norm(pn4)*norm(pn3));
figure
plot(abs(y))
title('cross-34')

%%
pn=[pn1;pn2;pn3;pn4];
save('pn.mat','pn')