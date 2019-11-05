%% De-Modulation Reference Signals (DMRS) patterns
%%
clc
clear all
close all
%%
plot_dmrs=1;
ref_sig_val=1;
%% Pattern 1- CUE
p1=zeros(12,14);
indx=[37:47];
p1(indx)=ref_sig_val;
%% Pattern 2 - CUE
p2=zeros(12,14);
indx=[37:47];
indx=[indx+84];
p2(indx)=ref_sig_val;

%% Pattern 3 - CUE
p3=zeros(12,14);
indx=[1:11];
p3(indx)=ref_sig_val;

%% Pattern 4- D2D

p4=zeros(12,14);
indx=[85:95];
p4(indx)=ref_sig_val;

%%
dmrs_patterns=zeros(4,12,14);
dmrs_patterns(1,:,:)=p1;
dmrs_patterns(2,:,:)=p2;
dmrs_patterns(3,:,:)=p3;
dmrs_patterns(4,:,:)=p4;
%%
tmp=1;
while tmp
    ch=menu('DMRS patterns in one slot','Pattern 1-CUE','Pattern 2-CUE','Pattern 3-CUE','Pattern 4-D2D','Exit');
    if ch==5
        tmp=0;
        clc
        close all
    else
        p=squeeze(dmrs_patterns(ch,:,:));
        figure
        for i=1:12
            for j=1:14
                if p(i,j)~=0
                    plot(j,i,'rs','linewidth',22)
                    hold on
                else
                    plot(j,i,'gs','linewidth',22)
                    hold on
                end
            end
        end
        axis([0.5 14.5 0.5 12.5])
    end
end
save('dmrs_patterns.mat','dmrs_patterns')