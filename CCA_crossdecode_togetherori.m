clear all;
weightfile = 'p9_V1LM_CCA_loc_together_allz.mat';
spikefile = 'p9_V1LM_CCA_stat_together_allz.mat';
load(weightfile);
Aweight = A;
Bweight = B;
cweight = c;
clearvars -except Aweight Bweight cweight weightfile spikefile
load(spikefile);
cspike = c;
clear A B r p;
r = [];
p = [];
nsub = 5;
winwidth = 150;
winstep = 10;

un_oris = [0 30 60 90];
for h = 1:10
    for i = 1:size(V_matonez,3)
        for j = (-100/winstep):(100/winstep)
            if i+j > 0 && i+j < size(V_matonez,3)
                testsize = sum(test(c,h));
                trainsize = sum(training(c,h));
                Utrain = transpose(squeeze(Aweight(h,:,1:nsub,i,j+(100/winstep)+1))) * squeeze(V_matonez(:,training(c,h),i));
                Vtrain = transpose(squeeze(Bweight(h,:,1:nsub,i,j+(100/winstep)+1))) * squeeze(LM_matonez(:,training(c,h),i+j));
                UVtrain = transpose(vertcat(Utrain,Vtrain));
                Utest = transpose(squeeze(Aweight(h,:,1:nsub,i,j+(100/winstep)+1))) * squeeze(V_matonez(:,test(c,h),i));
                Vtest = transpose(squeeze(Bweight(h,:,1:nsub,i,j+(100/winstep)+1))) * squeeze(LM_matonez(:,test(c,h),i+j));
                UVtest = transpose(vertcat(Utest,Vtest));
                cc = corrcoef(Utest(1,:),Vtest(1,:));
                cc(1,2);
                r(h,i,j+(100/winstep)+1) = cc(1,2);
%                 if r(h,i,j+21) < 0
%                     r(h,i,j+21) = NaN;
%                 end
                Umdllda = fitcdiscr(transpose(Utrain),ori3(training(c,h)),'DiscrimType','pseudolinear');
                ptemp = predict(Umdllda, transpose(Utest));
                pu(h,i,j+(100/winstep)+1) = sum(transpose(ptemp) == ori3(test(c,h)))/testsize;

                Vmdllda = fitcdiscr(transpose(Vtrain),ori3(training(c,h)),'DiscrimType','pseudolinear');
                ptemp = predict(Vmdllda, transpose(Vtest));
                pv(h,i,j+(100/winstep)+1) = sum(transpose(ptemp) == ori3(test(c,h)))/testsize;  

                UVmdllda = fitcdiscr(UVtrain,ori3(training(c,h)),'DiscrimType','pseudolinear');
                ptemp = predict(UVmdllda, UVtest);
                p(h,i,j+(100/winstep)+1) = sum(transpose(ptemp) == ori3(test(c,h)))/testsize; 
                p(h,i,j+(100/winstep)+1)
            else 
                r(h,i,j+(100/winstep)+1) = NaN;
                p(h,i,j+(100/winstep)+1) = NaN;
                pu(h,i,j+(100/winstep)+1) = NaN;
                pv(h,i,j+(100/winstep)+1) = NaN;
            end
        end
    end
end

ravg = squeeze(mean(r,1,'omitnan'));
ravg(ravg < 0) = NaN;
pavg = squeeze(mean(p,1,'omitnan'));

%%
figure;
h = heatmap(squeeze(ravg(:,:)),'Colormap',parula,'GridVisible','off');
h.YDisplayData = flipud(h.YDisplayData);
h.XDisplayLabels = nan(size(h.XDisplayData));
h.YDisplayLabels = nan(size(h.YDisplayData));
h.XDisplayLabels(1:(size(h.XDisplayLabels)-1)/4:size(h.XDisplayLabels)) = num2cell(-100:50:100);
h.YDisplayLabels(1:10:size(h.YDisplayLabels)) = num2cell(800:-(50):0);
h.XLabel = 'Delay';
h.YLabel = 'Time from Stimulus Onset';
h.Title = ['LGN-V1 CCA Correlation Coefficiencts (Locomoting Subspace, Stationary Activity, orientations calculated together'];
% savefig('../figs/together/LGNV1_CCAcorrelations_loconstat_together_allz')

figure;
h = heatmap(squeeze(pavg(:,:)),'Colormap',parula,'GridVisible','off');
h.YDisplayData = flipud(h.YDisplayData);
h.XDisplayLabels = nan(size(h.XDisplayData));
h.YDisplayLabels = nan(size(h.YDisplayData));
h.XDisplayLabels(1:(size(h.XDisplayLabels)-1)/4:size(h.XDisplayLabels)) = num2cell(-100:50:100);
h.YDisplayLabels(1:10:size(h.YDisplayLabels)) = num2cell(800:-(50):0);
h.XLabel = 'Delay';
h.YLabel = 'Time from Stimulus Onset';
h.Title = 'LGN-V1 CCA Decoding Heatmap (Locomoting Subspace, Stationary Activity, orientations calculated together';
% savefig('../figs/together/LGNV1_CCAperformance_loconstat_together_allz')


save('p9_V1LM_CCA_statonloc_together_allz.mat','r','p','pu','pv','c')


%%
figure;
t = 0:5:845;
plot(t,plgnstat(:,21),'b')
hold on;
plot(t,plgnloc(:,21),'r')
hold on;
plot(t,plgnlos(:,21),'g')
hold on;
plot(t,plgnsol(:,21),'m')
hold on;
xlabel('Time from stimulus onset (ms)')
ylabel('Decoding performance')
title('LGN')
legend('stat','loc','loc on stat','stat on loc');



figure;
t = 0:5:845;
plot(t,pvstat(:,21),'b')
hold on;
plot(t,pvloc(:,21),'r')
hold on;
plot(t,pvlos(:,21),'g')
hold on;
plot(t,pvsol(:,21),'m')
hold on;
xlabel('Time from stimulus onset (ms)')
ylabel('Decoding performance')
title('V1 (from LGN-V1)')
legend('stat','loc','loc on stat','stat on loc');

mydata = [mean(plgnstat(:,21),'omitnan') mean(plgnsol(:,21),'omitnan') mean(plgnloc(:,21),'omitnan') mean(plgnlos(:,21),'omitnan')];
color= {'b','m','r','g'};
figure, hold on
% % if data is more than colors then colors will be repeated
m = length(color);
for k = 1:length(mydata)
    i = mod(k-1,m); %%i is remainder after division of k-1 by m
    i = i+1;    
    h=bar(k,mydata(k));
    set(h,'FaceColor',color{i});
end
h=gca; h.XAxis.TickLength = [0 0];
xticks(1:4)
xticklabels({'Stationary activity + subspace','Stationary activity + Locomoting subspace','Locomoting activity + subspace','Locomoting activity + Stationary subspace'});
ylim([0.25 0.40])
% legend('Locomoting activity + subspace','Locomoting activity + Stationary subspace','Stationary activity + subspace','Stationary activity + Locomoting subspace');
title('LGN decoding')

mydata = [mean(pvstat(:,21),'omitnan') mean(pvsol(:,21),'omitnan') mean(pvloc(:,21),'omitnan') mean(pvlos(:,21),'omitnan')];
color= {'b','m','r','g'};
figure, hold on
% % if data is more than colors then colors will be repeated
m = length(color);
for k = 1:length(mydata)
    i = mod(k-1,m); %%i is remainder after division of k-1 by m
    i = i+1;    
    h=bar(k,mydata(k));
    set(h,'FaceColor',color{i});
end
h=gca; h.XAxis.TickLength = [0 0];
xticks(1:4)
xticklabels({'Stationary activity + subspace','Stationary activity + Locomoting subspace','Locomoting activity + subspace','Locomoting activity + Stationary subspace'});
ylim([0.25 0.40])
% legend('Locomoting activity + subspace','Locomoting activity + Stationary subspace','Stationary activity + subspace','Stationary activity + Locomoting subspace');
title('V1 (from LGN-V1) decoding')


%%
figure;
t = 0:5:845;
plot(t,plmstat(:,21),'b')
hold on;
plot(t,plmloc(:,21),'r')
hold on;
plot(t,plmlos(:,21),'g')
hold on;
plot(t,plmsol(:,21),'m')
hold on;
xlabel('Time from stimulus onset (ms)')
ylabel('Decoding performance')
title('LM')
legend('stat','loc','loc on stat','stat on loc');



figure;
t = 0:5:845;
plot(t,pvstat(:,21),'b')
hold on;
plot(t,pvloc(:,21),'r')
hold on;
plot(t,pvlos(:,21),'g')
hold on;
plot(t,pvsol(:,21),'m')
hold on;
xlabel('Time from stimulus onset (ms)')
ylabel('Decoding performance')
title('V1 (V1-LM)')
legend('stat','loc','loc on stat','stat on loc');

mydata = [mean(plmstat(:,21),'omitnan') mean(plmsol(:,21),'omitnan') mean(plmloc(:,21),'omitnan') mean(plmlos(:,21),'omitnan')];
color= {'b','m','r','g'};
figure, hold on
% % if data is more than colors then colors will be repeated
m = length(color);
for k = 1:length(mydata)
    i = mod(k-1,m); %%i is remainder after division of k-1 by m
    i = i+1;    
    h=bar(k,mydata(k));
    set(h,'FaceColor',color{i});
end
h=gca; h.XAxis.TickLength = [0 0];
xticks(1:4)
xticklabels({'Stationary activity + subspace','Stationary activity + Locomoting subspace','Locomoting activity + subspace','Locomoting activity + Stationary subspace'});
ylim([0.25 0.33])
% legend('Locomoting activity + subspace','Locomoting activity + Stationary subspace','Stationary activity + subspace','Stationary activity + Locomoting subspace');
title('LM decoding')

mydata = [mean(pvstat(:,21),'omitnan') mean(pvsol(:,21),'omitnan') mean(pvloc(:,21),'omitnan') mean(pvlos(:,21),'omitnan')];
color= {'b','m','r','g'};
figure, hold on
% % if data is more than colors then colors will be repeated
m = length(color);
for k = 1:length(mydata)
    i = mod(k-1,m); %%i is remainder after division of k-1 by m
    i = i+1;    
    h=bar(k,mydata(k));
    set(h,'FaceColor',color{i});
end
h=gca; h.XAxis.TickLength = [0 0];
xticks(1:4)
xticklabels({'Stationary activity + subspace','Stationary activity + Locomoting subspace','Locomoting activity + subspace','Locomoting activity + Stationary subspace'});
ylim([0.25 0.40])
% legend('Locomoting activity + subspace','Locomoting activity + Stationary subspace','Stationary activity + subspace','Stationary activity + Locomoting subspace');
title('V1 (from V1-LM) decoding')

%%
figure;
t = 0:5:845;
plot(t,plgnstat(:,21),'b')
hold on;
plot(t,plgnloc(:,21),'r')
hold on;
plot(t,plgnshuffstat(:,21),'black')
hold on;
xlabel('Time from stimulus onset (ms)')
ylabel('Decoding performance')
title('LGN')
legend('stationary','locomoting','shuffled stationary');
%%
mydata = [mean(pv2stat(:,21),'omitnan') mean(pv2shuffstat(:,21),'omitnan') mean(pv2loc(:,21),'omitnan') mean(pv2shuff(:,21),'omitnan')];
color= {'b','black','r','g'};
figure, hold on
% % if data is more than colors then colors will be repeated
m = length(color);
for k = 1:length(mydata)
    i = mod(k-1,m); %%i is remainder after division of k-1 by m
    i = i+1;    
    h=bar(k,mydata(k));
    set(h,'FaceColor',color{i});
end
h=gca; h.XAxis.TickLength = [0 0];
xticks(1:4)
xticklabels({'Stationary','Shuffled Stationary','Locomoting','Shuffled Locomoting'});
ylim([0.25 0.4])
% legend('Locomoting activity + subspace','Locomoting activity + Stationary subspace','Stationary activity + subspace','Stationary activity + Locomoting subspace');
title('V1 (from V1-LM) decoding')