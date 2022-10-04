clear all;
weightfile = 'p9_LGNV1_CCA_stat_together_allz.mat';
spikefile = 'p9_LGNV1_CCA_stat_together_allz.mat';

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
pu = [];
pv = [];
winstep = 10;
winwidth = 150;
un_oris = [0 30 60 90];
nsub = 5;
for g = 1:nsub
    for h = 1:10
        for i = 1:size(V_matonez,3)
            for j = (-100/winstep):(100/winstep)
                if i+j > 0 && i+j < size(V_matonez,3)
                    testsize = sum(test(c,h));
                    trainsize = sum(training(c,h));
                    if g == 1
                        Utrain = transpose(squeeze(Aweight(h,:,1:g,i,j+(100/winstep)+1)) * squeeze(LGN_matonez(:,training(c,h),i)));
                        Vtrain = transpose(squeeze(Bweight(h,:,1:g,i,j+(100/winstep)+1)) * squeeze(V_matonez(:,training(c,h),i+j)));
                        UVtrain = horzcat(Utrain,Vtrain);
                        Utest = transpose(squeeze(Aweight(h,:,1:g,i,j+(100/winstep)+1)) * squeeze(LGN_matonez(:,test(c,h),i)));
                        Vtest = transpose(squeeze(Bweight(h,:,1:g,i,j+(100/winstep)+1)) * squeeze(V_matonez(:,test(c,h),i+j)));
                        UVtest = horzcat(Utest,Vtest);
                    else
                        Utrain = transpose(transpose(squeeze(Aweight(h,:,1:g,i,j+(100/winstep)+1))) * squeeze(LGN_matonez(:,training(c,h),i)));
                        Vtrain = transpose(transpose(squeeze(Bweight(h,:,1:g,i,j+(100/winstep)+1))) * squeeze(V_matonez(:,training(c,h),i+j)));
                        UVtrain = horzcat(Utrain,Vtrain);
                        Utest = transpose(transpose(squeeze(Aweight(h,:,1:g,i,j+(100/winstep)+1))) * squeeze(LGN_matonez(:,test(c,h),i)));
                        Vtest = transpose(transpose(squeeze(Bweight(h,:,1:g,i,j+(100/winstep)+1))) * squeeze(V_matonez(:,test(c,h),i+j)));
                        UVtest = horzcat(Utest,Vtest);
                    end
                    cc = corrcoef(Utest(:,g),Vtest(:,g));
                    cc(1,2);
                    r(g,h,i,j+(100/winstep)+1) = cc(1,2);
                    Umdllda = fitcdiscr(Utrain,ori3(training(c,h)),'DiscrimType','pseudolinear');
                    ptemp = predict(Umdllda, Utest);
                    pu(g,h,i,j+(100/winstep)+1) = sum(transpose(ptemp) == ori3(test(c,h)))/testsize;

                    Vmdllda = fitcdiscr(Vtrain,ori3(training(c,h)),'DiscrimType','pseudolinear');
                    ptemp = predict(Vmdllda, Vtest);
                    pv(g,h,i,j+(100/winstep)+1) = sum(transpose(ptemp) == ori3(test(c,h)))/testsize;

                    UVmdllda = fitcdiscr(UVtrain,ori3(training(c,h)),'DiscrimType','pseudolinear');
                    ptemp = predict(UVmdllda, UVtest);
                    p(g,h,i,j+(100/winstep)+1) = sum(transpose(ptemp) == ori3(test(c,h)))/testsize;
                    p(g,h,i,j+(100/winstep)+1)
                else
                    r(g,h,i,j+(100/winstep)+1) = NaN;
                    pu(g,h,i,j+(100/winstep)+1) = NaN;
                    pv(g,h,i,j+(100/winstep)+1) = NaN;
                    p(g,h,i,j+(100/winstep)+1) = NaN;
                end
            end
        end
    end
end

ravg = squeeze(mean(r,1,'omitnan'));
ravg(ravg < 0) = NaN;
%%
% figure;
% h = heatmap(squeeze(ravg(:,:)),'Colormap',parula,'GridVisible','off');
% h.YDisplayData = flipud(h.YDisplayData);
% h.XDisplayLabels = nan(size(h.XDisplayData));
% h.YDisplayLabels = nan(size(h.YDisplayData));
% h.XDisplayLabels(1:(size(h.XDisplayLabels)-1)/4:size(h.XDisplayLabels)) = num2cell(-100:50:100);
% h.YDisplayLabels(1:10:size(h.YDisplayLabels)) = num2cell(800:-(50):0);
% h.XLabel = 'Delay';
% h.YLabel = 'Time from Stimulus Onset';
% h.Title = ['LGN-V1 CCA Correlation Coefficiencts (Locomoting Subspace, Stationary Activity, orientations calculated together'];
% % savefig('../figs/together/LGNV1_CCAcorrelations_loconstat_together_allz')
% % 
% %%
% figure;
% h = heatmap(squeeze(mean(pvfulltrain,1)),'Colormap',parula,'GridVisible','off');
% h.YDisplayData = flipud(h.YDisplayData);
% h.XDisplayLabels = nan(size(h.XDisplayData));
% h.YDisplayLabels = nan(size(h.YDisplayData));
% h.XDisplayLabels(1:(size(h.XDisplayLabels)-1)/4:size(h.XDisplayLabels)) = num2cell(-100:50:100);
% h.YDisplayLabels(1:10:size(h.YDisplayLabels)) = num2cell(800:-(50):0);
% h.XLabel = 'Delay';
% h.YLabel = 'Time from Stimulus Onset';
% h.Title = 'LGN-V1 CCA Decoding Heatmap (Locomoting Subspace, Stationary Activity, orientations calculated together';
% savefig('../figs/together/LGNV1_CCAperformance_loconstat_together_allz')
%%

save('p11_LGNV1_CCAsubspacenum_stat_together_individual.mat','r','pu','pv','p','c')
%%
clear all;
weightfile = 'p11_V1LM_CCA_stat_together_allz.mat';
spikefile = 'p11_V1LM_CCA_stat_together_allz.mat';

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
pu = [];
pv = [];
un_oris = [0 30 60 90];
for g = 1:nsub
    for h = 1:10
        for i = 1:size(V_matonez,3)
            for j = (-100/winstep):(100/winstep)
                if i+j > 0 && i+j < size(V_matonez,3)
                    testsize = sum(test(c,h));
                    trainsize = sum(training(c,h));
                    if g == 1
                        Utrain = transpose(squeeze(Aweight(h,:,1:g,i,j+(100/winstep)+1)) * squeeze(V_matonez(:,training(c,h),i)));
                        Vtrain = transpose(squeeze(Bweight(h,:,1:g,i,j+(100/winstep)+1)) * squeeze(LM_matonez(:,training(c,h),i+j)));
                        UVtrain = horzcat(Utrain,Vtrain);
                        Utest = transpose(squeeze(Aweight(h,:,1:g,i,j+(100/winstep)+1)) * squeeze(V_matonez(:,test(c,h),i)));
                        Vtest = transpose(squeeze(Bweight(h,:,1:g,i,j+(100/winstep)+1)) * squeeze(LM_matonez(:,test(c,h),i+j)));
                        UVtest = horzcat(Utest,Vtest);
                    else
                        Utrain = transpose(transpose(squeeze(Aweight(h,:,1:g,i,j+(100/winstep)+1))) * squeeze(V_matonez(:,training(c,h),i)));
                        Vtrain = transpose(transpose(squeeze(Bweight(h,:,1:g,i,j+(100/winstep)+1))) * squeeze(LM_matonez(:,training(c,h),i+j)));
                        UVtrain = horzcat(Utrain,Vtrain);
                        Utest = transpose(transpose(squeeze(Aweight(h,:,1:g,i,j+(100/winstep)+1))) * squeeze(V_matonez(:,test(c,h),i)));
                        Vtest = transpose(transpose(squeeze(Bweight(h,:,1:g,i,j+(100/winstep)+1))) * squeeze(LM_matonez(:,test(c,h),i+j)));
                        UVtest = horzcat(Utest,Vtest);
                    end
                    cc = corrcoef(Utest(:,g),Vtest(:,g));
                    cc(1,2);
                    r(g,h,i,j+(100/winstep)+1) = cc(1,2);
                    Umdllda = fitcdiscr(Utrain,ori3(training(c,h)),'DiscrimType','pseudolinear');
                    ptemp = predict(Umdllda, Utest);
                    pu(g,h,i,j+(100/winstep)+1) = sum(transpose(ptemp) == ori3(test(c,h)))/testsize;

                    Vmdllda = fitcdiscr(Vtrain,ori3(training(c,h)),'DiscrimType','pseudolinear');
                    ptemp = predict(Vmdllda, Vtest);
                    pv(g,h,i,j+(100/winstep)+1) = sum(transpose(ptemp) == ori3(test(c,h)))/testsize;

                    UVmdllda = fitcdiscr(UVtrain,ori3(training(c,h)),'DiscrimType','pseudolinear');
                    ptemp = predict(UVmdllda, UVtest);
                    p(g,h,i,j+(100/winstep)+1) = sum(transpose(ptemp) == ori3(test(c,h)))/testsize;
                    p(g,h,i,j+(100/winstep)+1)

                else
                    r(g,h,i,j+(100/winstep)+1) = NaN;
                    pu(g,h,i,j+(100/winstep)+1) = NaN;
                    pv(g,h,i,j+(100/winstep)+1) = NaN;
                    p(g,h,i,j+(100/winstep)+1) = NaN;
                end
            end
        end
    end
end

ravg = squeeze(mean(r,1,'omitnan'));
ravg(ravg < 0) = NaN;
%%
% figure;
% h = heatmap(squeeze(ravg(:,:)),'Colormap',parula,'GridVisible','off');
% h.YDisplayData = flipud(h.YDisplayData);
% h.XDisplayLabels = nan(size(h.XDisplayData));
% h.YDisplayLabels = nan(size(h.YDisplayData));
% h.XDisplayLabels(1:(size(h.XDisplayLabels)-1)/4:size(h.XDisplayLabels)) = num2cell(-100:50:100);
% h.YDisplayLabels(1:10:size(h.YDisplayLabels)) = num2cell(800:-(50):0);
% h.XLabel = 'Delay';
% h.YLabel = 'Time from Stimulus Onset';
% h.Title = ['V1-LM CCA Correlation Coefficiencts (Locomoting Subspace, Stationary Activity, orientations calculated together'];
% % savefig('../figs/together/VV1_CCAcorrelations_loconstat_together_allz')
% 
% figure;
% h = heatmap(squeeze(pavg(:,:)),'Colormap',parula,'GridVisible','off');
% h.YDisplayData = flipud(h.YDisplayData);
% h.XDisplayLabels = nan(size(h.XDisplayData));
% h.YDisplayLabels = nan(size(h.YDisplayData));
% h.XDisplayLabels(1:(size(h.XDisplayLabels)-1)/4:size(h.XDisplayLabels)) = num2cell(-100:50:100);
% h.YDisplayLabels(1:10:size(h.YDisplayLabels)) = num2cell(800:-(50):0);
% h.XLabel = 'Delay';
% h.YLabel = 'Time from Stimulus Onset';
% h.Title = 'V-V1 CCA Decoding Heatmap (Locomoting Subspace, Stationary Activity, orientations calculated together';
% savefig('../figs/together/VV1_CCAperformance_loconstat_together_allz')


save('p11_V1LM_CCAsubspacenum_stat_together_individual.mat','r','pu','pv','p','c')



%%
xaxis = [0:5:849];
figure;
plot(xaxis, mean(plgnfulltrain(:,:,21),1),'b')
hold on;
plot(xaxis, mean(plgnfull(:,:,21),1),'black')
hold on;
plot(xaxis, mean(putrain(:,:,21),1),'g')
hold on;
plot(xaxis, mean(pu(:,:,21),1),'r')
hold on;
xlim([0 845])
legend('LGN fullpop training','LGN fullpop test','LGN CCA training','LGN CCA test')
xlabel('Time to Stimulus Onset (ms)')
ylabel('Decoding Performance')
title('LGN decoding CCA vs full population training vs test')


%%
xaxis = [0:5:849];
figure;
plot(xaxis, mean(pvfulltrain(:,:,21),1),'b')
hold on;
plot(xaxis, mean(pvfull(:,:,21),1),'black')
hold on;
plot(xaxis, mean(pvtrain(:,:,21),1),'g')
hold on;
plot(xaxis, mean(pv(:,:,21),1),'r')
hold on;
xlim([0 845])
legend('V1 fullpop training','V1 fullpop test','V1 CCA training','V1 CCA test')
xlabel('Time to Stimulus Onset (ms)')
ylabel('Decoding Performance')
title('V1 decoding CCA vs full population training vs test')

%%
xaxis = [0:5:849];
figure;
plot(xaxis, mean(plmfulltrain(:,:,21),1),'b')
hold on;
plot(xaxis, mean(plmfull(:,:,21),1),'black')
hold on;
plot(xaxis, mean(pvtrain(:,:,21),1),'g')
hold on;
plot(xaxis, mean(pv(:,:,21),1),'r')
hold on;
xlim([0 845])
legend('LM fullpop training','LM fullpop test','LM CCA training','LM CCA test')
xlabel('Time to Stimulus Onset (ms)')
ylabel('Decoding Performance')
title('LM decoding CCA vs full population training vs test')

%%
ravg = squeeze(mean(r,2));
pavg = squeeze(mean(p,2));
puavg = squeeze(mean(pu,2));
pvavg = squeeze(mean(pv,2));
rmax = [];
for i = 1:nsub
    figure;
    rmax(i) = mean(max(squeeze(ravg(i,:,:)),[],2));
    h = heatmap(squeeze(ravg(i,:,:)),'Colormap',parula,'GridVisible','off');
    h.YDisplayData = flipud(h.YDisplayData);
    h.XDisplayLabels = nan(size(h.XDisplayData));
    h.YDisplayLabels = nan(size(h.YDisplayData));
    h.XDisplayLabels(1:(size(h.XDisplayLabels)-1)/4:size(h.XDisplayLabels)) = num2cell(-100:50:100);
    h.YDisplayLabels(5:(50/winstep):size(h.YDisplayLabels)) = num2cell((1000-winwidth-50):-(50):0);
    h.XLabel = 'Delay';
    h.YLabel = 'Time from Stimulus Onset';
    h.Title = strcat('LGN-V1 subspace ', num2str(i), ' correlation ',' (Locomoting)');
end
%%
close all
m = [];
mu = [];
mv = [];
for i = 1:5
    figure;
    m(i) = mean(max(squeeze(pavg(i,:,:)),[],2));
    mu(i) = mean(max(squeeze(puavg(i,:,:)),[],2));
    mv(i) = mean(max(squeeze(pvavg(i,:,:)),[],2));
    h = heatmap(squeeze(pavg(i,:,:)),'Colormap',parula,'GridVisible','off');
    h.YDisplayData = flipud(h.YDisplayData);
    h.XDisplayLabels = nan(size(h.XDisplayData));
    h.YDisplayLabels = nan(size(h.YDisplayData));
    h.XDisplayLabels(1:(size(h.XDisplayLabels)-1)/4:size(h.XDisplayLabels)) = num2cell(-100:50:100);
    h.YDisplayLabels(5:(50/winstep):size(h.YDisplayLabels)) = num2cell((1000-winwidth-50):-(50):0);
    h.XLabel = 'Delay';
    h.YLabel = 'Time from Stimulus Onset';
    h.Title = strcat('LGN-V1 ',{' '}, num2str(i), ' subspaces decoding ',' (Locomoting)');
end

figure;
plot(m,'black');
hold on;
plot(mu,'b');
hold on;
plot(mv,'r');
hold on;
plot(rmax,'g')
legend('LGN + V1','LGN','V1','LGN-V1 Correlations')
xlabel('Number of subspaces')
ylabel('Decoding performance')
title('LGN/V1 Decoding by number of subspaces (stationary)')

%%
ravg = squeeze(mean(r,2));
pavg = squeeze(mean(p,2));
puavg = squeeze(mean(pu,2));
pvavg = squeeze(mean(pv,2));
rmax = [];
for i = 1:5
    figure;
    rmax(i) = mean(max(squeeze(ravg(i,:,:)),[],2));
    h = heatmap(squeeze(ravg(i,:,:)),'Colormap',parula,'GridVisible','off');
    h.YDisplayData = flipud(h.YDisplayData);
    h.XDisplayLabels = nan(size(h.XDisplayData));
    h.YDisplayLabels = nan(size(h.YDisplayData));
    h.XDisplayLabels(1:(size(h.XDisplayLabels)-1)/4:size(h.XDisplayLabels)) = num2cell(-100:50:100);
    h.YDisplayLabels(1:10:size(h.YDisplayLabels)) = num2cell(800:-(50):0);
    h.XLabel = 'Delay';
    h.YLabel = 'Time from Stimulus Onset';
    h.Title = strcat('V1-LM subspace ', num2str(i), ' correlation ',' (Locomoting)');
end
%%
close all
m = [];
mu = [];
mv = [];
for i = 1:5
    figure;
    m(i) = max(squeeze(pavg(i,:,:)),[],'all');
    mu(i) = max(squeeze(puavg(i,:,:)),[],'all');
    mv(i) = max(squeeze(pvavg(i,:,:)),[],'all');
    h = heatmap(squeeze(pavg(i,:,:)),'Colormap',parula,'GridVisible','off');
    h.YDisplayData = flipud(h.YDisplayData);
    h.XDisplayLabels = nan(size(h.XDisplayData));
    h.YDisplayLabels = nan(size(h.YDisplayData));
    h.XDisplayLabels(1:(size(h.XDisplayLabels)-1)/4:size(h.XDisplayLabels)) = num2cell(-100:50:100);
    h.YDisplayLabels(5:(50/winstep):size(h.YDisplayLabels)) = num2cell((1000-winwidth-50):-(50):0);
    h.XLabel = 'Delay';
    h.YLabel = 'Time from Stimulus Onset';
    h.Title = strcat('V1-LM ',{' '}, num2str(i), ' subspaces decoding ',' (Locomoting)');
end

figure;
plot(m,'black');
hold on;
plot(mu,'b');
hold on;
plot(mv,'r');
hold on;
plot(rmax,'g')
legend('V1 + LM','V1','LM','V1-LM Correlations')
xlabel('Number of subspaces')
ylabel('Decoding performance')
title('V1/LM Decoding by number of subspaces (locomoting)')

%%
figure;
plot(mlmstat,'b');
hold on;
plot(rvlmstat,'c');
hold on;
plot(mlmloc,'r');
hold on;
plot(rvlmloc,'Color',[0.8500 0.3250 0.0980])
legend('V1 stationary decoding','stationary correlation','V1 locomotiong decoding','locomoting correlation')
xlabel('Number of subspaces')
ylabel('Decoding performance/correlation')
title('V1 (from LGN-V1) Decoding/correlation by number of subspaces')