clear all;
weightfile = 'p11_LGNV1_CCA_loc_together_allz.mat';
spikefile = 'p11_LGNV1_CCA_loc_together_allz.mat';

load(weightfile);
Aweight = A;
Bweight = B;
cweight = c;
clearvars -except Aweight Bweight cweight weightfile spikefile
load(spikefile);
cspike = c;
clear A B r p;
winstep = 10;
winwidth = 150;
r = [];
pu = [];
pv = [];
plgnfull = [];
pvfull = [];
plgnfulltrain = [];
pvfulltrain = [];
putrain = [];
pvtrain = [];
un_oris = [0 30 60 90];
for h = 1:10
    for i = 1:size(V_matonez,3)
        for j = (-100/winstep):(100/winstep)
            if i+j > 0 && i+j < size(V_matonez,3)
                testsize = sum(test(c,h));
                trainsize = sum(training(c,h));
                Utrain = transpose(squeeze(Aweight(h,:,1:5,i,j+(100/winstep)+1))) * squeeze(LGN_matonez(:,training(c,h),i));
                Vtrain = transpose(squeeze(Bweight(h,:,1:5,i,j+(100/winstep)+1))) * squeeze(V_matonez(:,training(c,h),i+j));
                UVtrain = transpose(vertcat(Utrain,Vtrain));
                Utest = transpose(squeeze(Aweight(h,:,1:5,i,j+(100/winstep)+1))) * squeeze(LGN_matonez(:,test(c,h),i));
                Vtest = transpose(squeeze(Bweight(h,:,1:5,i,j+(100/winstep)+1))) * squeeze(V_matonez(:,test(c,h),i+j));
                UVtest = transpose(vertcat(Utest,Vtest));
                cc = corrcoef(Utest(1,:),Vtest(1,:));
                cc(1,2);
                r(h,i,j+(100/winstep)+1) = cc(1,2);
%                 if r(h,i,j+(100/winstep)+1) < 0
%                     r(h,i,j+(100/winstep)+1) = NaN;
%                 end
                Umdllda = fitcdiscr(transpose(Utrain),ori3(training(c,h)),'DiscrimType','pseudolinear');
                ptemp = predict(Umdllda, transpose(Utest));
                pu(h,i,j+(100/winstep)+1) = sum(transpose(ptemp) == ori3(test(c,h)))/testsize;

                ptemp = predict(Umdllda, transpose(Utrain));
                putrain(h,i,j+(100/winstep)+1) = sum(transpose(ptemp) == ori3(training(c,h)))/trainsize; 

                Vmdllda = fitcdiscr(transpose(Vtrain),ori3(training(c,h)),'DiscrimType','pseudolinear');
                ptemp = predict(Vmdllda, transpose(Vtest));
                pv(h,i,j+(100/winstep)+1) = sum(transpose(ptemp) == ori3(test(c,h)))/testsize;  

                ptemp = predict(Vmdllda, transpose(Vtrain));
                pvtrain(h,i,j+(100/winstep)+1) = sum(transpose(ptemp) == ori3(training(c,h)))/trainsize; 

                LGNmdlldafull = fitcdiscr(transpose(LGN_matonez(:,training(c,h),i)),ori3(training(c,h)),'DiscrimType','pseudolinear');
                ptemp = predict(LGNmdlldafull, transpose(LGN_matonez(:,test(c,h),i)));
                plgnfull(h,i,j+(100/winstep)+1) = sum(transpose(ptemp) == ori3(test(c,h)))/testsize;
                plgnfull(h,i,j+(100/winstep)+1)


                ptemp = predict(LGNmdlldafull, transpose(LGN_matonez(:,training(c,h),i)));
                plgnfulltrain(h,i,j+(100/winstep)+1) = sum(transpose(ptemp) == ori3(training(c,h)))/trainsize;
                plgnfulltrain(h,i,j+(100/winstep)+1)

                Vmdlldafull = fitcdiscr(transpose(V_matonez(:,training(c,h),i)),ori3(training(c,h)),'DiscrimType','pseudolinear');
                ptemp = predict(Vmdlldafull, transpose(V_matonez(:,test(c,h),i)));
                pvfull(h,i,j+(100/winstep)+1) = sum(transpose(ptemp) == ori3(test(c,h)))/testsize;

                ptemp = predict(Vmdlldafull, transpose(V_matonez(:,training(c,h),i)));
                pvfulltrain(h,i,j+(100/winstep)+1) = sum(transpose(ptemp) == ori3(training(c,h)))/trainsize;

            else
                r(h,i,j+(100/winstep)+1) = NaN;
                pu(h,i,j+(100/winstep)+1) = NaN;
                pv(h,i,j+(100/winstep)+1) = NaN;
                pvtrain(h,i,j+(100/winstep)+1) = NaN;
                putrain(h,i,j+(100/winstep)+1) = NaN;
                pvfull(h,i,j+(100/winstep)+1) = NaN;
                plgnfull(h,i,j+(100/winstep)+1) = NaN;
                plgnfulltrain(h,i,j+(100/winstep)+1) = NaN;
                pvfulltrain(h,i,j+(100/winstep)+1) = NaN;

            end
        end
    end
end

ravg = squeeze(mean(r,1,'omitnan'));
ravg(ravg < 0) = NaN;

% plgnavg = squeeze(mean(pu(:,:,:),1,'omitnan'));
% plgnavg(plgnavg < 0) = NaN;

%%
figure;
h = heatmap(squeeze(ravg(:,:)),'Colormap',parula,'GridVisible','off');
h.YDisplayData = flipud(h.YDisplayData);
h.XDisplayLabels = nan(size(h.XDisplayData));
h.YDisplayLabels = nan(size(h.YDisplayData));
h.XDisplayLabels(1:(size(h.XDisplayLabels)-1)/4:size(h.XDisplayLabels)) = num2cell(-100:50:100);
h.YDisplayLabels(5:(50/winstep):size(h.YDisplayLabels)) = num2cell((1000-winwidth-50):-(50):0);
h.XLabel = 'Delay';
h.YLabel = 'Time from Stimulus Onset';
h.Title = [{'LGN-V1 CCA'},{'Correlations (Stationary)'}];
% savefig('../figs/together/LGNV1_CCAcorrelations_loconstat_together_allz')

% 
%%
figure;
h = heatmap(squeeze(mean(pvfulltrain,1)),'Colormap',parula,'GridVisible','off');
h.YDisplayData = flipud(h.YDisplayData);
h.XDisplayLabels = nan(size(h.XDisplayData));
h.YDisplayLabels = nan(size(h.YDisplayData));
h.XDisplayLabels(1:(size(h.XDisplayLabels)-1)/4:size(h.XDisplayLabels)) = num2cell(-100:50:100);
h.YDisplayLabels(5:(50/winstep):size(h.YDisplayLabels)) = num2cell((1000-winwidth-50):-(50):0);
h.XLabel = 'Delay';
h.YLabel = 'Time from Stimulus Onset';
h.Title = 'LGN-V1 CCA Decoding Heatmap (Locomoting Subspace, Stationary Activity, orientations calculated together';
% savefig('../figs/together/LGNV1_CCAperformance_loconstat_together_allz')
%%

save('p11_LGNV1_CCAtesttrain_loc_together_individual.mat','r','pu','pv','putrain','pvtrain','plgnfull','pvfull','plgnfulltrain','pvfulltrain','c')
%%
clear all;
weightfile = 'p11_V1LM_CCA_loc_together_allz.mat';
spikefile = 'p11_V1LM_CCA_loc_together_allz.mat';

load(weightfile);
Aweight = A;
Bweight = B;
cweight = c;
clearvars -except Aweight Bweight cweight weightfile spikefile
load(spikefile);
cspike = c;
clear A B r p;
r = [];
pu = [];
pv = [];
pvfull = [];
plmfull = [];
pvfulltrain = [];
plmfulltrain = [];
putrain = [];
pvtrain = [];
un_oris = [0 30 60 90];
for h = 1:10
    for i = 1:size(V_matonez,3)
        for j = (-100/winstep):(100/winstep)
            if i+j > 0 && i+j < size(V_matonez,3)
                testsize = sum(test(c,h));
                trainsize = sum(training(c,h));

                Utrain = transpose(squeeze(Aweight(h,:,1:5,i,j+(100/winstep)+1))) * squeeze(V_matonez(:,training(c,h),i));
                Vtrain = transpose(squeeze(Bweight(h,:,1:5,i,j+(100/winstep)+1))) * squeeze(LM_matonez(:,training(c,h),i+j));
                UVtrain = transpose(vertcat(Utrain,Vtrain));
                Utest = transpose(squeeze(Aweight(h,:,1:5,i,j+(100/winstep)+1))) * squeeze(V_matonez(:,test(c,h),i));
                Vtest = transpose(squeeze(Bweight(h,:,1:5,i,j+(100/winstep)+1))) * squeeze(LM_matonez(:,test(c,h),i+j));
                UVtest = transpose(vertcat(Utest,Vtest));
                cc = corrcoef(Utest(1,:),Vtest(1,:));
                cc(1,2);
                r(h,i,j+(100/winstep)+1) = cc(1,2);
%                 if r(h,i,j+(100/winstep)+1) < 0
%                     r(h,i,j+(100/winstep)+1) = NaN;
%                 end
                Umdllda = fitcdiscr(transpose(Utrain),ori3(training(c,h)),'DiscrimType','pseudolinear');
                ptemp = predict(Umdllda, transpose(Utest));
                pu(h,i,j+(100/winstep)+1) = sum(transpose(ptemp) == ori3(test(c,h)))/testsize;

                ptemp = predict(Umdllda, transpose(Utrain));
                putrain(h,i,j+(100/winstep)+1) = sum(transpose(ptemp) == ori3(training(c,h)))/trainsize; 

                Vmdllda = fitcdiscr(transpose(Vtrain),ori3(training(c,h)),'DiscrimType','pseudolinear');
                ptemp = predict(Vmdllda, transpose(Vtest));
                pv(h,i,j+(100/winstep)+1) = sum(transpose(ptemp) == ori3(test(c,h)))/testsize;  

                ptemp = predict(Vmdllda, transpose(Vtrain));
                pvtrain(h,i,j+(100/winstep)+1) = sum(transpose(ptemp) == ori3(training(c,h)))/trainsize; 

                Vmdlldafull = fitcdiscr(transpose(V_matonez(:,training(c,h),i)),ori3(training(c,h)),'DiscrimType','pseudolinear');
                ptemp = predict(Vmdlldafull, transpose(V_matonez(:,test(c,h),i)));
                pvfull(h,i,j+(100/winstep)+1) = sum(transpose(ptemp) == ori3(test(c,h)))/testsize;
                pvfull(h,i,j+(100/winstep)+1)


                ptemp = predict(Vmdlldafull, transpose(V_matonez(:,training(c,h),i)));
                pvfulltrain(h,i,j+(100/winstep)+1) = sum(transpose(ptemp) == ori3(training(c,h)))/trainsize;
                pvfulltrain(h,i,j+(100/winstep)+1)

                LMmdlldafull = fitcdiscr(transpose(LM_matonez(:,training(c,h),i)),ori3(training(c,h)),'DiscrimType','pseudolinear');
                ptemp = predict(LMmdlldafull, transpose(LM_matonez(:,test(c,h),i)));
                plmfull(h,i,j+(100/winstep)+1) = sum(transpose(ptemp) == ori3(test(c,h)))/testsize;

                ptemp = predict(LMmdlldafull, transpose(LM_matonez(:,training(c,h),i)));
                plmfulltrain(h,i,j+(100/winstep)+1) = sum(transpose(ptemp) == ori3(training(c,h)))/trainsize;

            else
                r(h,i,j+(100/winstep)+1) = NaN;
                pu(h,i,j+(100/winstep)+1) = NaN;
                pv(h,i,j+(100/winstep)+1) = NaN;
                pvtrain(h,i,j+(100/winstep)+1) = NaN;
                putrain(h,i,j+(100/winstep)+1) = NaN;
                plmfull(h,i,j+(100/winstep)+1) = NaN;
                pvfull(h,i,j+(100/winstep)+1) = NaN;
                plvfulltrain(h,i,j+(100/winstep)+1) = NaN;
                plmfulltrain(h,i,j+(100/winstep)+1) = NaN;

            end
        end
    end
end

ravg = squeeze(mean(r,1,'omitnan'));
ravg(ravg < 0) = NaN;
%%
figure;
h = heatmap(squeeze(ravg(:,:)),'Colormap',parula,'GridVisible','off');
h.YDisplayData = flipud(h.YDisplayData);
h.XDisplayLabels = nan(size(h.XDisplayData));
h.YDisplayLabels = nan(size(h.YDisplayData));
h.XDisplayLabels(1:(size(h.XDisplayLabels)-1)/4:size(h.XDisplayLabels)) = num2cell(-100:50:100);
h.YDisplayLabels(5:(50/winstep):size(h.YDisplayLabels)) = num2cell((1000-winwidth-50):-(50):0);
h.XLabel = 'Delay';
h.YLabel = 'Time from Stimulus Onset';
h.Title = ['V1-LM CCA Correlation Coefficiencts (Locomoting Subspace, Stationary Activity, orientations calculated together'];
% savefig('../figs/together/VV1_CCAcorrelations_loconstat_together_allz')
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


save('p11_V1LM_CCAtesttrain_loc_together_individual.mat','r','pu','pv','putrain','pvtrain','plmfull','pvfull','plmfulltrain','pvfulltrain','c')



%%
xaxis = [0:5:849];
figure;
plot(xaxis, mean(plgnfulltrain(:,:,11),1),'b')
hold on;
plot(xaxis, mean(plgnfull(:,:,11),1),'black')
hold on;
plot(xaxis, mean(putrain(:,:,11),1),'g')
hold on;
plot(xaxis, mean(pu(:,:,11),1),'r')
hold on;
xlim([0 845])
legend('LGN fullpop training','LGN fullpop test','LGN CCA training','LGN CCA test')
xlabel('Time to Stimulus Onset (ms)')
ylabel('Decoding Performance')
title('LGN decoding CCA vs full population training vs test')


%%
xaxis = [0:5:849];
figure;
plot(xaxis, mean(pvfulltrain(:,:,11),1),'b')
hold on;
plot(xaxis, mean(pvfull(:,:,11),1),'black')
hold on;
plot(xaxis, mean(pvtrain(:,:,11),1),'g')
hold on;
plot(xaxis, mean(pv(:,:,11),1),'r')
hold on;
xlim([0 845])
legend('V1 fullpop training','V1 fullpop test','V1 CCA training','V1 CCA test')
xlabel('Time to Stimulus Onset (ms)')
ylabel('Decoding Performance')
title('V1 decoding CCA vs full population training vs test')

%%
xaxis = [0:5:849];
figure;
plot(xaxis, mean(plmfulltrain(:,:,11),1),'b')
hold on;
plot(xaxis, mean(plmfull(:,:,11),1),'black')
hold on;
plot(xaxis, mean(pvtrain(:,:,11),1),'g')
hold on;
plot(xaxis, mean(pv(:,:,11),1),'r')
hold on;
xlim([0 845])
legend('LM fullpop training','LM fullpop test','LM CCA training','LM CCA test')
xlabel('Time to Stimulus Onset (ms)')
ylabel('Decoding Performance')
title('LM decoding CCA vs full population training vs test')