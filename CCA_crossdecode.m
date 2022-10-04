clear all;
weightfile = 'p11_V1LM_stat_CCA_allz_sepori.mat';
spikefile = 'p11_V1LM_loc_CCA_allz_sepori.mat';

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
un_oris = [0 30 60 90];
for h = 1:10
    for i = 1:size(V_matonez,3)
        for j = -20:20
            if i+j > 0 && i+j < size(V_matonez,3)
                oritrain = [];
                oritest = [];
                for g = 1:length(un_oris)
                    f = find(ori3 == un_oris(g));
                    Utrain = transpose(squeeze(Aweight(g,h,:,1:5,i,j+21))) * squeeze(V_matonez(:,f(training(c,h)),i));
                    Vtrain = transpose(squeeze(Bweight(g,h,:,1:5,i,j+21))) * squeeze(LM_matonez(:,f(training(c,h)),i+j));
                    UVtrain{g} = transpose(vertcat(Utrain,Vtrain));
                    Utest = transpose(squeeze(Aweight(g,h,:,1:5,i,j+21))) * squeeze(V_matonez(:,f(test(c,h)),i));
                    Vtest = transpose(squeeze(Bweight(g,h,:,1:5,i,j+21))) * squeeze(LM_matonez(:,f(test(c,h)),i+j));
                    UVtest{g} = transpose(vertcat(Utest,Vtest));
                    cc = corrcoef(Utest(1,:),Vtest(1,:));
                    cc(1,2);
                    r(g,h,i,j+21) = cc(1,2);
                    if r(g,h,i,j+21) < 0
                        r(g,h,i,j+21) = NaN;
                    end
                    oritrain = [oritrain ori3(f(training(c,h)))];
                    oritest = [oritest ori3(f(test(c,h)))];
                end
                UVmdllda = fitcdiscr(cat(1,UVtrain{:}),oritrain,'DiscrimType','pseudolinear');
                ptemp = predict(UVmdllda, cat(1,UVtest{:}));
                p(h,i,j+21) = sum(transpose(ptemp) == oritest);
                p(h,i,j+21)
            else
                r(1:length(un_oris),h,i,j+21) = NaN;
                p(h,i,j+21) = NaN;
            end
        end
    end
end

ravg = squeeze(mean(r,2,'omitnan'));
pavg = squeeze(sum(p,1)/length(ori3));
%%
figure;
h = heatmap(squeeze(mean(ravg(:,:,:),1)),'Colormap',parula,'GridVisible','off');
h.YDisplayData = flipud(h.YDisplayData);
h.XDisplayLabels = nan(size(h.XDisplayData));
h.YDisplayLabels = nan(size(h.YDisplayData));
h.XDisplayLabels(1:(size(h.XDisplayLabels)-1)/4:size(h.XDisplayLabels)) = num2cell(-100:50:100);
h.YDisplayLabels(1:10:size(h.YDisplayLabels)) = num2cell(800:-(50):0);
h.XLabel = 'Delay';
h.YLabel = 'Time from Stimulus Onset';
h.Title = ['V1-LM CCA Correlation Coefficiencts (Stationary Subspace, Locomoting Activity, orientations calculated seperately'];

figure;
h = heatmap(squeeze(pavg(:,:)),'Colormap',parula,'GridVisible','off');
h.YDisplayData = flipud(h.YDisplayData);
h.XDisplayLabels = nan(size(h.XDisplayData));
h.YDisplayLabels = nan(size(h.YDisplayData));
h.XDisplayLabels(1:(size(h.XDisplayLabels)-1)/4:size(h.XDisplayLabels)) = num2cell(-100:50:100);
h.YDisplayLabels(1:10:size(h.YDisplayLabels)) = num2cell(800:-(50):0);
h.XLabel = 'Delay';
h.YLabel = 'Time from Stimulus Onset';
h.Title = 'V1-LM CCA Decoding Heatmap (Stationary Subspace, Locomoting Activity, orientations calculated seperately';