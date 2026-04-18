% Select features for use in CCA

%% SETUP
%NYSPI desktop
datadir = '/path/to/data/';
scriptdir = '/path/to/scripts/FeatureSelection';
addpath(scriptdir)

cd(datadir)
X = readmatrix('Data_prep/clin_features_scaled.csv'); %clinical features
fid = fopen('Data_prep/rsFC_features_scaled.csv'); %rs-FC features
fmt = [repmat('%f',1,61776)];
Yreg = textscan(fid, fmt, 'collectoutput',true,'headerlines',1,'delimiter',',');
Yreg = Yreg{1,1};
fclose(fid);

%% feature selection by univariate correlations

%full discovery subsample
Pearsonthresh = 0.05;
allsubsample_selectedFeatsPearson = selectFeatsPearsonThresh(Yreg, X, Pearsonthresh);
num_feats = size(allsubsample_selectedFeatsPearson,1)

%within each cross-validation subset
cd('PartitionsFeatureSelection/')
partitions = dir('partitions_rep*')
CV_feats = cell(length(partitions),1); %initialize empty cell array
for p = 1:length(partitions)
    disp(partitions(p).name)
    thisSubset = table2array(readtable(partitions(p).name));
    %CV_feats{p} = selectFeatsSpearmanThresh(fullresidY(thisSubset,:), fullresidX(thisSubset,:), Spearmanthresh);
    %CV_feats{p} = selectFeatsPearsonNum(fullresidY(thisSubset,:), fullresidX(thisSubset,:), num_feats);
    CV_feats{p} = selectFeatsPearsonThresh(Yreg(thisSubset,:), X(thisSubset,:), Pearsonthresh);
end
cd('..')

% Export selected features
%writecell(CV_feats,'crossValidFeatsPearsonThresh05.csv')
%writematrix(alldiscovery_selectedFeatsSpearman,'alldiscovery_selectedFeatsPearson05.csv')
%writecell(CV_feats,'crossValidFeatsSpearman.csv')
%writematrix(alldiscovery_selectedFeatsSpearman,'alldiscovery_selectedFeatsSpearman.csv')

% Tabulate most stable features across subsets
CV_feats_concat = [];
for p = 1:length(partitions)
    CV_feats_concat = [CV_feats_concat; CV_feats{p,1}];
end

%CV_feats_concat = [CV_feats{1,1}; CV_feats{2,1}; CV_feats{3,1}; CV_feats{4,1}; CV_feats{5,1};
%    CV_feats{6,1}; CV_feats{7,1}; CV_feats{8,1}; CV_feats{9,1}; CV_feats{10,1}];
CV_feats_stability = tabulate(CV_feats_concat);
histogram(CV_feats_stability(:,2))
tabulate(CV_feats_stability(:,2))
stable_feats_95 = find(CV_feats_stability(:,2) >= 95);
%check that all of these stable features are also selected in the full
%discovery subsample
tabulate(ismember(stable_feats_95, allsubsample_selectedFeatsPearson))

%% Export stable features
cd(datadir)
cd('../results')

writematrix(stable_feats_95,'replication_stable_feats_95_Pearson05.csv')

%% Assess stability of discovery-selected features in held-out test subset
cd(datadir)
cd('../results')

X_test = readmatrix('clin_features_scaled_test.csv'); %clinical variables after zero-inflated poisson covariate regression (implemented in R, Main.Rmd)
fid = fopen('rsFC_features_scaled_test.csv'); %rs-FC variables after covariate regression (implemented in R, Main.Rmd)
fmt = [repmat('%f',1,61776)];
Yreg_test = textscan(fid, fmt, 'collectoutput',true,'headerlines',1,'delimiter',',');
Yreg_test = Yreg_test{1,1};

[test_corrs, test_corr_p] = corr(X_test, Yreg_test(:, stable_feats_95));
test_corr_minp = min(test_corr_p);
sum(test_corr_minp <= .05)
hist(test_corr_minp)

[test_corrs_all, test_corr_p_all] = corr(X_test, Yreg_test);
test_corr_minp_all = min(test_corr_p_all);
sum(test_corr_minp_all <= .05)
hist(test_corr_minp_all)


%% Assess stability of discovery-selected features in held-out replication subset
cd(datadir)
cd('replication')

X_replication = readmatrix('clin_features_scaled.csv'); %clinical features
fid = fopen('rsFC_features_scaled.csv'); %rs-FC features
fmt = [repmat('%f',1,61776)];
Yreg_replication = textscan(fid, fmt, 'collectoutput',true,'headerlines',1,'delimiter',',');
Yreg_replication = Yreg_replication{1,1};

[repl_corrs, repl_corr_p] = corr(X_replication, Yreg_replication(:, stable_feats_95));
repl_corr_minp = min(repl_corr_p);
sum(repl_corr_minp <= .05)





