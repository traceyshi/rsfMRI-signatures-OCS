function [feats] = selectFeatsPearsonNum(Y,X,num)
% Arguments:
%   Y: features to be selected among
%   X: features that Y features are to be correlated with for selection
%   num: number of features to return (sorted by lowest p-value)
%
% Returns:
%   feats: indices of selected features

[~, pvals] = corr(Y, X, 'type', 'Pearson'); 

minp = min(pvals, [], 2);

[B,I] = sort(minp);

feats = sort(I(1:num));
