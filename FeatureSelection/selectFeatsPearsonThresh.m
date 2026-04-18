function [feats] = selectFeatsPearsonThresh(Y,X,thresh)
% Arguments:
%   Y: features to be selected among
%   X: features that Y features are to be correlated with for selection
%   thresh: p-value threshold (features in Y are selected if they have
%       a Spearman rho with p-value < thresh with at least one X)
%
% Returns:
%   feats: indices of selected features

[~, pvals] = corr(Y, X, 'type', 'Pearson'); 

minp = min(pvals, [], 2);

feats = find(minp < thresh);