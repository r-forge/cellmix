####################################################
# Sample deconvolution analysis with DSA in CellMix
####################################################

# load benchmark data
x <- ExpressionMix('GSE19830', verbose=TRUE)
dim(x)
annotation(x)
# extract mixed samples
mix <- mixedSamples(x)

# load TIGER marker list
ml <- MarkerList('TIGER')
ml
names(ml)

# select markers for the tissues present in the mixture
basisnames(x)
ml <- ml[c('brain', 'liver', 'lung')]
summary(ml)

# convert to match annotations
mlx <- convertIDs(ml, mix, verbose=TRUE)
summary(mlx)

# QC on markers from their expression patterns in mixed samples
profplot(mlx[,1:10], mix)
# filter out poor markers using SCOREM (based on linear-scale expression)
mlsc <- extractMarkers(mlx, expb(mix, 2), method='SCOREM', alpha=10^-12)
summary(mlsc)
# expresison patterns are more correlated
profplot(mlsc[,1:10], mix)

# apply DSA using all markers
res <- ged(mix[mlsc,], mlsc, 'DSA', verbose=TRUE)

# plot against true proportions
profplot(mix, res)
