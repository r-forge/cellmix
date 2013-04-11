# Using the SCOREM approach one can filter marker genes that have consisten
# whose expression values in a given datasets.
# This is because the notion of marker can be very dependent on the data, biological
# conditions or platform used to define them.

# load kidney transplant dataset from Shen-Orr et al. (2010)
eset <- ExpressionMix('GSE20300', verbose=TRUE)

# load HaemAtlas marker gene list (on Illumina)
ml <- MarkerList('HaemAtlas')
ml
# convert marker IDs to Affy IDs from dataset
# using stringent one to one mapping
mla <- convertIDs(ml, eset, method='1:1', verbose=TRUE)
summary(mla)

# plot expression profile of each set of markers (only the first 10)
profplot(mla[, 1:10], eset, split=TRUE, lab='')

# filter out using SCOREM
sml <- extractMarkers(mla, eset, method='scorem', alpha=.005)
# plot selected markers
profplot(sml, eset, split=TRUE, lab='')
