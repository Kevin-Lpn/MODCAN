getwd();
workingDir = "";
setwd(workingDir);

library(WGCNA);
library(rhdf5);
options(stringsAsFactors = FALSE);

# ========================================================================
# load reference samples
# row: gene, columns: samples
exp = read.table("./exp_data.txt",
                 sep = "\t");
normalized_exp = as.data.frame(matrix(nrow = nrow(exp), ncol = ncol(exp)))
for (i in 1:ncol(exp)) {
  min_value = min(exp[[i]], na.rm = TRUE)
  max_value = max(exp[[i]], na.rm = TRUE)
  normalized_exp[[i]] = (exp[[i]] - min_value) / (max_value - min_value) 
}
normalized_exp = as.data.frame(t(normalized_exp))
gene = read.table("./exp_gene.txt",
                  sep = "\t")
samples_exp = read.table("./exp_sample.txt",
                    sep = "\t")
rownames(normalized_exp) = samples_exp[,1]
colnames(normalized_exp) = gene[,1]
num = dim(normalized_exp)[1];

# split reference samples into training set and test 
trainNum = round(num*0.8);
testNum = round(num*0.2);

datTrain = as.data.frame(normalized_exp[1:trainNum,]);
dim(datTrain);

gsg_train = goodSamplesGenes(datTrain, verbose=3);
gsg_train$allOK;

badSamples_train = which(!gsg_train$goodSamples)
cat("Train badSamples:", badSamples_train, "\n")

badGenes_train = which(!gsg_train$goodGenes)
cat("Train badGenes:", badGenes_train, "\n")

datTrain = datTrain[gsg_train$goodSamples, gsg_train$goodGenes]


datTest = as.data.frame(normalized_exp[trainNum+1:testNum,]);
dim(datTest)

gsg_test = goodSamplesGenes(datTest, verbose=3);
gsg_test$allOK;

badSamples_test = which(!gsg_test$goodSamples)
cat("Test badSamples:", badSamples_test, "\n")

badGenes_test = which(!gsg_test$goodGenes)
cat("Test badGenes:", badGenes_test, "\n")

datTest = datTest[gsg_test$goodSamples, gsg_test$goodGenes]

common_cols = intersect(colnames(datTrain), colnames(datTest))
datTrain = datTrain[, common_cols]
datTest = datTest[, common_cols]
normalized_exp = normalized_exp[, common_cols]

write.table(normalized_exp, file="./exp_data_normalized.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=TRUE)
write.table(colnames(normalized_exp), file="./exp_gene_normalized.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)


# ========================================================================
# determine parameters (Soft threshold)
powers = c(c(1:10), seq(from=12, to=20, by=2))
sft = pickSoftThreshold(datTrain, powerVector = powers, verbose=5)

sizeGrWindow(9,5)
par(mfrow=c(1,2))
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

beta = sft$fitIndices[,];


# ========================================================================
# calculate adjacency matrix
adjacency = adjacency(datTrain, power=8); # power需要按照上面的图进行调试，使得拟合和连通性都较好

# calculate weight matrix
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM;

# save weight matrix
write.csv(TOM, file="./df_tom_similarity_none.csv", row.names = FALSE)


# ========================================================================
# detect co-expressed modules
geneTree = hclust(as.dist(dissTOM), method = "average");
minModuleSize = 50;
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);

table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)

table(dynamicColors)

sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

genes = colnames(datTrain)
df_modules = data.frame(genes, dynamicColors)

# save module detection results
write.csv(df_modules, file = "./df_modules.csv")


# ========================================================================
# test module preservation using test set
setLabels = c("Train", "Test");
multiExpr = list(Train = list(data = datTrain), Test = list(data = datTest));
multiColor = list(Train = dynamicColors);

system.time( {
  mp = modulePreservation(multiExpr, multiColor,
                          referenceNetworks = 1,
                          nPermutations = 200,
                          randomSeed = 1,
                          quickCor = 0,
                          verbose = 3)
} );


ref = 1
test = 2

statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])

statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);

temp = cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2))

# save module preservation test results
write.csv(temp, file="./df_zsummary.csv")

print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

modColors = rownames(mp$preservation$observed[[ref]][[test]])

moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];

# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
text = modColors[plotMods];
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
mains = c("Preservation Median rank", "Preservation Zsummary");

sizeGrWindow(10, 5);
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  #labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08);
  # For Zsummary, add threshold lines
  if (p==2)
  {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=15, col = "darkgreen", lty = 2)
  }
}