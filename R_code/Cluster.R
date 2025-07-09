library(SNFtool)
data1 = read.csv("./community_cohesion_scores/df_ccs.csv")
rownames(data1) = data1[,1]
data1 = data1[, -c(1)]
colnames(data1) = gsub("\\.", "-", colnames(data1))

data2 = read.table("./exp_data.txt", sep = "\t")
gene2 = read.table("./exp_gene.txt",
                  sep = "\t")
samples_exp = read.table("./exp_sample.txt",
                         sep = "\t")
colnames(data2) = samples_exp[,1]
rownames(data2) = gene2[,1]

data3 = read.table("./met_data.txt", sep = "\t")
gene3 = read.table("./met_gene.txt",
                   sep = "\t")
samples_met = read.table("./met_sample.txt",
                         sep = "\t")
colnames(data3) = samples_met[,1]
rownames(data3) = gene3[,1]


gene_variances_exp <- apply(data2, 1, var)
top_genes_exp <- head(order(gene_variances_exp, decreasing = TRUE), 30)
top_gene_names_exp <- rownames(data2)[top_genes_exp]
data2 = data2[top_gene_names_exp,]

gene_variances_met <- apply(data3, 1, var)
top_genes_met <- head(order(gene_variances_met, decreasing = TRUE), 30)
top_gene_names_met <- rownames(data3)[top_genes_met]
data3 = data3[top_gene_names_met,]

data1_normalized = standardNormalization(data1)
data2_normalized = standardNormalization(data2)
data3_normalized = standardNormalization(data3)

data1_normalized = as.data.frame(t(data1_normalized))
data2_normalized = as.data.frame(t(data2_normalized))
data3_normalized = as.data.frame(t(data3_normalized))


Dist1 = (dist2(as.matrix(data1_normalized),as.matrix(data1_normalized)))^(1/2)
Dist2 = (dist2(as.matrix(data2_normalized),as.matrix(data2_normalized)))^(1/2)
Dist3 = (dist2(as.matrix(data3_normalized),as.matrix(data3_normalized)))^(1/2)

K = 
alpha = 
TM = 
W1 = affinityMatrix(Dist1, K, alpha)
W2 = affinityMatrix(Dist2, K, alpha)
W3 = affinityMatrix(Dist3, K, alpha)
W = SNF(list(W1,W2,W3), K, TM)

clusters = spectralClustering(W, K = )
clusters = as.data.frame(clusters)
rownames(clusters) = rownames(Dist1)
write.csv(clusters, file="./cluster_results_survival_analysis.csv")

