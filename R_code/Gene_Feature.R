setwd("")

library('assert')
library('igraph')

read_snv_matrix = function(cancer) {
  
  # snv matrix
  snv_matrix = read.table(paste0("./",cancer,"_snv_matrix.txt"))
  snv_sample = read.table(paste0("./",cancer,"_snv_sample.txt"))[,1]
  snv_gene = read.table(paste0("./",cancer,"_snv_gene.txt"))[,1]
  colnames(snv_matrix) = snv_sample
  rownames(snv_matrix) = snv_gene
  
  return(snv_matrix)
}

read_cnv_matrix = function(cancer) {
  
  # cnv data
  cnv_matrix = read.table(paste0("./",cancer,"_cnv_matrix.txt"))
  cnv_sample = read.table(paste0("./",cancer,"_cnv_sample.txt"))[,1]
  cnv_gene = read.table(paste0("./",cancer,"_cnv_gene.txt"))[,1]
  colnames(cnv_matrix) = cnv_sample
  rownames(cnv_matrix) = cnv_gene
  
  return(cnv_matrix)
}

read_mut_matrix = function(cancer) {
  
  # mut matrix
  mut_matrix = read.table(paste0("./",cancer,"/",cancer,"_mutation.txt"))
  mut_sample = read.table(paste0("./",cancer,"/",cancer,"_sample.txt"))[,1]
  mut_gene = read.table(paste0("./",cancer,"/",cancer,"_gene.txt"))[,1]
  colnames(mut_matrix) = mut_sample
  rownames(mut_matrix) = mut_gene
  
  return(mut_matrix)
}

read_tumor_exp = function(cancer) {
  
  # tumor exp data
  exp_data = read.table(paste0("./",cancer,"/",cancer,"_exp_data.txt"))
  exp_sample = read.table(paste0("./",cancer,"/",cancer,"_sample.txt"))[,1]
  exp_gene = read.table(paste0("./",cancer,"/",cancer,"_gene.txt"))[,1]
  colnames(exp_data) = exp_sample
  rownames(exp_data) = exp_gene
  
  return(exp_data)
}

read_tumor_met = function(cancer) {
  
  # tumor met data
  met_data = read.table(paste0("./",cancer,"/",cancer,"_met_data.txt"),sep = "\t")
  met_sample = read.table(paste0("./",cancer,"/",cancer,"_sample.txt"))[,1]
  met_gene = read.table(paste0("./",cancer,"/",cancer,"_gene.txt"))[,1]
  colnames(met_data) = met_sample
  rownames(met_data) = met_gene
  
  return(met_data)
}

read_normal_exp = function(cancer) {
  
  if (cancer == "TCGA-BRCA") {
    normal_exp = read.table(paste0("./",cancer,"/",cancer,"_exp_normal.txt"),
                             sep = "\t")
    normal_sample = read.table(paste0("./",cancer,"/",cancer,"_exp_sample.txt"),
                                sep = "\t")[,1]
    normal_gene = read.table(paste0("./",cancer,"/",cancer,"_gene.txt"),
                              sep = "\t")[,1]
  } else {
    normal_exp = read.table(paste0("./",cancer,"/",cancer,"_exp_normal.txt"),
                             sep = "\t")
    normal_sample = read.table(paste0("./",cancer,"/",cancer,"_sample.txt"),
                                sep = "\t")[,1]
    normal_gene = read.table(paste0("./",cancer,"/",cancer,"_gene.txt"),
                              sep = "\t")[,1]
  }
  rownames(normal_exp) = normal_gene
  
  return(normal_exp)
}

read_normal_met = function(cancer) {
  
  if (cancer == "TCGA-BRCA") {
    normal_met = read.table(paste0("./",cancer,"/",cancer,"_met_normal.txt"),
                             sep = "\t")
    normal_sample = read.table(paste0("./",cancer,"/",cancer,"_met_sample.txt"),
                                sep = "\t")[,1]
    normal_gene = read.table(paste0("./",cancer,"/",cancer,"_gene.txt"),
                              sep = "\t")[,1]
  } else {
    normal_met = read.table(paste0("./",cancer,"/",cancer,"_met_normal.txt"),
                             sep = "\t")
    normal_sample = read.table(paste0("./",cancer,"/",cancer,"_sample.txt"),
                                sep = "\t")[,1]
    normal_gene = read.table(paste0("./",cancer,"/",cancer,"_gene.txt"),
                              sep = "\t")[,1]
  }
  rownames(normal_met) = normal_gene
  
  return(normal_met)
}

read_filter_PPI = function(PPI_file,seita,is_filtered=F,if_write=T) {
  
  PPI = read.table(PPI_file,header = T)
  assert(var(PPI[,3]) > 0)
  if (!is_filtered) {
    threshold = as.numeric(quantile(PPI[,3],seita)) # seita的值表示quantile取多少分位数
    PPI = PPI[which(PPI[,3] >= threshold),]
    if (if_write) {
      write.table(PPI,paste0("./",PPI_name,"_",seita,".txt"),sep = "\t",
                  quote = F,col.names = T,row.names = F)
    }
  }
  
  return(PPI)
}

get_feature_in_tumor = function(cancer,snv_matrix,cnv_matrix,exp_data,met_data,
                                 normal_exp,normal_met,sample_cluster,if_write=T) {
  # gene feature
  cluster_num = max(sample_cluster[,1])
  genes_g = rownames(exp_data)
  sample_g = colnames(exp_data)
  ave_exp = apply(normal_exp,1,mean)
  ave_exp = data.frame(ave_exp)
  ave_met = apply(normal_met,1,mean)
  ave_met = data.frame(ave_met)
  gene_feature = list()
  
  # assert
  assert(sum(rownames(exp_data) == rownames(met_data)) == dim(exp_data)[1],
         sum(rownames(exp_data) == rownames(normal_exp)) == dim(exp_data)[1],
         sum(rownames(exp_data) == rownames(normal_met)) == dim(exp_data)[1])
  
  # feature for cluster_k
  for (k in 1:cluster_num) {
    cluster_k = rownames(sample_cluster)[which(sample_cluster[,1] == k)]
    # snv
    snv_cluster_k = intersect(colnames(snv_matrix),cluster_k)
    snv_ratio_k = rowSums(snv_matrix[,snv_cluster_k]) / length(cluster_k)
    snv_ratio_k = as.matrix(snv_ratio_k)
    # cnv
    cnv_cluster_k = intersect(colnames(cnv_matrix),cluster_k)
    cnv_ratio_k = rowSums(abs(cnv_matrix[,cnv_cluster_k])) / length(cluster_k)
    cnv_ratio_k = as.matrix(cnv_ratio_k)
    # exp
    exp_cluster_k = apply(exp_data[,cluster_k],1,mean) - ave_exp
    # met
    met_cluster_k = apply(met_data[,cluster_k],1,mean) - ave_met
    # union feature_k
    gene_feature_k = matrix(0,nrow = length(genes_g),ncol = 2)
    rownames(gene_feature_k) = genes_g
    gene_feature_k[rownames(snv_ratio_k),1] = snv_ratio_k
    gene_feature_k[rownames(cnv_ratio_k),2] = cnv_ratio_k
    gene_feature_k = cbind(gene_feature_k,exp_cluster_k,met_cluster_k)
    colnames(gene_feature_k) = paste0(c("snv_ratio","cnv_ratio","exp_diff","met_diff"),"_",k)
    gene_feature[[k]] = gene_feature_k
  }
  
  # feature for global
  # snv
  snv_ratio_g = rowSums(snv_matrix) / length(sample_g)
  snv_ratio_g = as.matrix(snv_ratio_g)
  # cnv
  cnv_ratio_g = rowSums(abs(cnv_matrix)) / length(sample_g)
  cnv_ratio_g = as.matrix(cnv_ratio_g)
  # exp
  exp_cluster_g = apply(exp_data,1,mean) - ave_exp
  # met
  met_cluster_g = apply(met_data,1,mean) - ave_met
  # union feature_g
  gene_feature_g = matrix(0,nrow = length(genes_g),ncol = 2)
  rownames(gene_feature_g) = genes_g
  gene_feature_g[rownames(snv_ratio_g),1] = snv_ratio_g
  gene_feature_g[rownames(cnv_ratio_g),2] = cnv_ratio_g
  gene_feature_g = cbind(gene_feature_g,exp_cluster_g,met_cluster_g)
  colnames(gene_feature_g) = paste0(c("snv_ratio","cnv_ratio","exp_diff","met_diff"),"_g")
  
  # union cluster global
  gene_feature[[k + 1]] = gene_feature_g
  gene_feature = data.frame(gene_feature)
  feature_name = colnames(gene_feature)
  
  # if write
  if (if_write) {
    write.table(gene_feature,paste0("./",cancer,"/",cancer,"_feature_tumor.txt"),
                sep = "\t",quote = F,col.names = F,row.names = F)
    write.table(feature_name,paste0("./",cancer,"/",cancer,"_feature_name_tumor.txt"),
                sep = "\t",quote = F,col.names = F,row.names = F)
    write.table(genes_g,paste0("./",cancer,"/",cancer,"_gene.txt"),
                sep = "\t",quote = F,col.names = F,row.names = F)
  }
  
  return(gene_feature)
}

calculate_clustering_coefficient = function(Graph,com_gene,delta_k,delta_m) {
  
  # clustering coefficient
  assert(delta_k > 0)
  clu_coef = c()
  for (i in 1:length(com_gene)) {
    # local network
    neighbor_gene = c(com_gene[i])
    count_num = 1
    for (k in 1:delta_k) {
      len_num = length(neighbor_gene)
      if (count_num > len_num) {
        break
      }
      for (j in count_num:len_num) {
        tmp_gene = neighbor_gene[j]
        neighbor_k_j = rownames(data.frame(Graph[[tmp_gene]]))
        neighbor_gene = unique(c(neighbor_gene,neighbor_k_j))
      }
      count_num = j + 1
    }
    neighbor_gene = neighbor_gene[-1]
    Graph_local = subgraph(Graph,neighbor_gene)
    
    # calculate clustering coefficient
    g_local = gorder(Graph_local)
    e_local = gsize(Graph_local)
    if (g_local >= delta_m) {
      clu_coef = c(clu_coef, 2 * e_local / (g_local * (g_local - 1)))
    } else {
      clu_coef = c(clu_coef, 0)
    }
  }
  names(clu_coef) = com_gene
  
  return(clu_coef)
}

get_feature_in_PPI = function(cancer,mut_gene,PPI,PPI_name,if_write=T) {
  
  # gene feature
  gene_feature = matrix(0,nrow = length(mut_gene),ncol = 4)
  colnames(gene_feature) = c("deg_cent","clu_coef","clo_cent","bet_cent")
  rownames(gene_feature) = mut_gene
  
  # Graph
  Graph = graph_from_data_frame(PPI)
  Graph = as_undirected(Graph,mode ="collapse", edge.attr.comb ="sum")
  gene_degree = degree(Graph)
  com_gene = intersect(mut_gene,names(V(Graph)))
  
  # degree centrality
  max_degree = max(gene_degree)
  gene_feature[com_gene,"deg_cent"] = gene_degree[com_gene] / max_degree
  
  # clustering coefficient
  clu_coef = calculate_clustering_coefficient(Graph,com_gene,delta_k = 1,delta_m = 2)
  gene_feature[com_gene,"clu_coef"] = clu_coef
  
  # closeness centrality
  clo_cent = closeness(Graph,vids = V(Graph),weights = NULL,normalized = T)
  gene_feature[com_gene,"clo_cent"] = clo_cent[com_gene]
  
  # betweenness centrality
  bet_cent = betweenness(Graph,v = V(Graph),directed = F,weights = NULL,normalized = T)
  gene_feature[com_gene,"bet_cent"] = bet_cent[com_gene]
  feature_name = colnames(gene_feature)
  
  # if write
  if (if_write) {
    write.table(gene_feature,paste0("./",cancer,"/",cancer,"_feature_PPI.txt"),
                sep = "\t",quote = F,col.names = F,row.names = F)
    write.table(feature_name,paste0("./",cancer,"/",cancer,"_feature_name_PPI.txt"),
                sep = "\t",quote = F,col.names = F,row.names = F)
  }
  
  return(gene_feature)
}

union_feature = function(cancer,feature_in_tumor,feature_in_PPI,PPI_name) {
  
  # assert
  genes_g = rownames(feature_in_tumor)
  assert(sum(genes_g == rownames(feature_in_PPI)) == length(genes_g))
  
  # union feature in tumor PPI
  gene_feature = cbind(feature_in_tumor,feature_in_PPI)
  feature_name = colnames(gene_feature)
  write.table(gene_feature,paste0("./",cancer,"/",cancer,"_feature_",PPI_name,".txt"),
              sep = "\t",quote = F,col.names = F,row.names = F)
  write.table(feature_name,paste0("./",cancer,"/",cancer,"_feature_name_",PPI_name,".txt"),
              sep = "\t",quote = F,col.names = F,row.names = F)
  write.table(genes_g,paste0("./",cancer,"/",cancer,"_gene_",PPI_name,".txt"),
              sep = "\t",quote = F,col.names = F,row.names = F)
  
  return(cancer)
}

reshape_data = function(cancer,data,gene,data_name) {
  
  data = data[gene,]
  sample_id = colnames(data)
  out_path = paste0("./",cancer,"/",data_name,"/")
  dir.create(out_path)
  write.table(data,paste0(out_path,cancer,"_",data_name,".txt"),
              sep = "\t",quote = F,col.names = F,row.names = F)
  write.table(sample_id,paste0(out_path,cancer,"_sample.txt"),
              sep = "\t",quote = F,col.names = F,row.names = F)
  write.table(gene,paste0(out_path,cancer,"_gene.txt"),
              sep = "\t",quote = F,col.names = F,row.names = F)
  
  return(data)
}


# weighted PPI network
seita = 0.7
PPI_file = ""
PPI = read_filter_PPI(PPI_file,seita,is_filtered = F,if_write = T)

# cancer with normal >= 20
cancer_types = c("TCGA-BLCA", "TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", 
                  "TCGA-KIRC", "TCGA-KIRP", "TCGA-LUAD", "TCGA-LUSC", 
                  "TCGA-PRAD", "TCGA-UCEC")

for (i in 1:length(cancer_types)) {
  cancer = cancer_types[i]
  dir.create(paste0("./",cancer,"/"))
  print(paste0("Create gene features for ",cancer,"..."))
  
  # data
  snv_matrix = read_snv_matrix(cancer)
  cnv_matrix = read_cnv_matrix(cancer)
  exp_data = read_tumor_exp(cancer)
  met_data = read_tumor_met(cancer)
  normal_exp = read_normal_exp(cancer)
  normal_met = read_normal_met(cancer)
  normal_met[is.na(normal_met)] <- 0 
  sample_cluster = read.table(paste0("./",cancer,"/community_cohesion_scores/cluster_results.txt"))
  rownames(sample_cluster)  =  sample_cluster[,2]
  mut_matrix = read_mut_matrix(cancer)
  
  # gene feature in tumor
  feature_in_tumor = get_feature_in_tumor(cancer,snv_matrix,cnv_matrix,exp_data,met_data,
                                           normal_exp,normal_met,sample_cluster,if_write = F)
  
  # gene feature in PPI
  mut_gene = rownames(feature_in_tumor)
  feature_in_PPI = get_feature_in_PPI(cancer,mut_gene,PPI,PPI_name,if_write = F)
  
  # union feature
  a = union_feature(cancer,feature_in_tumor,feature_in_PPI,PPI_name)
  print(a)
  write.table(seita,paste0("./",cancer,"/",cancer,"_",PPI_name,"_",seita,".txt"))
  print(paste0("------Store gene features in dir: ./",cancer,"/"))
  
  # reshape data
  exp_data = reshape_data(cancer,exp_data,mut_gene,"tumor_exp")
  normal_exp = reshape_data(cancer,normal_exp,mut_gene,"normal_exp")
  met_data = reshape_data(cancer,met_data,mut_gene,"tumor_met")
  normal_met = reshape_data(cancer,normal_met,mut_gene,"normal_met")
  mut_matrix = reshape_data(cancer,mut_matrix,mut_gene,"mutation")
}
