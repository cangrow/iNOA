library(docopt)
"Usage: case_control_compare_deseq2.r [options] INPUT CASE_GROUP CONTROL_GROUP CASE CONTROL LOG2FC PVALUE OUTPUT
Options:
   -w --width=width    the width of viewport  [default: 12]
   -h --height=width   the height of viewport [default: 12]
Arguments:
  INPUT          the input file name
  CASE_GROUP     the case group name
  CONTROL_GROUP  the control group name
  CASE           the case group
  CONTROL        the control group
  LOG2FC         the log2 fold change
  PVALUE         the pvalue
  OUTPUT         the output directory name" -> doc


opts          <- docopt(doc)
input         <- opts$INPUT
case_group    <- opts$CASE_GROUP
control_group <- opts$CONTROL_GROUP
case          <- opts$CASE
control       <- opts$CONTROL
log2fc        <- as.numeric(opts$LOG2FC)
p_value       <- as.numeric(opts$PVALUE)
output        <- opts$OUTPUT



library(DESeq2)

x <- read.table(input, header=T, sep="\t", row.names=1, comment.char="", check.names=F)

control_count           <- sapply(unlist(strsplit(control, split=",")), function(t){x[[t]] })
str(control_count)
colnames(control_count) <- unlist(strsplit(control, split=","))

case_count              <- sapply(unlist(strsplit(case, split=",")), function(t){ x[[t]] })
str(case_count)
colnames(case_count)    <- unlist(strsplit(case, split=","))

count                   <- cbind(data.frame(case_count, check.names=F), data.frame(control_count, check.names=F))
rownames(count)         <- rownames(x)


countData <- count[apply(count, 1, sum) > 0 , ]

group = factor(c(rep(case_group, length(colnames(case_count))), rep(control_group, length(colnames(control_count)))  ), levels = c(control_group, case_group))  # levles， control要放在前面，保证后续差异分析时，以control为对照
names(group) = colnames(countData)
colData = data.frame(condition = group)
deseq_formula = as.formula('~ condition')

dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = deseq_formula)
dds <- DESeq(dds,betaPrior=FALSE)
cnt <- as.data.frame(counts(dds, normalized=TRUE))

# norm
norm_file <- paste(output, "/diff_norm.txt", sep="")
# write.table(cnt, norm_file, sep="\t", quote=F,  col.names = NA)
# tmp = data.frame(ID = rownames(cnt), cnt)  ### 纯数字样本名会加上X前缀
tmp = data.frame(ID = rownames(cnt))
tmp = cbind(tmp, cnt)
write.table(tmp, norm_file, sep="\t", quote=F, col.names = T, row.names = F)

# result
sample_num <- ncol(cnt)
cnt$baseMeanA <- apply( sapply(unlist(strsplit(case, split=",")), function(t){cnt[[t]]}), 1, mean )
cnt$baseMeanB <- apply( sapply(unlist(strsplit(control, split=",")), function(t){cnt[[t]]}), 1, mean )
res <- as.data.frame(results(dds))
res <- cbind(cnt, res)
res$type <- "Not DEG"
res$type[res$pvalue < p_value & res$log2FoldChange >= log2fc ] <- "Up"
res$type[res$pvalue < p_value & res$log2FoldChange <= -(log2fc)] <- "Down"
res$type <- factor(res$type, levels = c("Up", "Down", "Not DEG"))
diff_file <- paste(output, "/diff.txt", sep="")
# write.table(res, diff_file, sep="\t", quote=F,  col.names = NA)
# tmp1 = data.frame(ID = rownames(res), res)
tmp1 = data.frame(ID = rownames(res))
tmp1 = cbind(tmp1, res)
write.table(tmp1, diff_file, sep="\t", quote=F, col.names = T, row.names = F)

# correlation plot
library(ggplot2)

correlation_file <- paste(output, "/correlation.pdf", sep="")
pdf(correlation_file)
ggplot(res, aes(x=log10(baseMeanA+ 0.00000000001),y= log10(baseMeanB+0.00000000001))) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, color="red") + 
  scale_x_continuous(limits=c(-2.5, 5)) + 
  scale_y_continuous(limits=c(-2.5, 5)) + 
  labs(x=case_group, y=control_group)
dev.off()

# MA plot
ma_file <- paste(output, "/ma.pdf", sep="")
pdf(ma_file)
ggplot(res, aes(x = baseMean, y = log2FoldChange , colour = type)) + 
  geom_point() + 
  scale_x_continuous(limits = c(0, 2e+04)) + 
  scale_y_continuous(limits = c(-20, 20)) + 
  theme(legend.title = element_blank()) + 
  labs(x = "Mean Expression", y="log2(FC)", tilte="MA plot")
dev.off()

# valcano plot
library(ggrepel)

up_gene <- res[which(res[, "type"] == "Up"), ]
down_gene <- res[which(res[, "type"] == "Down"), ]
if(nrow(up_gene) < 5)
{
	top5_up_gene <- up_gene
}else
{
	up_gene <- up_gene[order(up_gene$pvalue), ]
	top5_up_gene <- up_gene[1:5, ]
}

if(nrow(down_gene) < 5)
{
	top5_down_gene <- down_gene
}else
{
	down_gene <- down_gene[order(down_gene$pvalue), ]
	top5_down_gene <- down_gene[1:5, ]
}
top_gene <- rbind(top5_up_gene, top5_down_gene)
valcano_file <- paste(output, "/valcano.pdf", sep="")
pdf(valcano_file)
ggplot(res) + 
  geom_point(aes(x = log2FoldChange, y = -log10(res$pvalue), color =type)) +
  geom_text_repel(data = top_gene, aes(x = top_gene$log2FoldChange, y = -log10(top_gene$pvalue), label=rownames(top_gene)), segment.color = "black", segment.size = 0.5, force = 1, box.padding = unit(1.0, "lines"), point.padding = unit(0.5, "lines"), segment.alpha = 1, show.legend = FALSE, max.overlaps = 10000) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = bquote(paste(log[2],"(fold change)",sep="")), y = bquote(paste(-log[10],"(p value)",sep="")), title = paste("valcano(",case_group,"/",control_group,")",sep="")) +
  scale_x_continuous(limits=c(-10,10)) +
  geom_vline(xintercept = c(-log2fc, log2fc), lty = 4, col = "grey", lwd = 0.8) +
  geom_hline(yintercept = -log10(p_value),lty = 4, col = "grey", lwd = 0.8) +
  annotate("text", x = (log2fc+10)/2, y = Inf, label = sum(res$type=="Up"), color = "#f8766d", vjust = 2) +
  annotate("text", x = (log2fc-10)/2, y = Inf, label = sum(res$type=="Down"), color = "#00ba38", vjust = 2)
dev.off()

# heatmap plot
library(pheatmap)
myheatcol = colorRampPalette(c('green','black','red'))(100)
data <- data.matrix(res[res$type != "Not DEG", 1:sample_num])

# 样本数大于50，字号调整
fontsize_row <- 10
if(sample_num > 50){fontsize_row <- 500/sample_num}

if(nrow(data) > 1)
{
    # 基因水平聚类
    heatmap_file <- paste(output, "/heatmap.pdf", sep="")
    pdf(heatmap_file , width=as.numeric(opts$w), height = as.numeric(opts$h) )
    pheatmap(data,
             scale = 'row',
             margins = c(8,10),
             cluster_rows = T,
             cluster_cols = F,
             color = myheatcol,
             fontsize_row = fontsize_row,
             show_rownames = F,
    )
    dev.off()

    # 基因+样本水平聚类
    heatmap_file <- paste(output, "/heatmap_sample_cluster.pdf", sep="")
    pdf(heatmap_file , width=as.numeric(opts$w), height = as.numeric(opts$h) )
    pheatmap(data,
             scale = 'row',
             margins = c(8,10),
             cluster_rows = T,
             cluster_cols = T,
             color = myheatcol,
             fontsize_row = fontsize_row,
             show_rownames = F,
    )
    dev.off()

    # Top50 heatmap plot

    res <- res[res$type != 'Not DEG', ]
    if( nrow(res) < 50 )
    {
        data_top <- data.matrix(res[ , 1:sample_num])          # 小于50,全画
    }else
    {
        res <- res[order(res$pvalue), ]
        data_top <- data.matrix(res[1:50, 1:sample_num])           #  相对丰度最高的50个
    }

    # 基因水平聚类
    heatmap_file <- paste(output, "/heatmap_top50.pdf", sep="")
    pdf(heatmap_file , width=as.numeric(opts$w), height = as.numeric(opts$h) )
    pheatmap(data_top,
             scale = 'row',
             margins = c(8,10),
             cluster_rows = T,
             cluster_cols = F,
             color = myheatcol,
             fontsize_row = fontsize_row,
    )
    dev.off()

    # 基因+样本水平聚类
    heatmap_file <- paste(output, "/heatmap_top50_sample_cluster.pdf", sep="")
    pdf(heatmap_file , width=as.numeric(opts$w), height = as.numeric(opts$h) )
    pheatmap(data_top,
             scale = 'row',
             margins = c(8,10),
             cluster_rows = T,
             cluster_cols = T,
             color = myheatcol,
             fontsize_row = fontsize_row,
    )
    dev.off()
}else
{
    message("heatmap skipped. 差异表达量不足1个")
}
