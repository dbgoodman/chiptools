setwd('/Users/dbgoodman/Dropbox/chiptools/receq')
df_rc <- as.data.frame(read.csv('grid.csv', header=TRUE, 
    row.names=1, col.names=1))
df_rc <- as.data.frame(lapply(df_rc, as.logical))

convert.magic(df_rc, "logical")


df_rc2 <- as.matrix(read.csv('grid.csv', header=TRUE, row.names=1, colClasses="logical"))

read.csv('grid.csv')

df_rc3 <- as.matrix(read.csv('grid.csv', header=TRUE, row.names=1))

df_rc3 <- read.csv('grid.csv', header=TRUE, row.names=1)
df_rc3[df_rc3 == "False"] <- "0"
df_rc3[df_rc3 == "True"] <- "1"
df_rc3[df_rc3 == ""] <- "NA"
df_good <- data.matrix(df_rc3)
heatmap(t(df_good), Colv=NA, scale="none", col=colorRampPalette(c("red","white","blue"))(256))

pdf('dendrogram_of_snps.pdf', height=10, width=70)
heatmap(t(df_good), Colv=NA, scale="none", col=colorRampPalette(c("red","white","blue"))(256))
dev.off()

pdf('dendrogram_of_snps.pdf', height=10, width=70)
heatmap.2(t(df_good), Colv=FALSE, scale="none", key=FALSE, trace="none", col=colorRampPalette(c("red","white","blue"))(256))
dev.off()


setwd('/Users/dbgoodman/Dropbox/chiptools/receq')
snpdf <- as.data.frame(read.csv('all_snp_data.csv', header=TRUE, row.names=1))
snpdf$VAR <- as.logical(snpdf$VAR)
snpdf$HET <- as.logical(snpdf$HET)
snpdf$AMBER <- as.logical(snpdf$AMBER)
snpdf$STRAIN <- as.character(snpdf$STRAIN)
snpdf$STRAIN_NUM <- as.numeric(snpdf$STRAIN_NUM)

# snpdf$PARENT <- as.character(snpdf$PARENT)
# snpdf$CHILD <- as.character(snpdf$CHILD)
# 
# get_PARENT_var <- function (i) {
#     parent_var <- snpdf$VAR[
#         snpdf$POS == snpdf$POS[i] & 
#         1:dim(snpdf)[1] %in% grep(snpdf$STRAIN,snpdf$PARENT[i], fixed=TRUE)]
#     ifelse(parent_var != FALSE | is.na(parent_var), parent_var, FALSE)
# }
# 
# get_CHILD_var <- function (i) {
#     child_var <- snpdf$VAR[
#         snpdf$POS == snpdf$POS[i] & 
#         1:dim(snpdf)[1] %in% grep(snpdf$STRAIN,snpdf$CHILD[i], fixed=TRUE)]
#     ifelse(child_var != FALSE | is.na(child_var), child_var, FALSE)
# }

# snpdf$IN_PARENT <-vapply(c(1:dim(snpdf)[1]),get_PARENT_var, 0)
# snpdf$IN_CHILD <-vapply(c(1:dim(snpdf)[1]),get_CHILD_var, 0)

snpdf$LABEL <- rep('',dim(snpdf)[1])
snpdf$LABEL[snpdf$VAR] <- 'NONE'
snpdf$LABEL[snpdf$EFF_SEV == 'MODERATE'] <- 'MODERATE'
snpdf$LABEL[snpdf$EFF_SEV == 'HIGH'] <- 'HIGH'
snpdf$LABEL[snpdf$HET] <- 'HET'
snpdf$LABEL[snpdf$AMBER] <- 'AMBER'

snpdf$ORIGIN <- rep('',dim(snpdf)[1])
snpdf$ORIGIN[snpdf$IN_CHILD != 0] <- 'TRANSFERRED'
snpdf$ORIGIN[snpdf$IN_PARENT == 0] <- 'NEW SNP'
snpdf$ORIGIN[snpdf$IN_CHILD == 0] <- 'NOT TRANSFERRED'
snpdf$ORIGIN[snpdf$IN_CHILD == 0] <- 'NEW, NOT TRANSFERRED'


snpdf$ORIGIN[snpdf$IN_CHILD != 0 & snpdf$IN_PARENT != 0] <- 'TRANSFERRED'
snpdf$ORIGIN[snpdf$IN_CHILD == 0 & snpdf$IN_PARENT != 0] <- 'OLD, LOST'
snpdf$ORIGIN[snpdf$IN_CHILD != 0 & snpdf$IN_PARENT == 0] <- 'NEW, KEPT'
snpdf$ORIGIN[snpdf$IN_CHILD == 0 & snpdf$IN_PARENT == 0] <- 'NEW, LOST'


snpdf$ORIGIN[snpdf$IN_PARENT == 1] <- 'IN_ONE_PARENT'
snpdf$ORIGIN[snpdf$IN_PARENT > 1] <- 'IN_BOTH_PARENTS'
snpdf$ORIGIN[snpdf$IN_PARENT != 0 & snpdf$IN_CHILD != 0] <- 'IN_BOTH'

snpdf$SAMPLE <- reorder(snpdf$SAMPLE, snpdf$STRAIN_NUM)

ggplot(snpdf[snpdf$HET == FALSE,], 
    aes(x=POS, y=SAMPLE)) +
    geom_point(aes(
        color=interaction(as.factor(IN_CHILD),as.factor(IN_PARENT)), 
        size=as.factor(AMBER), shape=as.factor(HET)))
        
ggplot(snpdf[snpdf$HET == FALSE,], 
    aes(x=POS, y=SAMPLE)) +
    geom_point(aes(
        color=ORIGIN, 
        size=as.factor(AMBER), shape=as.factor(HET)))

ggplot(snpdf, 
    aes(x=POS, y=SAMPLE, label=EFF_GENE)) +
    geom_point(aes(
        color=ORIGIN, 
        size=as.factor(AMBER), shape=as.factor(HET))) +
    geom_text(angle = 45, size=1)

ggsave('dendrogram_of_snps_GENE_NAME.pdf', width=200,height=20)  


ggplot(snpdf[snpdf$AMBER == FALSE & snpdf$HET == FALSE,], 
    aes(x=POS, y=SAMPLE)) +
    geom_point(aes(color=EFF_SEV))

ggsave('dendrogram_of_snps_EFF_SEV.pdf', width=100,height=10)  
    



