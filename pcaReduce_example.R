b <-read.table("Pollen2014.txt", sep=',', header = T,row.names=1)

lb <-read.table("SupplementaryLabels.txt", sep=',', header = T)

D <- log2(as.matrix(b) + 1) # log transformation of count data

Input <- t(D) # data matrix, cells in rows, genes in columns

#====REPRODUCE PAPER FIGURE 3A==================================
true_tissue_cls <- lb[,6] # true data partition K=4
true_cell_cls <- lb[,4]   # true data partition K=11

library("pcaMethods")
Y <- prep(Input, scale="none", center=TRUE)
pca_out <- pca(Y, method="svd", center = FALSE, nPcs=2)
x <- pca_out@scores
l <- pca_out@loadings

pdf("Paper_figure2A.pdf")
C <- c("skyblue", "mediumvioletred", " cornflowerblue", "olivedrab4", "olivedrab3", "aquamarine4", "salmon1",  "black", "orangered3",   "royalblue3","darkseagreen3"  )
plot(x[which(lb[,3]=="K562"),1], x[which(lb[,3]=="K562"),2], type="p",pch=21, bg = "salmon1", col= "salmon2", cex=1.7, xlab="PCA1", ylab="PCA2",ylim=c(-150,200), xlim=c(-220,320),  cex.lab=1.5, cex.axis=1.7, font.lab=2, bty='n')
points(x[which(lb[,3]=="HL60"),1], x[which(lb[,3]=="HL60"),2], type="p",pch=21, bg = "orangered3", col="orangered4",cex=1.7)
points(x[which(lb[,3]=="2339"),1], x[which(lb[,3]=="2339"),2], type="p",pch=21, bg = "plum3", col = "plum4", cex=1.7)

points(x[which(lb[,3]=="iPS"),1], x[which(lb[,3]=="iPS"),2], type="p",pch=21, bg = "gray61", col = "gray46", cex=1.7)

points(x[which(lb[,3]=="2338"),1], x[which(lb[,3]=="2338"),2], type="p",pch=21, bg = "steelblue2", col= "steelblue3", cex=1.7)
points(x[which(lb[,3]=="BJ"),1], x[which(lb[,3]=="BJ"),2], type="p",pch=21, bg = "skyblue1", col= "skyblue2", cex=1.7)

points(x[which(lb[,3]=="NPC"),1], x[which(lb[,3]=="NPC"),2], type="p",pch=21, bg = "olivedrab2", col= "olivedrab3", cex=1.7)
points(x[which(lb[,3]=="GW21"),1], x[which(lb[,3]=="GW21"),2], type="p", pch=21, bg = "darkolivegreen4", col= "darkolivegreen", cex=1.7)

points(x[which(lb[,3]=="GW21+3"),1], x[which(lb[,3]=="GW21+3"),2], type="p", pch=21, bg = "darkseagreen", col= "darkseagreen4", cex=1.7)

points(x[which(lb[,3]=="Kera"),1], x[which(lb[,3]=="Kera"),2], type="p",pch=3, col = "royalblue3", cex=2)
points(x[which(lb[,3]=="GW16"),1], x[which(lb[,3]=="GW16"),2], type="p", pch=4, col = "aquamarine3", cex=2)

temp <- legend("topright",  border = NULL,legend = c("K562","HL60", "2339","", "iPS"), col = c("salmon2", "orangered4", "plum4","white","gray46"),
text.width = strwidth("00,00000,0000000"),
lty = c(0), pch=c(21,21,21,21,21), pt.bg=c("salmon1","orangered3","plum3","white","gray61"), xjust = 0, yjust = 0, bty = "n")


temp <- legend("topright",  border = NULL,legend = c("2338","BJ", "Kera", "", "NPC", "GW21", "GW21+3", "GW16"), col = c("steelblue3", "skyblue2", "royalblue3", "white", "olivedrab3", "darkolivegreen", "darkseagreen4", "aquamarine3" ),
text.width = strwidth("00,000,0"),
lty = c(0), pch=c(21,21,3,3,21,21,21,4), pt.bg=c("steelblue2","skyblue2", "white", "white", "olivedrab2", "darkolivegreen4", "darkseagreen"), xjust = 0, yjust = 0,bty = "n")
dev.off()


#====RUN pcaReduce==============================================
library("pcaReduce")

Output_S <- PCAreduce(Input, nbt=100, q=30, method='S')
# will produce a list, where each element in the list is a matrix/sample
# (i.e. we will have 100 runs of pcaReduce algorithm, where merging was achieved based on sampling);
# each row in each sample matrix corresponds to a clustering; first row contains data partition into K=q+1 clusters, and last row will contain K=2 clusters.

Output_M <- PCAreduce(Input, nbt=100, q=30, method='M')
# output is similar list, however the merging was achieved based on largest probability value.

#====PRODUCE PAPER FIGURE 3B==================================

N <- length(Output_S)
M <- dim(Output_S[[1]])[2]

K11 <- c()
K4 <- c()
for (n in 1:N){
    cls_cell <- c()
    cls_tissue <- c()
    labels <- c()
    
    for (m in 1:M){
        cls_cell <- c(cls_cell, adjustedRandIndex(Output_S[[n]][,m], true_cell_cls))
        cls_tissue <- c(cls_tissue, adjustedRandIndex(Output_S[[n]][,m], true_tissue_cls))
        labels <- c(labels, length(unique(Output_S[[n]][,m])))
    }
    
    K11 <- cbind(K11, cls_cell)
    K4 <- cbind(K4, cls_tissue)
}

pdf("Paper_figure2B.pdf")
plot(labels, K11[,1], col="cornflowerblue",type="l", lty=3, lwd=0.5, main="", xlab="Number of clusters", ylab="ARANDI score", ylim=c(0,1), cex.lab=1.7, cex.axis=1.5, font.lab=2, bty='n')
for (i in 1:N){
    lines(labels, K11[,i], col="cornflowerblue", type="l", lty=3, lwd=1.3)
}
for (i in 1:N){
    lines(labels, K4[,i],  col="yellowgreen",type="l", lty=3, lwd=1.3)
}
for (i in 1:N){
    points(labels[which(labels==11)], K11[which(labels==11),i], cex=1.4, pch=21, bg="cornflowerblue", col="royalblue3")
    points(labels[which(labels==4)],  K4[which(labels==4),i],   cex=1.4, pch=21, bg="yellowgreen", col="olivedrab4" )
}
temp <- legend("topright",  border = NULL,legend = c("K = 4", "K = 11"), col = c("olivedrab4", "royalblue3"),
text.width = strwidth("00,000"),
lty = c(0,0), pch=c(21,21), pt.bg=c("yellowgreen", "cornflowerblue"), xjust = 0, yjust = 0, bty = "n")
dev.off()

#====ADDITIONAL FIGURE ==================================
# code below will produce figure similar to paper figure 3B
# however here the merging is based on largest probability.

N <- length(Output_M)
M <- dim(Output_M[[1]])[2]

K11 <- c()
K4 <- c()
for (n in 1:N){
    cls_cell <- c()
    cls_tissue <- c()
    labels <- c()
    
    for (m in 1:M){
        cls_cell <- c(cls_cell, adjustedRandIndex(Output_M[[n]][,m], true_cell_cls))
        cls_tissue <- c(cls_tissue, adjustedRandIndex(Output_M[[n]][,m], true_tissue_cls))
        labels <- c(labels, length(unique(Output_M[[n]][,m])))
    }
    
    K11 <- cbind(K11, cls_cell)
    K4 <- cbind(K4, cls_tissue)
}

pdf("Additional_figure.pdf")
plot(labels, K11[,1], col="lightsalmon2",type="l", lty=3, lwd=0.5, main="", xlab="Number of cluster", ylab="ARANDI score", ylim=c(0,1), cex.lab=1.7, cex.axis=1.5, font.lab=2, bty='n')
for (i in 1:N){
    lines(labels, K11[,i], col="lightsalmon2", type="l", lty=3, lwd=1.3)
}
for (i in 1:N){
    lines(labels, K4[,i],  col="yellow3",type="l", lty=3, lwd=1.3)
}
for (i in 1:N){
    points(labels[which(labels==11)], K11[which(labels==11),i], cex=1.4, pch=21, bg="lightsalmon2", col="lightsalmon3")
    points(labels[which(labels==4)],  K4[which(labels==4),i],   cex=1.4, pch=21, bg="yellow3", col="wheat3" )
}
temp <- legend("topright",  border = NULL,legend = c("K = 4", "K = 11"), col = c("wheat3", "lightsalmon3"),
text.width = strwidth("00,000"),
lty = c(0,0), pch=c(21,21), pt.bg=c("yellow3", "lightsalmon2"), xjust = 0, yjust = 0, bty = "n")
dev.off()
