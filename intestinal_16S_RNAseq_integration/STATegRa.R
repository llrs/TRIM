library("STATegRa")
library("gridExtra")
library("integration")

today <- format(Sys.time(), "%Y%m%d")

A <- readRDS("TRIM.RDS")
expr <- t(A[["RNAseq"]])
colnames(expr) <- seq_len(ncol(expr))
rownames(meta_i) <- seq_len(nrow(meta_i))
B1 <- createOmicsExpressionSet(Data = expr, 
                               pData = meta_i)

microbiota <- t(A[["16S"]])
colnames(microbiota) <- seq_len(ncol(microbiota))
B2 <- createOmicsExpressionSet(Data = microbiota, 
                               pData = meta_i)
cc <- selectCommonComps(X = expr, Y = microbiota, Rmax = 3)

pdf(paste0("Figures/", today, "_STATegRa.pdf"))
discoRes <- omicsCompAnalysis(Input=list("expr" = B1, "micro" = B2), 
                              Names=c("expr", "micro"),
                              method="DISCOSCA", Rcommon=2, Rspecific=c(2, 2),
                              center=TRUE, scale=TRUE, weight=TRUE)

plotRes(object=discoRes, comps=c(1, 2), what="scores", type="common",
        combined=FALSE, block="", color="Time", shape=NULL, labels=NULL,
        background=TRUE, palette=NULL, pointSize=1, labelSize=NULL,
        axisSize=NULL, titleSize=NULL) 

plotRes(object=discoRes, comps=c(1, 2), what="scores", type="common",
        combined=FALSE, block="", color="HSCT_responder", shape=NULL, labels=NULL,
        background=TRUE, palette=NULL, pointSize=1, labelSize=NULL,
        axisSize=NULL, titleSize=NULL) 

p1 <- plotRes(object=discoRes, comps=c(1, 2), what="scores", type="individual",
              combined=FALSE, block="expr", color="Time", shape=NULL,
              labels=NULL, background=TRUE, palette=NULL, pointSize=1,
              labelSize=NULL, axisSize=NULL, titleSize=NULL)
# Associated to second block
p2 <- plotRes(object=discoRes, comps=c(1, 2), what="scores", type="individual",
              combined=FALSE, block="micro", color="Time", shape=NULL,
              labels=NULL, background=TRUE, palette=NULL, pointSize=1,
              labelSize=NULL, axisSize=NULL, titleSize=NULL)

p1
p2
p1 <- plotRes(object=discoRes, comps=c(1, 2), what="loadings", type="common",
              combined=FALSE, block="expr", color="HSCT_responder", shape=NULL,
              labels=NULL, background=TRUE, palette=NULL, pointSize=1,
              labelSize=NULL, axisSize=NULL, titleSize=NULL)
p2 <- plotRes(object=discoRes, comps=c(1, 2), what="loadings", type="individual",
              combined=FALSE, block="micro", color="HSCT_responder", shape=NULL,
              labels=NULL, background=TRUE, palette=NULL, pointSize=1,
              labelSize=NULL, axisSize=NULL, titleSize=NULL)

grid.arrange(arrangeGrob(p1+theme(legend.position="none"),
                         p2+theme(legend.position="none"), nrow=1), 
             heights=c(6/7, 1/7))

plotRes(object=discoRes, comps=c(1, 2), what="loadings", type="common",
        combined=TRUE, block="", color="HSCT_responder", shape=NULL,
        labels=NULL, background=TRUE, palette=NULL, pointSize=1,
        labelSize=NULL, axisSize=NULL, titleSize=NULL)

p1 <- plotRes(object=discoRes, comps=c(1, 2), what="loadings", type="individual",
              combined=FALSE, block="expr", color="classname", shape=NULL,
              labels=NULL, background=TRUE, palette=NULL, pointSize=1,
              labelSize=NULL, axisSize=NULL, titleSize=NULL)
p2 <- plotRes(object=discoRes, comps=c(1, 2), what="loadings", type="individual",
              combined=FALSE, block="micro", color="classname", shape=NULL,
              labels=NULL, background=TRUE, palette=NULL, pointSize=1,
              labelSize=NULL, axisSize=NULL, titleSize=NULL)
grid.arrange(arrangeGrob(p1+theme(legend.position="none"),
                         p2+theme(legend.position="none"), nrow=1), 
             heights=c(6/7, 1/7))

plotRes(object=discoRes, comps=c(1, 1), what="loadings", type="individual",
        combined=TRUE, block="", color="classname", shape=NULL,
        labels=NULL, background=TRUE, palette=NULL, pointSize=1,
        labelSize=NULL, axisSize=NULL, titleSize=NULL)
p1 <- plotRes(object=discoRes, comps=c(1, 1), what="loadings", type="both",
              combined=TRUE, block="expr", color="Time", shape=NULL,
              labels=NULL, background=TRUE, palette=NULL, pointSize=1,
              labelSize=NULL, axisSize=NULL, titleSize=NULL)
p2 <- plotRes(object=discoRes, comps=c(1, 1), what="loadings", type="both",
              combined=TRUE, block="micro", color="Time", shape=NULL,
              labels=NULL, background=TRUE, palette=NULL, pointSize=1,
              labelSize=NULL, axisSize=NULL, titleSize=NULL)
grid.arrange(arrangeGrob(p1+theme(legend.position="none"),
                         p2+theme(legend.position="none"), nrow=1), 
             heights=c(6/7, 1/7)) 

biplotRes(object=discoRes, type="common", comps=c(1, 2), block="",
          title=NULL, colorCol="HSCT_responder", sizeValues=c(2, 4),
          shapeValues=c(17, 0), background=TRUE, pointSize=1,
          labelSize=NULL, axisSize=NULL, titleSize=NULL) 
biplotRes(object=discoRes, type="common", comps=c(1, 2), block="",
          title=NULL, colorCol="Time", sizeValues=c(2, 4),
          shapeValues=c(17, 0), background=TRUE, pointSize=1,
          labelSize=NULL, axisSize=NULL, titleSize=NULL) 


# results <- omicsNPC(dataInput=sapply(A, function(x){Biobase::ExpressionSet(t(x))}, simplify = FALSE), dataTypes=c("count", "count"),
#                     combMethods="Fisher", numPerms=100,
#                     numCores=1, verbose=TRUE)
dev.off()
