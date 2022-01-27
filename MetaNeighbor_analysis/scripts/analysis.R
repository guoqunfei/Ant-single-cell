library(ggplot2)
library(pvclust)
library(reshape2)
library(RColorBrewer)
		
data<-read.table("pseudo.bulk.roboust.exp.txt")
data<-log(data*10+1)
result <- pvclust(data, method.dist="cor", method.hclust="ward.D", nboot=1000, parallel=TRUE)
pdf("Amel.KC.marker.pvclust.pdf")
plot(result)
dev.off()
data<-t(scale(t(data),scale=F, center=T))
data<-as.matrix(data)
cor.matrix<-cor(data)
cor.data<-melt(cor.matrix)

plot_theme<-theme(panel.background = element_blank(),axis.line = element_line(size=0.1),axis.ticks = element_line(size=0.1),axis.text = element_text(size=3),axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),axis.ticks.length = unit(1, "pt"),plot.title = element_text(size = 3, face = "bold",margin=margin(0,0,4,0)),plot.margin = unit(c(0,0,0,0),"cm"),legend.key.size = unit(0.2, 'cm'),legend.key.height = unit(0.2, 'cm'),legend.key.width = unit(0.2, 'cm'),legend.title = element_text(size=4),legend.text = element_text(size=4))

cor.data$Var1<-factor(cor.data$Var1, levels=rev(c("Amel_c5","Amel_c8","Amel_c4","Amel_c6","Amel_c2","Amel_c0")))
cor.data$Var2<-factor(cor.data$Var2, levels=c("Amel_c5","Amel_c8","Amel_c4","Amel_c6","Amel_c2","Amel_c0"))
#cor.data$value<-ifelse(cor.data$value< -0.4, -0.4, cor.data$value)
write.table(cor.data, "cor.txt",sep="\t",quote=F, row.names=F)
#cor.data$value<-ifelse(cor.data$value>0.2,0.2,cor.data$value)

pdf("Amel.KC.cor.heatmap.2.pdf",height=4/2.54,width=4/2.54)
#ggplot(cor.data, aes(Var1,Var2, fill = value))+ geom_tile(color = "white")+scale_fill_gradientn(colours=c("#0000FF","#FFFFFF","#FFCCCC","#FF0000"))+coord_fixed()+labs(x="",y="")+plot_theme
ggplot(cor.data, aes(Var1,Var2, fill = value))+ geom_tile(color = "white")+scale_fill_gradient2(low="#0000FF",mid="#FFFFFF",high="#FF0000",midpoint=0,breaks=c(-0.8,-0.4,0,0.4,0.8))+coord_fixed()+labs(x="",y="")+plot_theme
dev.off()

#kcdata<-subset(data,seurat_clusters %in% c(0,2,4,5,6,8))
#kcdata@active.ident <- factor(kcdata@active.ident,levels=c(2,0,6,4,8,5))
#DotPlot(kcdata,"RNA",c("Mblk-1","Camkii","Tk","LOC725542","E74"),cols = c("lightgrey", "red"),col.min= -0.5,col.max=1.5)+labs(x="",y="")+scale_x_discrete(labels=c("Mblk-1","CaMKII","Tk","mKast","E74"))+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

#DotPlot(kcdata,"RNA",c("Mblk-1","Camkii","LOC725542","Tk","E74"),cols = c("lightgrey", "red"),col.min= 0)+labs(x="",y="")+scale_x_discrete(labels=c("Mblk-1","CaMKII","mKast","Tk","E74","Ecr","HR38"))+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+plot_theme
