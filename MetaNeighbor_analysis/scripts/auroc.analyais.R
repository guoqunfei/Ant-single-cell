require(ggplot2)
require(RColorBrewer)
require(reshape2)
require(cowplot)

#source("/ldfssz1/ST_DIVERSITY/P18Z10200N0107/zhangpei/Ant/2.Project1/3.expre/single_cell/2.Cell/MetaNeighbor-master/2017-08-28-runMN-US.R")
#source("/ldfssz1/ST_DIVERSITY/P18Z10200N0107/zhangpei/Ant/2.Project1/3.expre/single_cell/2.Cell/MetaNeighbor-master/2017-08-28-runMN-US.pearson.R")
#data<-read.table("merged.cell.exp.marker.z.txt")
#pheno<-read.table("merged.pheno.cell.txt",header=T)
#celltypes<-read.table("merged.celltypes.list",heade=F)
#celltypes<-as.character(celltypes[,1])
#var.genes<-read.table("marker.list",header=F)
#var.genes<-as.character(var.genes[,1])
#AUROC.matrix.s=run_MetaNeighbor_s(var.genes, data, celltypes, pheno)
#AUROC.matrix.p=run_MetaNeighbor_p(var.genes, data, celltypes, pheno)
#AUROC.data.s<-melt(AUROC.matrix.s,value.name = "AUROCs")
#AUROC.data.p<-melt(AUROC.matrix.p,value.name = "AUROCp")
#res<-merge(AUROC.data.p,AUROC.data.s,by=c("Var1","Var2"))
#write.table(res,"AUROC.z.txt",sep="\t",quote=F,row.names=F)


data<-read.table("AUROC.z.txt",header=T)

plot_theme<-theme(panel.background = element_blank(),axis.line = element_line(size=0.1),axis.ticks = element_line(size=0.1),axis.text = element_text(size=1.5),axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),axis.ticks.length = unit(1, "pt"),plot.title = element_text(size = 3, face = "bold",margin=margin(0,0,4,0)),plot.margin = unit(c(0,0,0,0),"cm"),legend.key.size = unit(0.2, 'cm'),legend.key.height = unit(0.2, 'cm'),legend.key.width = unit(0.2, 'cm'),legend.title = element_text(size=3),legend.text = element_text(size=3))


MphaDmel<-subset(data, grepl("Mpha",Var1) & grepl("Dmel",Var2))
MphaDmel$AUROCp2<-ifelse(MphaDmel$AUROCp<=0.5, "", round(MphaDmel$AUROCp,2))
Dmel_order<-read.table("Dmel.order.list",header=F)

MphaDmel$Var1<-factor(MphaDmel$Var1,levels=c("Mpha_c15","Mpha_c16","Mpha_c27","Mpha_c6","Mpha_c25","Mpha_c20"))
MphaDmel$Var2<-factor(MphaDmel$Var2,levels=c(Dmel_order[,1]))


p1<-ggplot(MphaDmel, aes(Var2,Var1, fill = AUROCp))+ geom_tile(color = "white")+scale_fill_distiller(palette = "RdBu", limits = c(0,1))+geom_text(aes(label=AUROCp2),size=0.5)+coord_fixed()+labs(x="",y="")+plot_theme
p2<-ggplot(MphaDmel, aes(Var2,Var1, fill = AUROCs))+ geom_tile(color = "white")+scale_fill_distiller(palette = "RdBu", limits = c(0,1))+coord_fixed()+labs(x="",y="")+plot_theme
ggsave(filename="OL.MphaDmelNeuron.z.auroc.pdf", plot=plot_grid(p1,p2,ncol = 1),height=6/2.54,width=8/2.54)

