#Task 1 Visualizing high dimensional data (PCA)


load("D:\\course\\Advanced bioinformatic tools\\E4\\exprData.Rdata")
ls()# list the contents of the object 
data.red #display the data
colnames(data.red)#display the name of each collunm in the dataset
nrow(data.red)#displays the number of rows
ncol(data.red)#displays the number of collunms 


df1<- t(data.red)#Transposition between the columns or rows
result<-prcomp(df1, scale=TRUE) #PCA is applied to the dataset with the variables being scaled.
summary(result)#Displays the information of each principal component
plot(result)#Displays the variance of each PC
str(result) #Check the structure of the result: the contents in x represent PC scores
result$x
df2<-cbind(df1,result$x[,1:2])
df2<-data.frame(df2)

install.packges("ggplot2")
library(ggplot2)


group<-c("low fat diet_day 3","low fat diet_day 3","low fat diet_day 3","low fat diet_day 3","low fat diet_day 3","high fat diet_day 3","high fat diet_day 3","high fat diet_day 3","high fat diet_day 3","high fat diet_day 3","high fat diet_day 3","low fat diet_day 28","low fat diet_day 28","low fat diet_day 28","low fat diet_day 28","low fat diet_day 28","low fat diet_day 28","high fat diet_day 28","high fat diet_day 28","high fat diet_day 28","high fat diet_day 28","high fat diet_day 28","high fat diet_day 28","low fat diet_day 3")      
# Create a vector named group, containing the information of each sample: low fat diet or high fat diet.
df3<-cbind(df2,group) #Add the information to the observation dataset.
ggplot(df3,aes(PC1,PC2, col=group))+geom_point()

library(factoextra)
fviz_contrib(result, choice="var", axes = 1, top = 5)#use function fviz_contrib to find the top 5 genes that contribute to PC1. 

#correlation between vars and PCs
cor(df3[,2],df3$PC1


#Task 2 Feature selection

library(limma)
data<-data.frame(data.red)#
group=c("low","low","low","low","low","high","high","high","high","high","high","low","low","low","low","low","low","high","high","high","high","high","high","low")

modelmatirx=model.matrix(~ 0 + group)
data.fit = lmFit(data,modelmatirx)
data.fit$coefficients[1:10,]

contrast.matrix = makeContrasts(high-low,levels=c("high","low"))
data.fit.con = contrasts.fit(data.fit,contrast.matrix)
data.fit.eb = eBayes(data.fit.con)
ls(data.fit.eb)
data.fit.eb$t[1:10,]
data.fit.eb$p.value[1:10,]
data.fit.eb$coefficients[1:1000,]

data.combined = cbind (data.fit.eb$p.value,data.fit.eb$coefficients)

topgenes = data.combined[data.combined[,1] < 0.001, ]
top = topgenes[(topgenes[, 2] > 1.9)|(topgenes[, 2] < -1.9), ] 
length(topgenes)
length(top)

top2 = topgenes[(topgenes[, 2] > 4.9)|(topgenes[, 2] < -4.9), ] 
length(top2)
data["Zkscan3",]
data$Zkscan3
data2<-t(data[c("Myl3","Myh7"),])

group=c("low","low","low","low","low","high","high","high","high","high","high","low","low","low","low","low","low","high","high","high","high","high","high","low")

data2<-data.frame(data2)
data3<-cbind(data2,group)
ggplot(data3,aes(Myl3,Myh7, col=group))+geom_point()

#Task 3 KNN classiﬁcation
library(class)#The knn.cv function from class package is based on the leave one out cross validation. 

library (class)
top5<- topgenes[(topgenes[, 2] > 4.65)|(topgenes[, 2] < -4.65), ]
top10 <- topgenes[(topgenes[, 2] > 3.5)|(topgenes[, 2] < -3.5), ]  #Previously, the top 2 and top 50 genes with the highest FCs were identified. Now, the top 5 and top 10 genes with the highest FCs are identified.
listtop2<-rownames(top2)
listtop5<-rownames(top5)
listtop10<-rownames(top10)
listtop50<-rownames(top50)  #Create name lists of the top genes in FC. 
data.t<-t(data.red)  #Transpose the data matrix.
group=c("low","low","low","low","low","high","high","high","high","high","high","low","low","low","low","low","low","high","high","high","high","high","high","low")  #Create a vector containing group name.
data.group<-data.frame(cbind(data.t,group))  #Create the data frame with samples classification information.
data.top2<- data.group[,c(listtop2,"group")]
data.top5<- data.group[,c(listtop5,"group")]
data.top10<- data.group[,c(listtop10,"group")]
data.top50<- data.group[,c(listtop50,"group")] # Data frame with different sets of top genes are created. The group names are also added. For unknown reason, the gene H2-Q10 can not be recognized in the data frame, data.group, a new top 50 gene list without gene H2-Q10 is created as follows.
top51 = topgenes[(topgenes[, 2] > 1.99)|(topgenes[, 2] < -1.99), ]
listtop51<-rownames(top51)  
newtop50 <-listtop51[-26] 
data.newtop50 <-data.group[,c(newtop50,"group")]


knn.top2.1<-knn.cv(train = data.top2 [,c(1:2)], cl =as.factor(data.top2[,3]), k = 1, prob = TRUE, use.all = TRUE)
knn.top2.3<-knn.cv(train = data.top2 [,c(1:2)], cl =as.factor(data.top2[,3]), k = 3, prob = TRUE, use.all = TRUE)
knn.top2.5<-knn.cv(train = data.top2 [,c(1:2)], cl =as.factor(data.top2[,3]), k = 5, prob = TRUE, use.all = TRUE)
knn.top5.1<-knn.cv(train = data.top5 [,c(1:5)], cl =as.factor(data.top5[,6]), k = 1, prob = TRUE, use.all = TRUE)
knn.top5.3<-knn.cv(train = data.top5[,c(1:5)], cl =as.factor(data.top5[,6]), k = 3, prob = TRUE, use.all = TRUE)
knn.top5.5<-knn.cv(train = data.top5[,c(1:5)], cl =as.factor(data.top5[,6]), k = 5, prob = TRUE, use.all = TRUE)
knn.top10.1<-knn.cv(train = data.top10[,c(1:10)], cl =as.factor(data.top10[,11]), k = 1, prob = TRUE, use.all = TRUE)
knn.top10.3<-knn.cv(train = data.top10[,c(1:10)], cl =as.factor(data.top10[,11]), k = 3, prob = TRUE, use.all = TRUE)
knn.top10.5<-knn.cv(train = data.top10[,c(1:10)], cl =as.factor(data.top10[,11]), k = 5, prob = TRUE, use.all = TRUE)
knn.top50.1<-knn.cv(train = data.newtop50[,c(1:50)], cl =as.factor(data.newtop50 [,51]), k = 1, prob = TRUE, use.all = TRUE)
knn.top50.3<-knn.cv(train = data.newtop50[,c(1:50)], cl =as.factor(data.newtop50 [,51]), k = 3, prob = TRUE, use.all = TRUE)
knn.top50.5<-knn.cv(train = data.newtop50[,c(1:50)], cl =as.factor(data.newtop50 [,51]), k = 5, prob = TRUE, use.all = TRUE)


accuracy<-function(factor, KNN){
table<-table(factor, KNN)
accuracy<-sum(diag(table))/sum(table)
return(accuracy)}
accuracy.top2.1<-accuracy(as.factor(data.top2[,3]), knn.top2.1)
accuracy.top2.3<-accuracy(as.factor(data.top2[,3]), knn.top2.3)
accuracy.top2.5<-accuracy(as.factor(data.top2[,3]), knn.top2.5)
accuracy.top5.1<-accuracy(as.factor(data.top5[,6]), knn.top5.1)
accuracy.top5.3<-accuracy(as.factor(data.top5[,6]), knn.top5.3)
accuracy.top5.5<-accuracy(as.factor(data.top5[,6]), knn.top5.5)
accuracy.top10.1<-accuracy(as.factor(data.top10[,11]), knn.top10.1)
accuracy.top10.3<-accuracy(as.factor(data.top10[,11]), knn.top10.3)
accuracy.top10.5<-accuracy(as.factor(data.top10[,11]), knn.top10.5)
accuracy.top50.1<-accuracy(as.factor(data.newtop50[,51]), knn.top50.1)
accuracy.top50.3<-accuracy(as.factor(data.newtop50[,51]), knn.top50.3)
accuracy.top50.5<-accuracy(as.factor(data.newtop50[,51]), knn.top50.5)

#Create a accuracy table
accuracy.table<-data.frame(rbind(accuracy.top2.1, accuracy.top2.3,accuracy.top2.5, accuracy.top5.1, accuracy.top5.3,accuracy.top5.5, accuracy.top10.1,accuracy.top10.3,accuracy.top10.5, accuracy.top50.1,accuracy.top50.3,accuracy.top50.5))
colnames(accuracy.table)="accuracy"
accuracy.table #The first 5 combinations have the best accuracy


ConfusionTable.top2.1<-table(as.factor(data.top2[,3]), knn.top2.1)
ConfusionTable.top2.3<-table(as.factor(data.top2[,3]), knn.top2.1)
ConfusionTable.top2.5<-table(as.factor(data.top2[,3]), knn.top2.1)
ConfusionTable.top5.1<-table(as.factor(data.top5[,6]), knn.top5.1)
ConfusionTable.top5.3<-table(as.factor(data.top5[,6]), knn.top5.3)
ConfusionTable.top5.5<-table(as.factor(data.top5[,6]), knn.top5.5)

#Task 4 Linear classiﬁcation

library(MASS)
data.top2
data.top5
data.top10
data.newtop50 

lda.top2<-lda(group~ .,data.top2, CV = TRUE)
lda.top5<-lda(group~ .,data.top5, CV = TRUE) 
lda.top10<-lda(group~ .,data.top10, CV = TRUE) 
lda.newtop50 <-lda(group~ .,data.newtop50 , CV = TRUE) 

accuracy.lda.top2<-accuracy(data.top2$group,lda.top2$class)
accuracy.lda.top5<-accuracy(data.top5$group,lda.top5$class)
accuracy.lda.top10<-accuracy(data.top10$group,lda.top10$class)
accuracy.lda.newtop50<-accuracy(data.newtop50$group,lda.newtop50$class)

accuracy.table<-data.frame(rbind(accuracy.lda.top2, accuracy.lda.top5,accuracy.lda.top10, accuracy.lda.newtop50))
colnames(accuracy.table)="accuracy"
accuracy.table 

table.top10<-table(data.top10$group,lda.top10$class)



#Task 5 Classiﬁcation based on days exposed to diet

df.lowfat<-t(data.red)[c(1:5,24,12:17),] #Extract the samples with low fat diet.
result.lowfat<-prcomp(df.lowfat, scale=TRUE)      #PCA is applied to the dataset with the variables being scaled.
result.lowfat$x
group.lowfat<-rep(c("day3","day28"), times=c(6,6)) #Create a vercot containing group information
df.lowfat.PC<-cbind(df.lowfat,result.lowfat$x[,1:2]) #Extract PC1 and PC2 from the PCA result, and combine it with the observation dataset.
df.lowfat.PC<-data.frame(df.lowfat.PC)
ggplot(df.lowfat.PC,aes(PC1,PC2, col=group.lowfat))+geom_point()      #Use ggplot to plot the graph with PC1 and PC2 as x and y axis. The factor for grouping is created previously (group.lowfat).

library(limma)

model.matrix<-model.matrix(~0+ group.lowfat)
data.fit<-lmFit(t(df.lowfat),model.matrix)#Transpose the previous data frame which contains the data of low fat diet sample. The mean of the expression level of each gene for each group is calculated based on previously created model matrix.  
contrast.matrix = makeContrasts(group.lowfatday28-group.lowfatday3,levels=c("group.lowfatday28","group.lowfatday3"))
data.fit.con = contrasts.fit(data.fit,contrast.matrix)
data.fit.eb = eBayes(data.fit.con)
data.combined = cbind (data.fit.eb$p.value,data.fit.eb$coefficients)
topgenes = data.combined[data.combined[,1] < 0.05, ]
top51 = topgenes[(topgenes[, 2] > 0.95)|(topgenes[, 2] < -0.95), ] 
length(top51)
top2 = topgenes[(topgenes[, 2] > 2.4)|(topgenes[, 2] < -2.4), ] 
length(top2)
 
top5<- topgenes[(topgenes[, 2] > 1.9)|(topgenes[, 2] < -1.9), ]
top11 <- topgenes[(topgenes[, 2] > 1.5)|(topgenes[, 2] < -1.5), ]  #Previously, the top 2 and top 51 genes with the highest FCs were identified. Now, the top 5 and top 10 genes with the highest FCs are identified.
listtop2<-rownames(top2)
listtop5<-rownames(top5)
listtop11<-rownames(top11)
listtop51<-rownames(top51)  #Create name lists of the top genes in FC. 
data.t<-t(data.red)[c(1:5,24,12:17),]  #Transpose the data matrix.
group.lowfat<-rep(c("day3","day28"), times=c(6,6)) #Create a vector containing group information
data.group<-data.frame(cbind(data.t,group.lowfat))  #Create the data frame with samples classification information.
data.top2<- data.group[,c(listtop2,"group.lowfat")]
data.top5<- data.group[,c(listtop5,"group.lowfat")]
data.top10<- data.group[,c(listtop11[-4],"group.lowfat")]
data.top49<- data.group[,c(listtop51[c(-19,-32)],"group.lowfat")] # Data frame with different sets of top genes are created. The group names are also added. For unknown reason, the gene H2-Q10 and 8430408G22Rik can not be recognized in the data frame, data.group, thus they are removed from the gene list.

# KNN classification with LOOCV is applied with the 12 parameter combinations 

knn.top2.1<-knn.cv(train = data.top2 [,c(1:2)], cl =as.factor(data.top2[,3]), k = 1, prob = TRUE, use.all = TRUE)
knn.top2.3<-knn.cv(train = data.top2 [,c(1:2)], cl =as.factor(data.top2[,3]), k = 3, prob = TRUE, use.all = TRUE)
knn.top2.5<-knn.cv(train = data.top2 [,c(1:2)], cl =as.factor(data.top2[,3]), k = 5, prob = TRUE, use.all = TRUE)
knn.top5.1<-knn.cv(train = data.top5 [,c(1:5)], cl =as.factor(data.top5[,6]), k = 1, prob = TRUE, use.all = TRUE)
knn.top5.3<-knn.cv(train = data.top5[,c(1:5)], cl =as.factor(data.top5[,6]), k = 3, prob = TRUE, use.all = TRUE)
knn.top5.5<-knn.cv(train = data.top5[,c(1:5)], cl =as.factor(data.top5[,6]), k = 5, prob = TRUE, use.all = TRUE)
knn.top10.1<-knn.cv(train = data.top10[,c(1:10)], cl =as.factor(data.top10[,11]), k = 1, prob = TRUE, use.all = TRUE)
knn.top10.3<-knn.cv(train = data.top10[,c(1:10)], cl =as.factor(data.top10[,11]), k = 3, prob = TRUE, use.all = TRUE)
knn.top10.5<-knn.cv(train = data.top10[,c(1:10)], cl =as.factor(data.top10[,11]), k = 5, prob = TRUE, use.all = TRUE)
knn.top49.1<-knn.cv(train = data.top49[,c(1:49)], cl =as.factor(data.top49[,50]), k = 1, prob = TRUE, use.all = TRUE)
knn.top49.3<-knn.cv(train = data.top49[,c(1:49)], cl =as.factor(data.top49[,50]), k = 3, prob = TRUE, use.all = TRUE)
knn.top49.5<-knn.cv(train = data.top49[,c(1:49)], cl =as.factor(data.top49[,50]), k = 5, prob = TRUE, use.all = TRUE)


accuracy<-function(factor, KNN){
table<-table(factor, KNN)
accuracy<-sum(diag(table))/sum(table)
return(accuracy)}  #A function is created to calculate the accuracy, the input parameters are KNN result and corresponding factor.
#Calculate the accuracy for the classification with different combinations
accuracy.top2.1<-accuracy(as.factor(data.top2[,3]), knn.top2.1)
accuracy.top2.3<-accuracy(as.factor(data.top2[,3]), knn.top2.3)
accuracy.top2.5<-accuracy(as.factor(data.top2[,3]), knn.top2.5)
accuracy.top5.1<-accuracy(as.factor(data.top5[,6]), knn.top5.1)
accuracy.top5.3<-accuracy(as.factor(data.top5[,6]), knn.top5.3)
accuracy.top5.5<-accuracy(as.factor(data.top5[,6]), knn.top5.5)
accuracy.top10.1<-accuracy(as.factor(data.top10[,11]), knn.top10.1)
accuracy.top10.3<-accuracy(as.factor(data.top10[,11]), knn.top10.3)
accuracy.top10.5<-accuracy(as.factor(data.top10[,11]), knn.top10.5)
accuracy.top49.1<-accuracy(as.factor(data.top49[,50]), knn.top49.1)
accuracy.top49.3<-accuracy(as.factor(data.top49[,50]), knn.top49.3)
accuracy.top49.5<-accuracy(as.factor(data.top49[,50]), knn.top49.5)
#Create an accuracy table
accuracy.table<-data.frame(rbind(accuracy.top2.1, accuracy.top2.3,accuracy.top2.5, accuracy.top5.1, accuracy.top5.3,accuracy.top5.5, accuracy.top10.1,accuracy.top10.3,accuracy.top10.5, accuracy.top49.1,accuracy.top49.3,accuracy.top49.5))
colnames(accuracy.table)="accuracy"
accuracy.table #The first 5 combinations have the best accuracy


#Create confusion table for the best combination. 
ConfusionTable.top49.3<-table(as.factor(data.top49[,50]), knn.top49.3)
ConfusionTable.top49.3

t<-rep(c("day3","day28"), times=c(6,6)) #Create a vercot containing group information
da<-data.frame(cbind(df.lowfat.PC,t))
ggplot(df.lowfat.PC,aes(Cidea,Hspa1a, col=t))+geom_point()  #A 2-D scatter plot is produced with the two 


