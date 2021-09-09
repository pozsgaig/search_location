#### Loading necessary packages and functions ####

library(data.table) # for character issues and calculating AADP
library(scales) # for rescaling data
library(olsrr) # for testing heteroscedasticity
library(fitdistrplus) # for checking dsitribution
library(car) # for boxcox transformation
library(yarrr) # for pirate plots
library(BiodiversityR) # quick plotting and vegdist
library(reshape2) # for melt
library(ggplot2) # plotting
library(ggpubr) # plotting
library(sp) #for convex hull areas
library(welchADF) # for Welch ADF


### function replacing NAs with 0
NA20<- function(x)
{
        x[is.na(x)]<-0
        x}


#### Data loading and data preparation ####
search_hits<-read.table("ecol_search_hits.csv", header = T, quote = "\"")

# Centering data
search_hits_group_means<- data.table(search_hits)


search_hits_group_means[, AAD:=abs(Number_hits-mean(Number_hits)),  
                        by = c("Keyword_complexity", "Search_engine", "Browser", "Cache")]


search_hits_group_means[, meanHits:=mean(Number_hits),  
                        by = c("Keyword_complexity", "Search_engine", "Browser", "Cache")]

search_hits_group_means[, sumHits:=sum(Number_hits),  
                        by = c("Keyword_complexity", "Search_engine", "Browser", "Cache")]




search_hits_group_means[, AADP := abs(Number_hits-mean(Number_hits))/mean(Number_hits), 
                        by = c("Keyword_complexity","Search_engine")]

# non-absolute average deviation
search_hits_group_means[, NAADP := (Number_hits-mean(Number_hits))/mean(Number_hits), 
                        by = c("Keyword_complexity","Search_engine")]


search_hits_group_means[, scaledHits_SD := sd(Number_hits)/mean(Number_hits), 
                        by = c("Keyword_complexity","Search_engine")]

search_hits_group_means[, logscaledHits := log(1+abs(Number_hits-mean(Number_hits))/mean(Number_hits)), 
                        by = c("Keyword_complexity","Search_engine")]


search_hits_group_means$loghits<-log(1+search_hits_group_means$Number_hits)

# The percentage one institution got from the max value of hits per search expressions


search_hits_group_means[, gmedianhits := max(Number_hits), 
                        by = c("Search_engine","Keyword_complexity", "Topic")]

search_hits_group_means$maxperc<-search_hits_group_means$Number_hits/search_hits_group_means$gmedianhits

search_hits_group_means$negAADP<-1-search_hits_group_means$AADP

# add scaled (0,1) AADP (lehetne (0.00000001, 0.99999999) is betareg-hez)
search_hits_group_means$SAADP<-rescale(search_hits_group_means$AADP, to=c(0.000001, 0.999999))

# add grouping factor
search_hits_group_means$gr_fact<-as.factor(paste0(search_hits_group_means$Search_engine,
                                                  search_hits_group_means$Keyword_complexity, 
                                                  search_hits_group_means$Browser, 
                                                  search_hits_group_means$Cache))



### Preparation for multivariate analysis

ecol_lines<-read.table("ecol_search_lines.csv", header = T, quote = "\"")

ecol_searches<-read.table("ecol_search_twenty.csv", header =  T, quote = "\"")

com_mat<-with(ecol_searches,
              tapply(First_auth, list(paste("SL", Search_line, Affiliation, Computer, Replica, sep="_"), num_id),  length)
)



com_mat<-NA20(com_mat)
com_mat<-as.data.frame(com_mat)

sort(rowSums(com_mat))
sort(colSums(com_mat))


env_mat<-with(ecol_searches,
              aggregate(Author, list(Search_line=Search_line, Affiliation = Affiliation, Computer=Computer, Replica=Replica), length)
)
env_mat[,5]<-NULL

nrow(env_mat)
nrow(com_mat)


env_mat<-merge(env_mat, ecol_lines, by= "Search_line")
rownames(env_mat)<-paste("SL", env_mat$Search_line, env_mat$Affiliation, env_mat$Computer, env_mat$Replica, sep="_")

com_mat<-com_mat[rownames(env_mat),]

rownames(env_mat)
rownames(com_mat)


####################################

#### Testing normality and homoscedasticity ####

### Checking distributions with fitdistrplus package
### https://stats.stackexchange.com/questions/132652/how-to-determine-which-distribution-fits-my-data-best
par(ask=T)
with(search_hits_group_means, 
     by(search_hits_group_means, gr_fact, function(x) descdist(x$AADP, discrete = F))
)
descdist(search_hits_group_means$AADP, discrete = F)

### Testing heteroscedasticity with olsrr package (Breusch Pagan Test)
### https://cran.r-project.org/web/packages/olsrr/vignettes/heteroskedasticity.html
lm1<-lm(AADP ~ Search_engine*(Keyword_complexity+Browser+Cache), data=search_hits_group_means)

ols_test_breusch_pagan(lm1, rhs = TRUE)

with(search_hits_group_means, 
     by(search_hits_group_means, Keyword_complexity, function(x) {
             lm1<-lm(log(0.11+AADP) ~ Search_engine*Browser*Cache, data=x)
             m<-ols_test_breusch_pagan(lm1, rhs = TRUE)
             return(m)}
     )
)

### Trying log transformation
search_hits_group_means[, logscaledHits := log(1+abs(Number_hits-mean(Number_hits))/mean(Number_hits)), 
                        by = c("Keyword_complexity","Search_engine", "Topic", "Browser", "Cache")]

### Arcsin transformation does not work because of the values greater than 1
### Trying square root transformation
search_hits_group_means[, sqrtscaledHits := sqrt(AADP), 
                        by = c("Keyword_complexity","Search_engine", "Topic", "Browser", "Cache")]

### Trying boxcox transformations
# with car package
search_hits_group_means[, boxcoxgHits := boxCoxVariable(AADP+0.0000001),  
                        by = c("Keyword_complexity","Search_engine", "Topic", "Browser", "Cache")]



###################################

#### Welch tests ####


omnibus_LSM <- welchADF.test(search_hits_group_means, response = "negAADP", between.s =
                                     c("Search_engine", "Keyword_complexity", "Browser", "Cache"), contrast = "omnibus")


omnibus_trimmed <- update(omnibus_LSM, trimming = TRUE)

omnibus_trimmed_boot <- update(omnibus_trimmed, bootstrap = TRUE, seed = 12345)

summary(omnibus_LSM)

summary(omnibus_trimmed_boot)


pairwise_trimmed <- welchADF.test(AADP ~ Search_engine, data = temp, effect = "Search_engine",
                                  contrast = "all.pairwise", trimming = TRUE, effect.size = TRUE)

pairwise_trimmed_boot <- update(pairwise_trimmed, bootstrap = TRUE, seed = 12345)
summary(pairwise_trimmed)





##################################

#### Basic stats and graphs ####
# Hit number means
tapply(search_hits$Number_hits, list(search_hits$Keyword_complexity,
                                       search_hits$Search_engine), function(x) sd(x)/mean(x))


tapply(search_hits_group_means$Number_hits, 
       list(search_hits_group_means$Keyword_complexity,
            search_hits_group_means$Search_engine), sd)


### Tables for paper 

grouped_hits_mean<-with(search_hits,
                tapply(Number_hits, list(paste(Search_engine, Browser, Cache, sep="_"), Keyword_complexity), mean))

grouped_hits_sd<-with(search_hits,
              tapply(Number_hits, list(paste(Search_engine, Browser, Cache, sep="_"), Keyword_complexity), sd))

# AADP
tapply(search_hits_group_means_$AADP, 
       list(search_hits_group_means$Keyword_complexity,
            search_hits_group_means$Search_engine), mean)
tapply(search_hits_group_means$AADP, 
       list(search_hits_group_means$Keyword_complexity,
            search_hits_group_means$Search_engine), sd)



### Plotting

par(mfrow=c(1,1), las=1)
pirateplot(logscaledHits~Search_engine, data=search_hits_group_means)

par(mfrow=c(1,1), las=1)
pirateplot(AADP~Search_engine, data=search_hits_group_means)

par(mfrow=c(1,1), las=1)
pirateplot(NAADP~Search_engine, data=search_hits_group_means)

par(mfrow=c(2,2), las=1)

for (i in unique(search_hits_group_means$Topic))
{
        for (j in unique(search_hits_group_means$Keyword_complexity))
        {
                x<- search_hits_group_means[search_hits_group_means$Topic==i & search_hits_group_means$Keyword_complexity==j,]
                pirateplot(log(0.001+AADP) ~ Search_engine, data=x, main= paste(i, j))  
        }
        
        
}

### GGPLOT version

ggboxplot(data=search_hits_group_means, x="Search_engine", y="AADP", color="Search_engine",
          palette = "jco", add = "jitter", add.params = list(size = 1, alpha = 0.2), lwd=0.7,
          facet.by = "Keyword_complexity", 
          xlab = "Search engine", ylab = "AADP", legend.title="") +
        rotate_x_text()


par(mfrow=c(2,2), las=1)
for (i in unique(search_hits_group_means$Topic))
{
        for (j in unique(search_hits_group_means$Keyword_complexity))
        {
                x<- search_hits_group_means[search_hits_group_means$Topic==i & search_hits_group_means$Keyword_complexity==j,]
                pirateplot(AADP ~ Search_engine*Browser, data=x, main= paste(i, j))  
        }
        
        
}



pirateplot(AADP ~ Search_engine*Cache, data=search_hits_group_means)





par(mfrow=c(2,1), las=2)
for (n in unique(env_mat$Keyword_complexity))
{
        sim_data<-NULL
        for (k in unique(env_mat$Search_engine))
        {
                e_mat<-env_mat[as.character(env_mat$Keyword_complexity) == n & as.character(env_mat$Search_engine) == k,]
                c_mat<-com_mat[rownames(com_mat) %in% rownames(e_mat),]
                c_mat<-c_mat[rowSums(c_mat)>0,colSums(c_mat)>0]
                
                e_mat<-e_mat[rownames(e_mat) %in% rownames(c_mat),]
                c_mat<-c_mat[order(rownames(c_mat)),]
                e_mat<-e_mat[order(rownames(e_mat)),]
                
                c_mat <-as.matrix(vegdist(c_mat, "jaccard", binary=F))
                
                
                sim_df<-subset(melt(c_mat))
                sim_df$Search_engine<-rep(k, nrow(sim_df))
                sim_data<-rbind(sim_data, sim_df)
        }
        colnames(sim_data)[3]<-"Similarity"
        plottitle<-paste(substring(as.character(e_mat[1,"Keyword_complexity"]), 
                                   1, nchar(as.character(e_mat[1,"Keyword_complexity"]))-4))
        pirateplot(Similarity~Search_engine, 
                   data=sim_data, main = plottitle,
                   xlab = "")
}




#################################

#### Multivariate for unique papers ####

### Convex hulls


article_dist<-vegdist(com_mat, distance='jaccard')
Ordination.model1 <- monoMDS(article_dist, 
                             data=env_mat)

Ordination.model1 <- capscale(com_mat ~ Search_engine*(Browser +Cache), 
                              data=env_mat, distance='jaccard', sqrt.dist=F, add=F)

plot1 <- ordiplot(Ordination.model1, choices=c(1,2), main=n)
plot(Ordination.model1)
check.ordiscores(com_mat, Ordination.model1, check.species=T)
summary(Ordination.model1, scaling='sites')
eigenvals(Ordination.model1)
RsquareAdj(Ordination.model1)
deviance(Ordination.model1)
vif.cca(Ordination.model1)

permutest(Ordination.model1, permutations=100)
anova.cca(Ordination.model1, step=100, by='terms')
anova.cca(Ordination.model1, step=100, by='margin')

par(mfrow=c(1,1), cex=1)
plot1 <- ordiplot(Ordination.model1, type='none',choices=c(1,2), 
                  scaling='sites')
attach(env_mat)
summary(ordiellipse(plot1, groups=Search_engine, conf=0.9, kind='sd', 
                    show.groups = env_mat[env_mat$Keyword_complexity=="simple_KWE",]))

par (mfrow=c(2,1))
all_hulls<-NULL
for (n in unique(env_mat$Keyword_complexity))
{       
        
        cat(n)
        cat("\n")
        
        #for testing
        # n="medicine_complex_KWE"
        e_mat<-env_mat[env_mat$Keyword_complexity == n,]
        
        nrow(c_mat)
        c_mat<-com_mat[rownames(com_mat) %in% rownames(e_mat),]
        # c_mat<-c_mat[rowSums(c_mat)>0,colSums(c_mat)>0]
        
        e_mat<-e_mat[rownames(e_mat) %in% rownames(c_mat),]
        c_mat<-c_mat[order(rownames(c_mat)),]
        e_mat<-e_mat[order(rownames(e_mat)),]
        c_mat <- as.matrix(vegdist(c_mat, "jaccard", binary=T))
        
        melted_c_mat <- melt(c_mat)
        ggplot(data = melted_c_mat, aes(x=Var1, y=Var2, fill=value)) + 
                geom_tile()
        Ordination.model1 <- capscale(c_mat ~ Search_engine*(Browser +Cache), data=e_mat)
        Ordination.model1$points<-rescale(Ordination.model1$points, to=c(1, 100))
        plot1 <- ordiplot(Ordination.model1, type='none',choices=c(1,2), main=n)
        
        env_hull<-with(e_mat, ordihull(plot1, groups=Search_engine), col=1:4)
        # env_hull<-lapply(env_hull, FUN= function(x) rescale(x, to=c(1, 100)))
        
      
        adonis_result<-with(e_mat,
                            adonis2(c_mat ~ Search_engine*Browser*Cache, data=e_mat, method='jaccard',
                                   permutations=100)
        )  
        
        hullsizes<-sapply(names(env_hull), function(x) {chull.poly <- Polygon(env_hull[x], hole=F)
        chull.area <- chull.poly@area
        return(chull.area)})
        
        hullsizes<-as.data.frame(t(hullsizes))
        rownames(hullsizes)<-n
        all_hulls<-rbind(all_hulls, hullsizes)
        
        print(hullsizes)
        cat("\n")
        print(adonis_result)
        cat("\n")
        cat("******************")
        
}


