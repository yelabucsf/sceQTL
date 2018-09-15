library(lme4)
library(data.table)

cd52.expr=read.table('cd52.txt')
cd52.covs=fread('cd52.covs.txt', sep=',',header=T)

#lets put all of our variables in here
df=data.frame(cd52.covs, cd52=cd52.expr)

#just get the monocytes for now, and just SLE
df.use=df[intersect(grep('CD14', df$ct_cov), which(df$disease_cov=='sle')), ]


genos=read.table('cd52.geno.txt')
inds=sapply(strsplit(colnames(genos), 'X'), '[',2)
genos=as.numeric(genos)

df.use$geno=genos[match(df.use$ind_cov, inds)]

#lm is working.. thats a good sign! but it probably works tooooo well because we're not accounting for effects w/in an individual

print('Linear Model')
print(coef(summary(lm(V1 ~ geno , data=df.use))))

print('Linear Mixed Model')
print(coef(summary(lmer(V1 ~ geno + (1|ind_cov), data=df.use))))


#try a number of group sizes for pseudobulk to see how t value changes this needs the raw counts.. maybe try this later

#how many cells do you need??
cells=seq(1000, nrow(df.use), by=1000)
tvals=c()
for(c in cells){
  print(c)
  t=coef(summary(lmer(V1 ~ geno + (1|ind_cov), data=df.use[sample( nrow(df), c), ])))[2, 3]
  tvals=c(tvals, t)
}
pdf('tvals.cells.pdf')
plot(cells, tvals)
abline(h=16.5754120606288, col='red')
dev.off()

#sample cells to get an estimate of the std error of the variance  or t value estimate)

#if you seq'd 40,000 cells

n.iter=100
tvals.50k=c()
for(i in 1:n.iter){
  if(i %%10==0){print(i)}
  t=coef(summary(lmer(V1 ~ geno + (1|ind_cov), data=df.use[sample(nrow(df), 50000),])))[2, 3]
  tvals.50k=c(tvals.50k, t)

}


n.iter=100
tvals.80k=c()
for(i in 1:n.iter){
  if(i %%10==0){print(i)}
  t=coef(summary(lmer(V1 ~ geno + (1|ind_cov), df.use[sample(nrow(df), 80000),])))[2, 3]
  tvals.80k=c(tvals.80k, t)

}

df=data.frame(vals=c(tvals.50k, tvals.80k), cells=c(rep('50k', 100), rep('80k', 100)))

pdf('dist.cells.pdf')
ggplot(df, aes(vals, fill=cells))+ geom_density(alpha=0.5)
dev.off()
