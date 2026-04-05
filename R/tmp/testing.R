library(MLMusingR)
data(sch29)
head(sch29)
g1 <- lm(math ~ male + mses + mhmwk, data = sch29)
re = wcb2(math ~ male + mses + mhmwk, cluster = 'schid', dat = sch29,
       B = 9999
        )
re$res
aa = wildcb(math ~ male + mses + mhmwk, cluster = 'schid', dat = sch29,
            B = 9999, param = 'male')
aa$res
aa$setype
aa$wctype

library(fwildclusterboot)
bb = boottest(g1, clustid = 'schid', param = 'male', B = 9999,
              impose_null = T)
bb
plot(bb)

####

library(microbenchmark)
microbenchmark(
  a = wcb2(math ~ male + mses + mhmwk, cluster = 'schid', dat = sch29,
           B = 9999
  ),
  b = wildcb(math ~ male + mses + mhmwk, cluster = 'schid', dat = sch29,
             B = 9999),
  c = boottest(g1, clustid = 'schid', param = 'male', B = 9999,
               impose_null = F),
  times = 100
)

library(clusterSEs)
g2 <- lm(math ~ male + mses + mhmwk, data = sch29)

# cluster.bs.glm(g2, dat = sch29, cluster = ~schid)
cluster.wild.glm(g2, dat = sch29, cluster = ~schid)
gg = wildcb(math ~ male + mses + mhmwk, cluster = 'schid',
       dat = sch29, param = 'mhmwk')
gg$res

vcovBS(g2, type = 'wild', cluster = sch29$schid)
res3 <- ClusterBootstrap::clusbootglm(g2, data = sch29, clusterid = schid)
summary(res3)

citation(
  "clusterSEs"
)

