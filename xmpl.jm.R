# Example for HEV seminar

# load data
load(file = "data/xmpl.jm.RData")

library(JM)
# survival part
surv.jm <- coxph(Surv(time.event, event) ~ sex + dx + I(log2(uprot+0.001)) + I(age/10), 
                 data = d.jm.surv, x= TRUE )
summary(surv.jm)

# the comparator survival model
surv.only <- coxph(Surv(time.event, event) ~ sex + dx + I(log2(uprot+0.001)) + I(age/10) + egfr, 
                 data = d.jm.surv, x= TRUE )
summary(surv.only)
ci.surv.only <- cbind(surv.only$coefficients, confint(surv.only))
ci.surv.only

# LMM part, prediction only work for unstructured (default) or Diagonal random effects.
lme.jm <- lme(fixed = egfr ~ ns(time, df = 3) + sex + I(age/10) + dx + log2(uprot+0.001),
              random = list(pid = pdDiag(form = ~ ns(time,df = 3))),
              control = lmeControl(maxIter=50, msMaxIter=50, msVerbose = FALSE, returnObject = TRUE), 
              data = d.jm 
)
summary(lme.jm)

# plot the predicted trajectory per patient
library(ggplot2)
# select first 50 patients
list.patients <- split(d.jm, d.jm$pid)
list.patients <- list.patients[c(1:50)]
d.plot <- do.call(rbind, list.patients)
d.plot <- d.plot[order(d.plot$pid, d.plot$time),]
#
d.jm$fitted.gfr <- predict(lme.jm, data = d.plot)
plot <- ggplot(data = d.plot) +
  geom_line(aes(x = time, y = fitted.gfr, group = pid), size = 0.1)
plot + coord_cartesian(ylim = c(0,150)) + scale_y_continuous(breaks = seq(from = 0, to = 150, by = 15)) +
  coord_cartesian(xlim = c(0,6)) + scale_x_continuous(breaks = seq(from = 0,to = 6,by = 1)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", fill = "white"))
 
# join using both time-updated eGFR value
jm <- jointModel(lme.jm, surv.jm, 
                      timeVar = "time",
                      param = "value"
)
summary(jm)

# Make sure that both files are in same order
d.jm <- d.jm[order(d.jm$pid,d.jm$time),]
d.jm.surv <- d.jm.surv[order(d.jm.surv$pid),]

jm.bayes <- jointModelBayes(lme.jm, surv.jm, 
                 timeVar = "time",
                 param = "td-value",
                 control = list(n.iter = 1000, n.burnin = 150, n.thin = 100, n.adapt = 150) 
)
summary(jm.bayes) 

ci <- round(confint(jm),3)
write.table(ci, file = "masterplan_jm.csv", quote=FALSE, sep =";", col.names = TRUE, row.names = TRUE)

# Dynamic C index
dynCJM(jm, newdata = d.jm, idVar = "pid", Dt=1, t.max = 5)

# Calibration

# survfitJM predicts survival up to a time s+t, given survival until s
# the last.time argument provides the landmark time s
# by default survival up to a sequence of 35 times is evaluated, by providing the argument survTimes a specific time t
# is requested instead
list.patients <- split(d.jm, d.jm$pid)
pred.jm <- function(x) { 
  pred <- survfitJM(jm,
                    newdata = list.patients[[x]], 
                    idVar = "pid", 
                    simulate = FALSE, # set TRUE if confidence interval are of interest
                    last.time = 1,
                    survTimes = 5
  )
  return(pred)
}
n <- c(1:length(list.patients))
pred <- lapply(n, pred.jm) 
# Extract and join the results
extract <- function(x) {
  ex <- data.frame(pred[[x]]$summaries[1])
}
d.pred <- lapply(n, extract)
for (i in 1:length(list.patients)) {
  names(d.pred[[i]]) <- c("horizon", "surv")
  d.pred[[i]]$pid <- names(list.patients[i])
  d.pred[[i]]$last.time <- pred[[i]]$last.time
}
d.pred.long<- do.call(rbind, d.pred)

# now get the observed survival and times from the landmark onward
obs <- data.frame(d.jm.surv$pid)
obs$time.event <- d.jm.surv$time.event
obs$event <- d.jm.surv$event
names(obs) <- c("pid", "time.event", "event") 
# Join observed to predicted outcome
join <- merge(obs, d.pred.long, by = "pid")
join$surv <- 1 - join$surv
# Keep only patients with survival up to the landmark
join$time.event <- join$time.event - join$last.time
join <- join[join$time.event>0,]

# store calibration data
names(join) <-  c("pid", "time.event", "event", "horizon", "risk", "landmark")
write.table(join, 
            file="calibration_jm.csv", 
            col.names=TRUE, row.names=FALSE, quote=FALSE, sep=",")

# clean up
rm(list = ci, d.plot, d.pred.long, join, obs, d.pred, i, list.patients, myvars, n, plot, pred, sex)

# get subjects
for (i in c(1002,1004,1031,1024,1020,1041)) {
  write.table(d.jm[d.jm$pid == i,],
              file = paste("subject",".",i,".csv"),
              col.names=TRUE, row.names=FALSE, quote=FALSE, sep=",")
}

# Run the Shiny app
library(JMbayes)
runDynPred()

