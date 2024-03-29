## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  dev = "png",
  collapse = TRUE,
  message =FALSE,
  warning =FALSE,
  fig.width = 7,
  comment = "#>"
)
if (capabilities(("cairo"))) {
  knitr::opts_chunk$set(dev.args = list(png = list(type = "cairo")))
}
options(rmarkdown.html_vignette.check_title = FALSE)
library(coveffectsplot)
library(mrgsolve)
library(ggplot2)
library(ggstance)
library(ggridges)
library(tidyr)
library(dplyr)
library(table1)
library(egg)
library(data.table)
library(patchwork)
library(GGally)
theme_set(theme_bw())
#utility function to simulate varying one covariate at a time keeping the rest at the reference
expand.modelframe <- function(..., rv, covcol="covname") {
  args <- list(...)
  df <- lapply(args, function(x) x[[1]])
  df[names(rv)] <- rv
  res <- lapply(seq_along(rv), function(i) {
    df[[covcol]] <- names(rv)[i]
    df[[names(rv)[i]]] <- args[[names(rv)[i]]]
    as.data.frame(df)
  })
  do.call(rbind, res)
}
round_pad <- function(x, digits = 2, round5up = TRUE) {
  eps <- if (round5up) x * (10^(-(digits + 3))) else 0
  formatC(round(x + eps, digits), digits = digits, format = "f", flag = "0")
}


derive.exposure <- function(time, CP) {
  n <- length(time)
  x <- c(
    Cmax = max(CP),
    Clast = CP[n],
    AUC = sum(diff(time) * (CP[-1] + CP[-n])) / 2
  )
  data.table(paramname=names(x), paramvalue=x)
}

cor2cov <- function (cor, sd) 
{
  if (missing(sd)) {
    sd <- diag(cor)
  }
  diag(cor) <- 1
  n <- nrow(cor)
  diag(sd, n) %*% cor %*% diag(sd, n)
}

stat_sum_df <- function(fun, geom="point", ...) {
  stat_summary(fun.data = fun,  geom=geom,  ...)
}

summary_conc <- function(CP) {
  x <- c(
    Cmed = median(CP),
    Clow = quantile(CP, probs = 0.05),
    Cup =  quantile(CP, probs = 0.95)
  )
  data.table(paramname=names(x), paramvalue=x)
}


nbsvsubjects <- 1000
nsim <- 10 # uncertainty replicates for vignette you might want a higher number


## ----pkmodel, collapse=TRUE---------------------------------------------------
codepkmodelcov <- '
$PARAM @annotated
KA    : 0.5   : Absorption rate constant Ka (1/h)
CL    : 4     : Clearance CL (L/h)
V     : 10    : Central volume Vc (L)
Vp    : 50    : Peripheral volume Vp (L)
Qp    : 10    : Intercompartmental clearance Q (L/h)
CLALB : -0.8  : Ablumin on CL (ref. 45 g/L)
CLSEX : 0.2   : Sex on CL (ref. Female)
CLWT  : 1     : Weight on CL (ref. 85 kg)
VSEX  : 0.07  : Sex on Vc (ref. Female)
VWT   : 1     : Weight on Vc (ref. 85 kg)

$PARAM @annotated // reference values for covariate
WT    : 85    : Weight (kg)
SEX   : 0     : Sex (0=Female, 1=Male)
ALB   : 45    : Albumin (g/L)

$PKMODEL cmt="GUT CENT PER", depot=TRUE, trans=11

$MAIN
double CLi = CL *
    pow((ALB/45.0), CLALB)*
    (SEX == 1.0 ? (1.0+CLSEX) : 1.0)*
    pow((WT/85.0), CLWT)*exp(nCL); 
double V2i = V *
    (SEX == 1.0 ? (1.0+VSEX) : 1.0)*
    pow((WT/85.0), VWT)*exp(nVC);  

double KAi = KA;
double V3i = Vp *pow((WT/85.0),    1);
double Qi = Qp *pow((WT/85.0), 0.75);

$OMEGA @annotated @block
nCL : 0.09       : ETA on CL
nVC : 0.01 0.09  : ETA on Vc

$TABLE
double CP   = CENT/V2i;

$CAPTURE CP KAi CLi V2i V3i Qi WT SEX ALB
'
modcovsim <- mcode("codepkmodelcov", codepkmodelcov)

partab <- setDT(modcovsim@annot$data)[block=="PARAM", .(name, descr, unit)]
partab <- merge(partab, melt(setDT(modcovsim@param@data), meas=patterns("*"), var="name"))
knitr::kable(partab)

## ----pksimulation, fig.width=7, message=FALSE---------------------------------
idata <- data.table(ID=1:nbsvsubjects, WT=85, SEX=0, ALB=45)

ev1 <- ev(time = 0, amt = 100, cmt = 1)
data.dose <- ev(ev1)
data.dose <- setDT(as.data.frame(data.dose))
data.all <- data.table(idata, data.dose)

outputsim <- modcovsim %>%
  data_set(data.all) %>%
  mrgsim(end = 24, delta = 0.25) %>%
  as.data.frame %>%
  as.data.table

outputsim$SEX <- factor(outputsim$SEX, labels="Female")

# Only plot a random sample of N=500
set.seed(678549)
plotdata <- outputsim[ID %in% sample(unique(ID), 500)]

# New facet label names for dose variable
albumin.labs <- c("albumin: 45 ng/mL")
names(albumin.labs) <- c("45")
wt.labs <- c("weight: 85 kg")
names(wt.labs) <- c("85")


p1 <- ggplot(plotdata, aes(time, CP, group = ID)) +
  geom_line(alpha = 0.2, size = 0.1) +
  facet_grid(~ WT + ALB + SEX,
             labeller =  labeller(ALB = albumin.labs,
                                  WT = wt.labs)) +
  labs(y = "Plasma Concentrations", x = "Time (h)")

p2 <- ggplot(plotdata, aes(time, CP, group = ID)) +
  geom_line(alpha = 0.2, size = 0.1) +
  facet_grid(~ WT + ALB + SEX,
             labeller =  labeller(ALB = albumin.labs,
                                  WT = wt.labs)) +
  scale_y_log10() +
  labs(y = "Plasma~Concentrations\n(logarithmic scale)", x = "Time (h)")

egg::ggarrange(p1, p2, ncol = 2)


## ----computenca, fig.width=7, message=FALSE-----------------------------------
derive.exposure <- function(time, CP) {
  n <- length(time)
  x <- c(
    Cmax = max(CP),
    Clast = CP[n],
    AUC = sum(diff(time) * (CP[-1] + CP[-n])) / 2
  )
  data.table(paramname=names(x), paramvalue=x)
}
refbsv <- outputsim[, derive.exposure(time, CP), by=.(ID, WT, SEX, ALB)]

p3 <- ggplot(refbsv, aes(
        x      = paramvalue,
        y      = paramname,
        fill   = factor(..quantile..),
        height = ..ndensity..)) +
  facet_wrap(~ paramname, scales="free", ncol=1) +
  stat_density_ridges(
    geom="density_ridges_gradient", calc_ecdf=TRUE,
    quantile_lines=TRUE, rel_min_height=0.001, scale=0.9,
    quantiles=c(0.05, 0.25, 0.5, 0.75, 0.95)) +
  scale_fill_manual(
    name   = "Probability",
    values = c("white", "#FF000050", "#FF0000A0", "#FF0000A0", "#FF000050", "white"),
    labels = c("(0, 0.05]", "(0.05, 0.25]",
               "(0.25, 0.5]", "(0.5, 0.75]",
               "(0.75, 0.95]", "(0.95, 1]")) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.y     = element_blank(),
    axis.ticks.y    = element_blank(),
    axis.title.y    = element_blank()) +
  labs(x="PK Parameters", y="") +
  scale_x_log10() +
  coord_cartesian(expand=FALSE)

# Obtain the standardized parameter value by dividing by the median.
refbsv[, stdparamvalue := paramvalue/median(paramvalue), by=paramname]

p4 <- ggplot(refbsv, aes(
        x      = stdparamvalue,
        y      = paramname,
        fill   = factor(..quantile..),
        height = ..ndensity..)) +
  facet_wrap(~ paramname, scales="free", ncol=1) +
  stat_density_ridges(
    geom="density_ridges_gradient", calc_ecdf=TRUE,
    quantile_lines=TRUE, rel_min_height=0.001, scale=0.9,
    quantiles=c(0.05, 0.25, 0.5, 0.75, 0.95)) +
  scale_fill_manual(
    name="Probability",
    values=c("white", "#FF000050", "#FF0000A0", "#FF0000A0", "#FF000050", "white"),
    labels = c("(0, 0.05]", "(0.05, 0.25]",
             "(0.25, 0.5]", "(0.5, 0.75]",
             "(0.75, 0.95]", "(0.95, 1]")) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.y     = element_blank(),
    axis.ticks.y    = element_blank(),
    axis.title.y    = element_blank()) +
  labs(x="Standardized PK Parameters", y="") +
  scale_x_log10() +
  coord_cartesian(expand=FALSE, xlim = c(0.3,3))

p3 + p4

## ----computebsvpk , fig.width=7 , message=FALSE-------------------------------
bsvranges <- refbsv[,list(
    P05 = quantile(stdparamvalue, 0.05),
    P25 = quantile(stdparamvalue, 0.25),
    P50 = quantile(stdparamvalue, 0.5),
    P75 = quantile(stdparamvalue, 0.75),
    P95 = quantile(stdparamvalue, 0.95)), by = paramname]
bsvranges

## ----covcomb, fig.width=7-----------------------------------------------------
reference.values <- data.frame(WT = 85, ALB = 45, SEX = 0)   
covdatasim$SEX<- ifelse(covdatasim$SEX==0,1,0)
covdatasim$SEX <- as.factor(covdatasim$SEX )
covdatasim$SEX <- factor(covdatasim$SEX,labels = c("Female","Male"))
covdatasimpairs <- covdatasim
covdatasimpairs$Weight <- covdatasimpairs$WT
covdatasimpairs$Sex <- covdatasimpairs$SEX
covdatasimpairs$Albumin <- covdatasimpairs$ALB

ggpairsplot <- GGally::ggpairs(covdatasimpairs,
                       columns = c("Weight","Sex","Albumin"),mapping = aes(colour=SEX),
                       diag= list(
                         continuous = GGally::wrap("densityDiag", alpha = 0.3,colour=NA),
                         discrete   = GGally::wrap("barDiag",  alpha =0.3, position = "dodge2")
                       ),
                       lower = list(
                         continuous = GGally::wrap("points", alpha = 0.2, size = 2),
                         combo = GGally::wrap("facethist", alpha =
                                                0.2, position = "dodge2")
                       ),
                       upper = list(
                         continuous = GGally::wrap("cor", size = 4.75, align_percent = 0.5),
                         combo = GGally::wrap("box_no_facet", alpha =0.3),
                         discrete = GGally::wrap("facetbar",  alpha = 0.3, position = "dodge2")
                       )
)
covdatasim$SEX <- as.numeric(covdatasim$SEX)-1
ggpairsplot + theme_bw(base_size = 12) +
  theme(axis.text = element_text(size=9))



## ----fig.width=7 ,message=FALSE, fig.height=5---------------------------------
idata <- data.table::copy(covdatasim)
idata$covname <-  NULL
ev1 <- ev(time=0, amt=100, cmt=1)
data.dose <- as.data.frame(ev1)
data.all <- data.table(idata, data.dose)

outcovcomb<- modcovsim %>%
  data_set(data.all) %>%
  zero_re() %>% 
  mrgsim(end=24, delta=0.25) %>%
  as.data.frame %>%
  as.data.table

outcovcomb$SEX <- as.factor(outcovcomb$SEX )
outcovcomb$SEX <- factor(outcovcomb$SEX, labels=c("Female", "Male"))

stat_sum_df <- function(fun, geom="ribbon", ...) {
  stat_summary(fun.data = fun, geom = geom, ...)
}
stat_sum_df_line <- function(fun, geom="line", ...) {
  stat_summary(fun.data = fun, geom = geom, ...)
}

fwt <- function(x, xcat, which, what, from, to, ...) {
  what <- sub("WT", "\nWeight", what)
  sprintf("%s %s [%s to %s[",
          which, what, signif_pad(from, 3, FALSE), signif_pad(to, 3, FALSE))
}

f <- function(x, xcat, which, what, from, to, ...) {
  what <- sub("ALB", "\nAlbumin", what)
  sprintf("%s %s [%s to %s[",
          which, what, signif_pad(from, 3, FALSE), signif_pad(to, 3, FALSE))
}

plotlines<- ggplot(outcovcomb, aes(time,CP,col=SEX ) )+
  geom_line(aes(group=ID),alpha=0.1,size=0.1)+
  facet_grid(table1::eqcut(ALB,2,f) ~ table1::eqcut(WT,2,fwt),labeller = label_value)+
  labs(colour="Sex",caption ="Simulation without Uncertainty\nFull
       Covariate Distribution\nwithout BSV/Uncertainty",
       x = "Time (h)", y="Plasma Concentrations")+
  coord_cartesian(ylim=c(0,3.5))

plotranges<- ggplot(outcovcomb, aes(time,CP,col=SEX,fill=SEX ) )+
  stat_sum_df(fun="median_hilow",alpha=0.2,
              mapping = aes(group=interaction(table1::eqcut(WT,2,fwt),
                                              SEX,
                                              table1::eqcut(ALB,2,f))
              ), colour = "transparent")+
  stat_sum_df_line(fun="median_hilow",size =2,
                   mapping = aes(linetype = SEX,
                                 group=interaction(table1::eqcut(WT,2,fwt),
                                                   SEX,table1::eqcut(ALB,2,f))))+
  facet_grid(table1::eqcut(ALB,2,f) ~ table1::eqcut(WT,2,fwt),
             labeller = label_value)+
  labs(linetype="Sex",colour="Sex",fill="Sex",
  caption ="Simulation with Full Covariate Distribution with BSV
  95% (Covariate Effects + BSV) Percentiles",
       x = "Time (h)", y="Plasma Concentrations")+
  coord_cartesian(ylim=c(0,3.5))
plotranges


## ----fig.width=7--------------------------------------------------------------
theta <- unclass(as.list(param(modcovsim)))
theta[c("WT", "SEX", "ALB")] <- NULL
theta <- unlist(theta)
as.data.frame(t(theta))

varcov <- cor2cov(
  matrix(0.2, nrow=length(theta), ncol=length(theta)),
  sd=theta*0.25)
rownames(varcov) <- colnames(varcov) <- names(theta)
as.data.frame(varcov)

## ----fig.width=7--------------------------------------------------------------
set.seed(678549)
# mvtnorm::rmvnorm is another option that can be explored
sim_parameters <- MASS::mvrnorm(nsim, theta, varcov, empirical=T) %>% as.data.table
head(sim_parameters)

## ----fig.width=7,fig.height=4, fig.height=5-----------------------------------
idata <- copy(covdatasim)
ev1       <- ev(time=0, amt=100, cmt=1)
data.dose <- as.data.frame(ev1)

iter_sims <- NULL
for(i in 1:nsim) {
  # you might want to resample your covariate database here
  # e.g. subsample from a large pool of patient
  # include uncertainty on your covariate distribution
  
  data.all  <- data.table(idata, data.dose, sim_parameters[i])
  out <- modcovsim %>%
    data_set(data.all) %>%
    #zero_re() %>%
    #omat(rxode2::cvPost(2000, matrix(c(0.09,0.01,0.01,0.09), 2, 2),
    #type = "invWishart")) %>%  # unc on bsv uncomment and increase nsim for CPT:PSP paper
    mrgsim(start=0, end=24, delta=0.25) %>%
    as.data.frame %>%
    as.data.table
  out[, rep := i]
  iter_sims <- rbind(iter_sims, out)
}
f <- function(x, xcat, which, what, from, to, ...) {
  what <- sub("ALB", "\nAlbumin", what)
  sprintf("%s %s [%s to %s[",
          which, what, signif_pad(from, 3, FALSE), signif_pad(to, 3, FALSE))
}
fwt <- function(x, xcat, which, what, from, to, ...) {
  what <- sub("WT", "\nWeight", what)
  sprintf("%s %s [%s to %s[",
          which, what, signif_pad(from, 3, FALSE), signif_pad(to, 3, FALSE))
}

iter_sims_summary_all <- iter_sims %>%
  mutate(WT=table1::eqcut(WT,2,fwt),ALB=table1::eqcut(ALB,2,f)) %>% 
  group_by(time,WT,ALB,SEX)%>%
  summarize( P50= median(CP) ,
             P05 = quantile(CP,0.05),
             P95= quantile(CP,0.95))


iter_sims_summary_all$SEX <- as.factor(iter_sims_summary_all$SEX )
iter_sims_summary_all$SEX <- factor(iter_sims_summary_all$SEX,labels = c("Female","Male"))
legendlabel<- "Median\n5%-95%"
plotrangesunc<- ggplot(iter_sims_summary_all,
                       aes(time,P50,col=SEX,fill=SEX,group=SEX,linetype=SEX) )+
  geom_ribbon(aes(ymin=P05,ymax=P95),alpha=0.3,linetype=0)+
  geom_line(size=1)+
  facet_grid(ALB ~ WT, labeller = label_value)+
  labs(linetype=legendlabel,colour=legendlabel,fill=legendlabel,
       caption ="Simulation with joint correlated covariate distributions
       with uncertainty and between subject variability",
       x = "Time (h)", y="Plasma Concentrations")+
  coord_cartesian(ylim=c(0,3.5))

plotrangesunc+
  theme(axis.title.y = element_text(size=12))+
  coord_cartesian(ylim=c(0,4))


## ----fig.width=7, include=TRUE, message=FALSE---------------------------------
out.df.parameters <- iter_sims[, derive.exposure(time, CP),
                                    by=.(rep, ID, WT, SEX, ALB)]

refvalues <- out.df.parameters[,.(medparam = median(paramvalue)), by=.(paramname,rep)]


## ----fig.width=7,fig.height=5 ,message=FALSE----------------------------------

setkey(out.df.parameters, paramname, rep)
out.df.parameters <- merge(out.df.parameters,refvalues)
out.df.parameters[, paramvaluestd := paramvalue/medparam]

out.df.parameters[, SEXCAT := ifelse( SEX==0,"Female","Male")]
out.df.parameters[, REF := "All Subjects"]
out.df.parameters[, WTCAT4 := table1::eqcut( out.df.parameters$WT,4,varlabel = "Weight")]
out.df.parameters[, ALBCAT3 := table1::eqcut( out.df.parameters$ALB,3,varlabel = "Albumin")]

nca.summaries.long <-  melt(out.df.parameters, measure=c("REF","WTCAT4","ALBCAT3","SEXCAT"),
                            value.name = "covvalue",variable.name ="covname" )

nca.summaries.long$covvalue <- as.factor( nca.summaries.long$covvalue)
nca.summaries.long$covvalue <- reorder(nca.summaries.long$covvalue,nca.summaries.long$paramvalue)

nca.summaries.long$covvalue <- factor(nca.summaries.long$covvalue,
                                      levels =c(  
                                        "1st tertile of Albumin: [31.0,44.0)"
                                        , "2nd tertile of Albumin: [44.0,46.0)"
                                        , "3rd tertile of Albumin: [46.0,54.0]"
                                        , "Male"  
                                        , "Female"
                                        , "All Subjects"
                                        , "1st quartile of Weight: [40.6,71.3)"
                                        , "2nd quartile of Weight: [71.3,85.0)"
                                        , "3rd quartile of Weight: [85.0,98.2)"
                                        ,"4th quartile of Weight: [98.2,222]"
                                      ))

nca.summaries.long$covvalue2 <- factor(nca.summaries.long$covvalue,
                                       labels =c(  
                                         "T1: [31.0,44.0)"
                                         , "T2: [44.0,46.0)"
                                         , "T3: [46.0,54.0]"
                                         , "Male"  
                                         , "Female"
                                         , "All Subjects"     
                                         , "Q1: [40.6,71.3)"
                                         , "Q2: [71.3,85.0)"
                                         , "Q3: [85.0,98.2)"
                                         , "Q4: [98.2,222]"
                                       ))

nca.summaries.long$covname<- as.factor(nca.summaries.long$covname)
nca.summaries.long$covname<- factor(nca.summaries.long$covname,
                                    levels =c("WTCAT4","SEXCAT","ALBCAT3","REF"),
                                    labels = c("Weight","Sex","Albumin","REF"))
func <- function(bob) c(min(bob), median(bob), max(bob))
boxplotMV<- ggplot(nca.summaries.long
                   , aes(x=covvalue2  , y=paramvalue ))+
  facet_grid (  paramname ~covname, scales="free", labeller=label_parsed,
                switch="both",space="free_x") +
  geom_boxplot(outlier.shape = NA) +
  theme_bw(base_size = 12)+
  theme(axis.title=element_blank(),
        strip.placement = "outside",
        axis.text.x = element_text(angle=20,vjust = 1, hjust = 1, face = "bold"),
        strip.text.y.left = element_text(angle= 0,vjust = 1, hjust = 1,face = "bold"))+
  scale_y_continuous(breaks = scales::pretty_breaks(n=4)  )
boxplotMV

## ----fig.width=7, fig.height=4 ,message=FALSE---------------------------------
ggridgesplot<- ggplot(nca.summaries.long,
                      aes(x=paramvaluestd,y=covvalue,
                          fill=factor(..quantile..),
                          height=..ndensity..))+
  facet_grid(covname~paramname,scales="free_y")+
  annotate( "rect",
            xmin = 0.8,
            xmax = 1.25,
            ymin = -Inf,
            ymax = Inf,
            fill = "gray",alpha=0.4
  )+
  stat_density_ridges(
    geom = "density_ridges_gradient", calc_ecdf = TRUE,
    quantile_lines = TRUE, rel_min_height = 0.01,scale=0.9,
    quantiles = c(0.05,0.5, 0.95))+
  geom_vline( aes(xintercept = 1),size = 1)+
  scale_fill_manual(
    name = "Probability", values = c("white","#0000FFA0", "#0000FFA0", "white"),
    labels = c("(0, 0.05]", "(0.05, 0.5]","(0.5, 0.95]", "(0.95, 1]")
  )+
  geom_vline(data=data.frame (xintercept=1),  aes(xintercept =xintercept  ),size = 1)+
  theme_bw()+
  theme(legend.position = "none")+
  labs(x="Effects Of Covariates on PK Parameter",y="")+
  scale_x_continuous(breaks=c(0.5,0.8,1/0.8,1/0.5,1/0.25),trans ="log" )+
  coord_cartesian(xlim=c(0.25,3))
ggridgesplot


## ----fig.width=7, fig.height=6------------------------------------------------
coveffectsdatacovrep <- nca.summaries.long %>% 
  dplyr::group_by(paramname,covname,covvalue) %>% 
  dplyr::summarize(
    mid= median(paramvaluestd),
    lower= quantile(paramvaluestd,0.05),
    upper = quantile(paramvaluestd,0.95)) %>% 
  dplyr::filter(!is.na(mid)) 


coveffectsdatacovrepbsv <- coveffectsdatacovrep[coveffectsdatacovrep$covname=="REF",]
coveffectsdatacovrepbsv$covname <- "BSV"
coveffectsdatacovrepbsv$covvalue <- "90% of patients"
coveffectsdatacovrepbsv$label <-    "90% of patients"
coveffectsdatacovrepbsv$lower <- bsvranges$P05
coveffectsdatacovrepbsv$upper <- bsvranges$P95
coveffectsdatacovrepbsv2 <- coveffectsdatacovrep[coveffectsdatacovrep$covname=="REF",]
coveffectsdatacovrepbsv2$covname <- "BSV"
coveffectsdatacovrepbsv2$covvalue <- "50% of patients"
coveffectsdatacovrepbsv2$label <-    "50% of patients"
coveffectsdatacovrepbsv2$lower <- bsvranges$P25
coveffectsdatacovrepbsv2$upper <- bsvranges$P75

coveffectsdatacovrepbsv<- rbind(coveffectsdatacovrep,coveffectsdatacovrepbsv2,
                                coveffectsdatacovrepbsv)
coveffectsdatacovrepbsv <- coveffectsdatacovrepbsv %>% 
  mutate(
    label= covvalue,
    LABEL = paste0(format(round(mid,2), nsmall = 2),
                   " [", format(round(lower,2), nsmall = 2), "-",
                   format(round(upper,2), nsmall = 2), "]"))
coveffectsdatacovrepbsv<- as.data.frame(coveffectsdatacovrepbsv)

coveffectsdatacovrepbsv$label <- gsub(": ", ":\n", coveffectsdatacovrepbsv$label)

coveffectsdatacovrepbsv$covname <-factor(as.factor(coveffectsdatacovrepbsv$covname ),
                                         levels = c("Weight","Sex","Albumin","REF","BSV"))
coveffectsdatacovrepbsv$label <- factor(coveffectsdatacovrepbsv$label,
                                        levels =c( "1st tertile of Albumin:\n[31.0,44.0)"
                                                   , "2nd tertile of Albumin:\n[44.0,46.0)"
                                                   , "3rd tertile of Albumin:\n[46.0,54.0]"
                                                   , "Male", "Female"
                                                   , "All Subjects","90% of patients","50% of patients"
                                                   , "1st quartile of Weight:\n[40.6,71.3)"
                                                   , "2nd quartile of Weight:\n[71.3,85.0)"
                                                   , "3rd quartile of Weight:\n[85.0,98.2)"
                                                   ,"4th quartile of Weight:\n[98.2,222]"
                                        ))
coveffectsdatacovrepbsv$label <- factor(coveffectsdatacovrepbsv$label,
                                        labels =c("T1:\n[31.0,44.0)"
                                                  , "T2:\n[44.0,46.0)"
                                                  , "T3:\n[46.0,54.0]"
                                                  , "Male", "Female"
                                                  , "All Subjects","90% of patients","50% of patients"
                                                  , "Q1:\n[40.6,71.3)"
                                                  , "Q2:\n[71.3,85.0)"
                                                  , "Q3:\n[85.0,98.2)"
                                                  , "Q4:\n[98.2,222]"
                                        ))

interval_legend_text <- "Median (points)\n90% intervals (horizontal lines) of joint effects:
covariate distributions, uncertainty
and between subject variability"
interval_bsv_text    <- "BSV (points)\nPrediction Intervals (horizontal lines)"
ref_legend_text      <- "Reference (vertical line)\nClinically relevant limits\n(gray area)"
area_legend_text     <- "Reference (vertical line)\nClinically relevant limits\n(gray area)"


#emf("Figure_PKdist_forest.emf",width= 15, height = 7.5)
png("./coveffectsplot_full.png",width =9.5 ,height = 8,units = "in",res=72)

forest_plot(coveffectsdatacovrepbsv[coveffectsdatacovrepbsv$covname!="REF"&
                                                      coveffectsdatacovrepbsv$covname!="BSV",
],
ref_area = c(0.8, 1/0.8),x_range = c(0.4,3),
strip_placement = "outside",base_size = 18,
y_label_text_size = 9,x_label_text_size = 10,
xlabel = "Fold Change Relative to Reference",
ref_legend_text =ref_legend_text,
area_legend_text =ref_legend_text ,
interval_legend_text = interval_legend_text,
plot_title = "",
interval_bsv_text = interval_bsv_text,
facet_formula = "covname~paramname",
facet_switch = "y",
table_facet_switch = "both",
reserve_table_xaxis_label_space = FALSE,
facet_scales = "free_y", facet_space = "free",
paramname_shape = FALSE,
table_position = "below",
show_table_yaxis_tick_label = TRUE,
table_text_size= 4,
plot_table_ratio = 1,
show_table_facet_strip = "both",
logxscale = TRUE,
major_x_ticks = c(0.5,0.8,1/0.8,1/0.5),
return_list = FALSE)
dev.off()


