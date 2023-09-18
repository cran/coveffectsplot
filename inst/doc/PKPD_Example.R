## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
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
library(patchwork)
library(ggh4x)
library(data.table)
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
cor2cov <- function (cor, sd) 
{
    if (missing(sd)) {
        sd <- diag(cor)
    }
    diag(cor) <- 1
    n <- nrow(cor)
    diag(sd, n) %*% cor %*% diag(sd, n)
}
nbsvsubjects <- 100
nsim <- 100 # uncertainty replicates for vignette you might want a higher number
round_pad <- function(x, digits = 2, round5up = TRUE) {
  eps <- if (round5up) x * (10^(-(digits + 3))) else 0
  formatC(round(x + eps, digits), digits = digits, format = "f", flag = "0")
}

## ----pkpdmodel, collapse=TRUE-------------------------------------------------
codepkpdmodelcov <- '
$PARAM @annotated
KA     : 0.5   : Absorption rate constant Ka (1/h)
CL     : 4     : Clearance CL (L/h)
V      : 10    : Central volume Vc (L)
Vp     : 50    : Peripheral volume Vp (L)
Qp     : 10    : Intercompartmental clearance Q (L/h)
CLALB  : -0.8  : Ablumin on CL (ref. 45 g/L)
CLSEX  : 0.2   : Sex on CL (ref. Female)
CLWT   : 1     : Weight on CL (ref. 85 kg)
VSEX   : 0.07  : Sex on Vc (ref. Female)
VWT    : 1     : Weight on Vc (ref. 85 kg)
KIN    : 3     : Zero-order Rate constant of biomarker production (amount/h)
KOUT   : 0.06  : First-order Rate constant of biomarker loss (1/h)
IC50   : 3     : Drug concentration producing 50% of maximum inhibition
IMAX   : 0.999 : Maximum Inhibition Response
gamma  : 0.55  : Sigmoidicity factor of the sigmoid Emax equation
KINWT  : 0.4   : Weight on KIN (ref. 85 kg)
KINAGE : -0.08 : Age on KIN (ref. 40 years)
KINHLTY: 1.5   : Weight on CL (ref. 85 kg)

$PARAM @annotated // reference values for covariate
WT     :  85    : Weight (kg)
SEX    :  0     : Sex (0=Female, 1=Male)
ALB    :  45    : Albumin (g/L)
AGE    :  40    : Age (years)
HEALTHY:  0     : Health Status (0=Diseased, 1=Healthy)

$CMT GUT CENT PER RESP
$GLOBAL
#define CP   (CENT/Vi)
#define CPER (PER/Vpi)
#define INH  (IMAX*pow(CP,gamma)/(pow(IC50,gamma)+pow(CP,gamma)))
#define PDRESP RESP

$MAIN
double KAi = KA;
double Vpi = Vp *pow((WT/70.0),    1);
double Qpi = Qp *pow((WT/70.0), 0.75);
double CLi = CL *
    pow((ALB/45.0), CLALB)*
    (SEX == 1.0 ? (1.0+CLSEX) : 1.0)*
    pow((WT/85.0), CLWT)*exp(ETA(1)); 
double Vi = V *
    (SEX == 1.0 ? (1.0+VSEX) : 1.0)*
    pow((WT/85.0), VWT)*exp(ETA(2));  
double KINi = KIN *
  pow((AGE/40), KINAGE)*
  (HEALTHY == 1.0 ? KINHLTY : 1.0)*
  pow((WT/85.0), KINWT)*exp(ETA(3));
double RESP_0 = KINi/KOUT;

$OMEGA
0.09 
0.01 0.09
$OMEGA
0.25

$ODE
dxdt_GUT    = -KAi *GUT;
dxdt_CENT   =  KAi *GUT  - (CLi+Qpi)*CP  + Qpi*CPER;
dxdt_PER    =                   Qpi*CP   - Qpi*CPER;
dxdt_RESP   =  KINi*(1-INH) - KOUT*RESP;

$CAPTURE CP PDRESP KAi CLi Vi Vpi Qpi WT SEX ALB AGE HEALTHY
'
modpkpdsim <- mcode("codepkpdmodelcov", codepkpdmodelcov)
partab <- setDT(modpkpdsim@annot$data)[block=="PARAM", .(name, descr, unit)]
partab <- merge(partab, melt(setDT(modpkpdsim@param@data), meas=patterns("*"), var="name"))
knitr::kable(partab)


## ----pkpdsimulation, fig.width=7,fig.height=4, message=FALSE------------------
idata <- data.table(ID=1:nbsvsubjects, WT=85, SEX=0, ALB=45, AGE=40, HEALTHY = 0)
ev1 <- ev(time = 0, amt = 100, cmt = 1, ii = 24, addl = 20)
data.dose <- ev(ev1)
data.dose <- setDT(as.data.frame(data.dose))
data.all <- data.table(idata, data.dose)

set.seed(678549)
outputpkpdsim <- modpkpdsim %>%
  data_set(data.all) %>%
  mrgsim(end = 28*24, delta = 0.25) %>%
  as.data.frame %>%
  as.data.table

outputpkpdsim$HEALTHY <- as.factor(outputpkpdsim$HEALTHY)

yvar_names <- c(
  'CP'="Plasma Concentrations",
  'RESP'="PD Values"
)
set.seed(678549)
outputpkpdsimlong <- outputpkpdsim[outputpkpdsim$ID %in%
sample(unique(outputpkpdsim$ID), 5), ] %>% 
  gather(key,value,CP,RESP)

ggplot(data =outputpkpdsimlong ,
       aes(time, value, group = ID)) +
  geom_line(alpha = 0.8, size = 0.3) +
  facet_grid(key ~ID,scales="free_y",switch="y",
             labeller = labeller(key=yvar_names)) +
  labs(y = "", color = "Sex", x = "Time (h)")+
  theme(strip.placement = "outside",
        axis.title.y=element_blank())

## ----computenca , fig.width=7,  message=FALSE---------------------------------
derive.exposure <- function(time, PDRESP) {
  x <- c(
    nadir = min(PDRESP, na.rm = TRUE),
    baselinepd = PDRESP[1L],
    deltapd = PDRESP[1L]-min(PDRESP, na.rm = TRUE)
  )
  data.table(paramname=names(x), paramvalue=x)
}
refbsv <- outputpkpdsim[, derive.exposure(time, PDRESP),
                        by=.(ID, WT, SEX, ALB, AGE, HEALTHY)]

refbsv[, stdparamvalue := paramvalue/median(paramvalue), by=paramname]

bsvranges <- refbsv[,list(
    P05 = quantile(stdparamvalue, 0.05),
    P25 = quantile(stdparamvalue, 0.25),
    P50 = quantile(stdparamvalue, 0.5),
    P75 = quantile(stdparamvalue, 0.75),
    P95 = quantile(stdparamvalue, 0.95)), by = paramname]
bsvranges

## ----covcomb , fig.width=7----------------------------------------------------
reference.values <- data.frame(WT = 85, ALB = 45, AGE = 40, SEX = 0, HEALTHY = 0)   
covcomb <- expand.modelframe(
  WT  = c(56,128), 
  AGE = c(20,60),
  ALB = c(40,50),
  SEX = c(1),#Refernce is for SEX =0
  HEALTHY = c(1),#Refernce is for HEALTHY =0
  rv = reference.values)

# Add the reference
covcomb <- rbind(covcomb, data.table(reference.values, covname="REF"))
covcomb$ID <- 1:nrow(covcomb)

covcomb

## ---- fig.width=7 ,message=FALSE, include=FALSE-------------------------------
idata <- data.table::copy(covcomb)
idata$covname <-  NULL
ev1 <- ev(time=0, amt=100, cmt=1, ii = 24, addl = 20)
data.dose <- as.data.frame(ev1)
data.all <- data.table(idata, data.dose)

outcovcomb<- modpkpdsim %>%
  data_set(data.all) %>%
  zero_re() %>% 
  mrgsim(start=0,end=24*28,delta=0.25)%>%
  as.data.frame %>%
  as.data.table

outcovcomb$SEX <- as.factor(outcovcomb$SEX )
outcovcomb$SEX <- factor(outcovcomb$SEX, labels=c("Female", "Male"))

outcovcomb$HEALTHY <- as.factor(outcovcomb$HEALTHY )


theta <- unclass(as.list(param(modpkpdsim)))
theta[c("WT", "SEX", "ALB","AGE","HEALTHY")] <- NULL
theta <- unlist(theta)
as.data.frame(t(theta))

varcov <- cor2cov(
  matrix(0.2, nrow=length(theta), ncol=length(theta)),
  sd=theta*0.25)
rownames(varcov) <- colnames(varcov) <- names(theta)
as.data.frame(varcov)

set.seed(678549)
# mvtnorm::rmvnorm is another option that can be explored
sim_parameters <- MASS::mvrnorm(nsim, theta, varcov, empirical=T) %>% as.data.table
head(sim_parameters)

idata <- data.table::copy(covcomb)
idata$covname <-  NULL
ev1 <- ev(time=0, amt=100, cmt=1, ii = 24, addl = 20)
data.dose <- as.data.frame(ev1)

iter_sims <- NULL
for(i in 1:nsim) {
  data.all  <- data.table(idata, data.dose, sim_parameters[i])
  out <- modpkpdsim %>%
    data_set(data.all) %>%
    zero_re() %>% 
    mrgsim(start=0, end=28*24, delta=0.25) %>%
    as.data.frame %>%
    as.data.table
  out[, rep := i]
  iter_sims <- rbind(iter_sims, out)
}
iter_sims$SEX <- as.factor(iter_sims$SEX )
iter_sims$SEX <- factor(iter_sims$SEX, labels=c("Female", "Male"))


## ---- fig.width=7, fig.height=6, message=FALSE, warning=FALSE-----------------
albumin.labs <- c("albumin: 40 ng/mL","albumin: 45 ng/mL","albumin: 50 ng/mL")
names(albumin.labs) <- c("40","45","50")
wt.labs <- c("weight: 85 kg","weight: 56 kg","weight: 128 kg")
names(wt.labs) <- c("85","56","128")

age.labs <- c("age: 20 years","age: 40 years","age: 60 years")
names(age.labs) <- c("20","40","60")

pdprofiles <- ggplot(iter_sims[iter_sims$rep<=10,], aes(time/24,PDRESP,col=factor(WT),linetype=factor(HEALTHY) ) )+
  geom_line(aes(group=interaction(ID,rep)),alpha=0.3,size=0.3)+
  geom_line(data=outcovcomb,aes(group=interaction(ID)),color="black")+
  facet_nested(ALB+SEX~ AGE+WT,  labeller = 
                 labeller( WT = wt.labs,
                           ALB = albumin.labs,
                           AGE = age.labs))+
  labs(linetype="Black Lines\nNo Uncertainty\nHealthy Status",
       colour="Colored Lines\nUncertainty\nReplicates\n(1 to 10)\nWeight (kg)",
       caption ="Simulation\nwith Uncertainty without BSV" ,
       x="Days", y = "PD Values")+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))

pdprofiles

## ---- fig.width=7,fig.height=6, include=FALSE, message=FALSE------------------

out.df.univariatecov.nca <- iter_sims[, derive.exposure(time, PDRESP),
                                      by=.(rep, ID, WT, SEX, ALB, AGE, HEALTHY)]
out.df.univariatecov.nca

refvalues <- out.df.univariatecov.nca[
  ALB==45 & WT==85 & SEX=="Female"& AGE==40 & HEALTHY==0,
  .(medparam = median(paramvalue)), by=paramname]
refvalues

covcomb <- as.data.table(covcomb)

covcomb[covname=="WT",  covvalue := paste(WT,"kg")]
covcomb[covname=="ALB", covvalue := paste(ALB,"g/L")]
covcomb[covname=="AGE", covvalue := paste(AGE,"years")]
covcomb[covname=="SEX", covvalue := "Male"]
covcomb[covname=="HEALTHY", covvalue := "Diseased"]
covcomb[covname=="REF", covvalue := "85 kg-Female-45 g/L-40 years-healthy"]
covcomb
covcomb[covname=="REF", covvalue := "85 kg-Female\n45 g/L-40 years\nhealthy"]
covcomb <- as.data.table(covcomb)
out.df.univariatecov.nca <- merge(
  out.df.univariatecov.nca,
  covcomb[, .(ID, covname, covvalue)])

setkey(out.df.univariatecov.nca, paramname)
setkey(refvalues, paramname)
out.df.univariatecov.nca <- merge(out.df.univariatecov.nca,refvalues)

out.df.univariatecov.nca[, paramvaluestd := paramvalue/medparam]

out.df.univariatecov.nca$covvalue <-factor(as.factor(out.df.univariatecov.nca$covvalue ),
                                               levels =  c("56 kg",
                                                           "85 kg",
                                                           "128 kg",
                                                           "Male",
                                                           "40 g/L",
                                                           "50 g/L",
                                                           "20 years",
                                                           "60 years",
                                                           "Diseased",
"85 kg-Female\n45 g/L-40 years\nhealthy")
)
out.df.univariatecov.nca$covname2 <- as.factor(out.df.univariatecov.nca$covname2)
out.df.univariatecov.nca$covname2 <- factor(out.df.univariatecov.nca$covname2,
                                            levels= c( "Weight", "Sex",
                                                    "Albumin","Age", "Healthy",
                                                    "Reference")
                                            )
boxplotdat <- out.df.univariatecov.nca
boxplotdat[covname=="WT",  covname2 := "Weight"]
boxplotdat[covname=="ALB", covname2 := "Albumin"]
boxplotdat[covname=="SEX", covname2 := "Sex"]
boxplotdat[covname=="AGE", covname2 := "Age"]
boxplotdat[covname=="HEALTHY", covname2 := "Healthy"]
boxplotdat[covname=="REF", covname2 := "Reference"]


## ---- fig.width=7, fig.height=5, message=FALSE--------------------------------
boxplotpd <- ggplot(boxplotdat,
       aes(x=covvalue ,y=paramvalue))+
  facet_grid(paramname ~covname2,scales="free",switch="both",
             labeller = label_parsed)+
    geom_boxplot()+
  theme(axis.title = element_blank(),strip.placement = "outside")+
  labs(y="PD Parameter Values",x="Covariate Value")
boxplotpd

## ---- fig.width=7 ,message=FALSE, include=FALSE-------------------------------
#  pdprofiles<- pdprofiles+theme(axis.title.y = element_text(size=15))+
#    guides(colour=guide_legend(override.aes = list(alpha=1,size=0.5)),
#           linetype=guide_legend(override.aes = list(size=0.5)))
# pdprofiles
# ggsave("pd3.png", device="png",type="cairo-png",width= 7, height = 5,dpi=72)
#  boxplotpd
#  ggsave("pd4.png", device="png",type="cairo-png",width= 7, height = 4,dpi=2*72)
#   png("Figure_S_PD_2.png", type="cairo-png",width= 2*7*72, height =5*72)
#  egg::ggarrange(pdprofiles,boxplotpd,nrow=1)
#  dev.off()

## ---- fig.width=7,fig.height=5,message=FALSE----------------------------------
pdggridges<- ggplot(out.df.univariatecov.nca,
       aes(x=paramvaluestd,y=covvalue,fill=factor(..quantile..),height=..ndensity..))+
  facet_grid(covname2~paramname,scales="free_y",space="free")+
  annotate( "rect",
            xmin = 0.5,
            xmax = 2,
            ymin = -Inf,
            ymax = Inf,
            fill = "gray",alpha=0.4
  )+
  stat_density_ridges(
    geom = "density_ridges_gradient", calc_ecdf = TRUE,
    quantile_lines = TRUE, rel_min_height = 0.001,scale=0.9,
    quantiles = c(0.05,0.5, 0.95)) +
  scale_fill_manual(
    name = "Probability", values = c("white", "#0000FFA0","#0000FFA0", "white"),
    labels = c("(0, 0.05]", "(0.05, 0.5]","(0.5, 0.95]", "(0.95, 1]")
  )+
  geom_vline( aes(xintercept = 1),size = 1)+
  theme_bw()+
  labs(x="Effects Relative to Parameter Reference value",y="")+
  scale_x_continuous(breaks=c(0.25,0.5,0.8,1/0.8,1/0.5,1/0.25))+
  scale_x_log10()
pdggridges

## ---- fig.width=7 ,message=FALSE, include=FALSE-------------------------------
# pdggridges+theme(legend.position = "none")
#  ggsave("Figure_S_PD_3.png", device="png",type="cairo-png",
#         width= 7, height = 5,dpi=72)


## ---- fig.width=7, fig.height=7 ,message=FALSE--------------------------------
coveffectsdatacovrep <- out.df.univariatecov.nca %>% 
  dplyr::group_by(paramname,covname,covvalue) %>% 
  dplyr::summarize(
    mid= median(paramvaluestd),
    lower= quantile(paramvaluestd,0.05),
    upper = quantile(paramvaluestd,0.95))

coveffectsdatacovreplabel<-   coveffectsdatacovrep %>%
  mutate(
    label= covvalue,
    LABEL = paste0(format(round(mid,2), nsmall = 2),
                   " [", format(round(lower,2), nsmall = 2), "-",
                   format(round(upper,2), nsmall = 2), "]"))


## ---- fig.width=7, fig.height=7 ,message=FALSE--------------------------------
setkey(bsvranges, paramname)
coveffectsdatacovrepbsv <- coveffectsdatacovrep[coveffectsdatacovrep$covname=="REF",]
coveffectsdatacovrepbsv$covname <- "BSV"
coveffectsdatacovrepbsv$covvalue <- "50% of patients"
coveffectsdatacovrepbsv$label <-    "50% of patients"
coveffectsdatacovrepbsv$lower <- bsvranges$P25
coveffectsdatacovrepbsv$upper <- bsvranges$P75

coveffectsdatacovrepbsv2 <- coveffectsdatacovrep[coveffectsdatacovrep$covname=="REF",]
coveffectsdatacovrepbsv2$covname <- "BSV"
coveffectsdatacovrepbsv2$covvalue <- "90% of patients"
coveffectsdatacovrepbsv2$label <-    "90% of patients"
coveffectsdatacovrepbsv2$lower <- bsvranges$P05
coveffectsdatacovrepbsv2$upper <- bsvranges$P95
coveffectsdatacovrepbsv<- rbind(coveffectsdatacovrep,coveffectsdatacovrepbsv,coveffectsdatacovrepbsv2)
coveffectsdatacovrepbsv <- coveffectsdatacovrepbsv %>% 
  mutate(
    label= covvalue,
    LABEL = paste0(format(round(mid,2), nsmall = 2),
                   " [", format(round(lower,2), nsmall = 2), "-",
                   format(round(upper,2), nsmall = 2), "]"))
coveffectsdatacovrepbsv<- as.data.frame(coveffectsdatacovrepbsv)
coveffectsdatacovrepbsv$label <- as.factor(coveffectsdatacovrepbsv$covvalue )
coveffectsdatacovrepbsv$label <- reorder(coveffectsdatacovrepbsv$label,
                                         coveffectsdatacovrepbsv$lower)


coveffectsdatacovrepbsv$covname <-factor(as.factor(coveffectsdatacovrepbsv$covname ),levels =c("WT","SEX","ALB","AGE","HEALTHY", "REF", "BSV"),
labels=  c("Weight","Sex","Albumin","Age","Healthy", "Reference", "BSV")
)

interval_legend_text <- "Median (points)\n90% CI (horizontal lines)"
interval_bsv_text <- "BSV (points)\nPrediction Intervals (horizontal lines)"
ref_legend_text <- "Reference\n(vertical line)\nClinically relevant limits\n(gray area)"
area_legend_text <- "Reference\n(vertical line)\nClinically relevant limits\n(gray area)"
png("./Figure_S_PD_4.png",width =12 ,height = 9,units = "in",res=72)
coveffectsplot::forest_plot(coveffectsdatacovrepbsv,
                            ref_area = c(0.5, 1/0.5),
                            x_range = c(0.25,4),
                            strip_placement = "outside",
                            base_size = 16,
                            y_label_text_size = 12,
                            y_label_text_width = 50,
                            xlabel = "Fold Change Relative to Reference",
                            ref_legend_text = ref_legend_text,
                            area_legend_text =area_legend_text,
                            interval_legend_text = interval_legend_text,
                            interval_bsv_text = interval_bsv_text,
                            facet_formula = "covname~paramname",
                            facet_switch = "y",
                            facet_scales = "free_y",
                            facet_space = "fixed",
                            paramname_shape = FALSE,
                            table_position = "below",
                            table_text_size=4,
                            plot_table_ratio = 1,
                            table_facet_switch = "both",
                            show_table_facet_strip = "both",
                            show_table_yaxis_tick_label = TRUE,
                            logxscale = TRUE,
                            major_x_ticks          = c(0.5, 1,  1/0.5),
                            major_x_labels         = c("1/2", "1", "2"),
                            table_margin = c(0,5.5,0,0),
                            plot_margin =c(0,5.5,0,0),
                            reserve_table_xaxis_label_space = FALSE,
                            return_list = FALSE)
dev.off()
# consider returning a list and editing the y axis label line breaks height
# theme(axis.text.y = element_text(lineheight = ))


