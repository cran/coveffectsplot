## ----setup, include = FALSE---------------------------------------------------
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
library(coveffectsplot)
library(ggplot2)
library(ggridges)
library(tidyr)
suppressPackageStartupMessages( library(dplyr) )
nuncertainty <- 10000
nbsvsubjects <- 100000


## ---- echo=TRUE, results='asis',fig.align = "center",fig.width = 6------------
set.seed(657687)
df <- data.frame(
MASS::mvrnorm(n = nuncertainty,
                mu = c(10,0.75,1.5),
                Sigma=matrix(c((10*0.15)^2,
                               0.001,0.001,0.001,(0.75*0.15)^2,
                               0.001,0.001,0.001,(1.5*0.15)^2),3,3,byrow = TRUE) 
))
names(df) <- c("POPCL","dWTdCL","dSexdCL")
knitr::kable(head(round(df,2),5))

## ---- echo=FALSE, results='asis',fig.align = "center",fig.width = 6-----------
dflong <- gather(df)

ggplot(dflong,aes(x=value))+
  geom_density(alpha=0.6, fill = "blue")+
    facet_wrap(~key,scales="free",ncol=1)+
    labs(fill="",x="Uncertainty Distribution (RSE 15%) of the Parameters")+
  theme(legend.position = "none",legend.background = 
                   element_rect(fill="transparent"),
                 axis.ticks.y = element_blank(),axis.text.y =element_blank())+
  guides(fill=guide_legend(reverse = FALSE))
set.seed(657687)
dfcov<- data.frame(
MASS::mvrnorm(n=nbsvsubjects,
                mu =c(65,75),
                Sigma=matrix(c(20^2,0.01,0.01,20^2),2,2,byrow = TRUE) 
))
names(dfcov)<- c("WTWOMAN","WTMAN")
dfcovlong <- gather(dfcov)
ggplot(dfcovlong,aes(x=value,linetype=key))+
  geom_density(alpha=0.2,fill="darkgreen",size=2)+
  labs(linetype="",x="Weight (kg)")+
  theme(legend.position = "right",legend.background = 
                   element_rect(fill="transparent"),axis.ticks.y = element_blank(),axis.text.y =element_blank())+
  guides(fill=guide_legend(reverse = FALSE))+
  theme_bw()


dfcovlongquantile<- as.data.frame(
  round(quantile(dfcovlong$value,probs=c(0.01,0.05,0.25,0.5,0.75,0.95,0.99)),0)
)
names(dfcovlongquantile)<- "Weightquantilevalue"
dfcovlongquantile$quantile<- rownames(dfcovlongquantile)

dfcovlongquantiletable<- t(dfcovlongquantile)
knitr::kable(dfcovlongquantiletable[1,,drop=FALSE],row.names=FALSE)

## ---- echo=TRUE,fig.align = "center",fig.width = 6----------------------------
set.seed(546789)
CLBSVdistribution <- data.frame(CL= 10*exp(rnorm(nbsvsubjects,0,sd=0.09^0.5)))
CLBSVdistribution$CLBSV<- CLBSVdistribution$CL/10

## ---- echo=FALSE,fig.align = "center",fig.width = 6 ,fig.height=4-------------
dfbsv<- as.data.frame(
  round( quantile(CLBSVdistribution$CLBSV,probs=c(0.01,0.05,0.25,0.5,0.75,0.95,0.99)),2))
names(dfbsv)<- "BSVquantilevalue"
dfbsv$quantile<- rownames(dfbsv)
CLBSVdistribution$paramname<- "CL"
bsvplot<-   ggplot(CLBSVdistribution, aes(
  x      = CLBSV,
  y      = paramname,
  fill   = factor(..quantile..),
  height = ..ndensity..)) +
  stat_density_ridges(
    geom="density_ridges_gradient", calc_ecdf=TRUE,
    quantile_lines=TRUE, rel_min_height=0.001, scale=0.9,
    quantiles=c(0.05, 0.25, 0.5, 0.75, 0.95)) +
  scale_fill_manual(
    name="BSV Ranges",
    values=c("white", "#FF000050", "#FF0000A0", "#FF0000A0", "#FF000050", "white"),
    labels = c("(0, 0.05]", "(0.05, 0.25]",
               "(0.25, 0.5]", "(0.5, 0.75]",
               "(0.75, 0.95]", "(0.95, 1]")) +
  theme_bw(base_size = 16) +
  theme(
    legend.position = "right",
    axis.text.y     = element_blank(),
    axis.ticks.y    = element_blank(),
    axis.title.y    = element_blank(),
    axis.text.x=element_text(size=12),
    axis.title.x=element_text(size=14)) +
  scale_x_continuous(breaks=c(0.61,0.82,1,1.22,1.63))+
  coord_cartesian(expand=FALSE,xlim = c(0.49,2.01))+
  labs(x="Standardized Individual Clearances with BSV",
                title="Illustrating 30.7% BSV")

bsvplot

dfbsvtable<- t(dfbsv)
knitr::kable(dfbsvtable[1,,drop=FALSE],row.names=FALSE)
#bayestestR::hdi(CLBSVdistribution$CLBSV,ci = c(.50, .90))
#0.55, 0.74, 1 , 1.12, 1.53

## ----fig.width= 7-------------------------------------------------------------
dfeffects <- df
dfeffects$REFrange <- dfeffects$POPCL/ median(dfeffects$POPCL)

ggplot(dfeffects)+
  geom_density(aes(x=REFrange,y=..scaled..,col="a.Uncertainty\nRSE=15%"))+
  geom_density(data=dfcovlong,
                        aes(x=(value/70)^0.75 ,
                                y=..scaled..,col="b.Weight\nMean=70 kg, sd=20 kg"))+
  geom_density(data=CLBSVdistribution ,aes(x=CLBSV,y=..scaled..,
  col="c.Between subject variability\nCV=30.1%"))+
  theme_bw(base_size = 16)+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  scale_x_continuous(breaks=c(0.25,0.5,0.8,1,1.25,1.5,2,3))+
  coord_cartesian(xlim=c(0.25,2))+
  labs(color="",x="Effects Standardized Relative\nto the Typical Value",y= "Scaled Density")+
  scale_color_manual(values=c("blue","green","red"))


## ----fig.width= 7-------------------------------------------------------------
dfeffects$REF <- 1
dfeffects$SEX_FEMALE_WT_50 <- (50/70)^dfeffects$dWTdCL
dfeffects$SEX_FEMALE_WT_90 <-  (90/70)^dfeffects$dWTdCL
dfeffects$SEX_Male_WT_70 <- dfeffects$dSexdCL
dfeffects$SEX_Male_WT_90 <- dfeffects$dSexdCL*(90/70)^dfeffects$dWTdCL
dfeffects$BSV_REF<-  sample(CLBSVdistribution$CLBSV, nuncertainty)

dfeffects<- dfeffects[,c("SEX_FEMALE_WT_50",
                         "SEX_FEMALE_WT_90",
                         "SEX_Male_WT_70",
                         "SEX_Male_WT_90",
                         "REF",
                         "BSV_REF")]


dflong <- tidyr::gather(dfeffects,key,value,-REF)
ggplot(dflong,aes(x=value,y=key,fill=factor(..quantile..)))+
ggridges::stat_density_ridges(
  geom = "density_ridges_gradient", calc_ecdf = TRUE,
  quantile_lines = TRUE, rel_min_height = 0.01,
  quantiles = c(0.05,0.5, 0.95)) +
  scale_fill_manual(
    name = "Probability", values = c("#FF0000A0", "white","white", "#0000FFA0"),
    labels = c("(0, 0.05]", "(0.05, 0.5]","(0.5, 0.95]", "(0.95, 1]")
  )+
      annotate(
        "rect",
        xmin = 0.8,
        xmax = 1.25,
        ymin = -Inf,
        ymax = Inf,
        fill = "gray",alpha=0.4
      )+
  geom_vline(
      aes(xintercept = 1),
      size = 1
    )+
  theme_bw()+
  labs(x="Effects Relative to parameter reference value",y="")

## -----------------------------------------------------------------------------
dfeffects$SEX_Male_WT_90<- NULL
dfeffects$REF <- 1
dfeffectslong<- gather(dfeffects)
dfeffectslong<- dplyr::group_by(dfeffectslong,key)
dfeffectslongsummaries<- dplyr::summarise(dfeffectslong,mid=quantile(value,0.5),
                                   lower=quantile(value,0.05),
                                   upper=quantile(value,0.95))

dfeffectslongsummaries$paramname <- "CL"
dfeffectslongsummaries$covname <- c("BSV","REF","Weight","Weight","Sex")
dfeffectslongsummaries$label <- c("95% of patients","70 kg/Woman",
                                  "50 kg/Woman", "90 kg/Woman","70 kg/Man")
dfeffectslongsummaries<- rbind(dfeffectslongsummaries,
data.frame(key=c("BSV","BSV"),
           mid=c(quantile(dfeffects$BSV,0.5), quantile(dfeffects$BSV,0.5)),
           lower = c(quantile(dfeffects$BSV,0.25), quantile(dfeffects$BSV,0.05)),
            upper = c(quantile(dfeffects$BSV,0.75), quantile(dfeffects$BSV,0.95)),
           paramname= "CL",
           covname=c("BSV","BSV"),
           label = c("50% of patients","90% of patients")
)
)
dfeffectslongsummaries<- dfeffectslongsummaries[c(2,6,7,3,4,5),]

plotdata <- dplyr::mutate(dfeffectslongsummaries,
          LABEL = paste0(format(round(mid,2), nsmall = 2),
                         " [", format(round(lower,2), nsmall = 2), "-",
                         format(round(upper,2), nsmall = 2), "]"))
plotdata<- as.data.frame(plotdata)
plotdata<- plotdata[,c("paramname","covname","label","mid","lower","upper","LABEL")]
knitr::kable(plotdata)


## ----fig.width=7--------------------------------------------------------------
plotdata$covname <- as.factor(plotdata$covname)
plotdata$covname <- reorder(plotdata$covname , c(3,4,4,2,1,1))

plotdata$label <- reorder(as.factor(plotdata$label) , c(1,3,2,4,5,6))
  ggplot(data = plotdata[plotdata$covname!="REF",], aes_string(
      y = "label",
      x = "mid",
      xmin = "lower",
      xmax = "upper"
    )) +
      geom_pointrange(
      aes(color = "90 %CI\nCovariate Effects"),
      size = 1,
      alpha = 1
    )+
  annotate("rect", xmin = min(0.8), 
      xmax = max(1.25), ymin = -Inf, ymax = Inf, fill = "gray",alpha=0.1)+
  geom_vline(aes(xintercept = 1,linetype="Reference"))+ 
  facet_grid(covname~.,scales="free_y",switch="y")+
  labs(y="",x="Effects Relative to Reference Value",
                colour="",linetype="")+
  theme_bw()+
    scale_color_manual(values=c("blue"))

## ----dpi = 72-----------------------------------------------------------------
png("./coveffectsplot.png",width =9 ,height = 6,units = "in",res=72)
 coveffectsplot::forest_plot(plotdata[plotdata$covname!="REF",],
            ref_area = c(0.8, 1/0.8),
            x_facet_text_size = 13,
            y_facet_text_size = 13,
            interval_legend_text = "Median (points)\n90% CI (horizontal lines)",
            ref_legend_text = "Reference (vertical line)\n+/- 20% ratios (gray area)",
            area_legend_text = "Reference (vertical line)\n+/- 20% ratios (gray area)",
            xlabel = "Fold Change Relative to a 70kg Woman",
            facet_formula = "covname~.",
            facet_switch = "both",
            facet_scales = "free",
            facet_space = "fixed",
            paramname_shape = TRUE,
            show_table_facet_strip = "none",
            table_position = "right",
            table_text_size=4,
            plot_table_ratio = 4,
            legend_space_x_mult = 0.5,
            return_list = FALSE)
dev.off()

## ----fig.width=7 ,eval=FALSE--------------------------------------------------
#  
#    plotdata<- plotdata[ c(3,2,1,4,5,6),]
#    plotly::plot_ly(plotdata) %>%
#    plotly::add_segments(
#    x = ~ round(lower, 2),
#    xend = ~ round(upper, 2),
#    y = ~ label,
#    yend = ~ label,
#    name = '90%CI',
#    line = list(color = plotly::toRGB("blue", alpha = 0.5), width = 5),
#    hoverinfo = "text",
#    text = ~ paste("</br> 90%CI: ",
#    paste(round(lower, 2), round(upper, 2)))
#    ) %>%
#    plotly::add_markers(
#    x = ~ round(mid, 2),
#    y = ~ label,
#    name = "Median",
#    marker  = list(
#    color = plotly::toRGB("black", alpha = 0.3),
#    size = 20,
#    symbol = "diamond"
#    ),
#    hoverinfo = "text",
#    text = ~ paste("</br> Median: ",
#    paste(round(mid, 2)))
#    ) %>%
#    plotly::layout(
#    xaxis = list(
#    title = 'Effects Relative to Reference',
#    ticks = "outside",
#    autotick = TRUE,
#    ticklen = 5,
#    gridcolor = plotly::toRGB("gray50"),
#    showline = TRUE
#    ) ,
#    yaxis = list (
#    title = '' ,
#    autorange = TRUE,
#    type = "category",
#    categoryorder = "trace",
#    ticks = "outside",
#    autotick = TRUE,
#    ticklen = 5,
#    gridcolor = plotly::toRGB("gray50"),
#    showline = TRUE
#    ),
#       shapes =list(
#        type = "rect",
#        x0 = 0.8,
#        x1 = 1.25,
#        xref = "x",
#        yref = "paper",
#        y0 = 0,
#        y1 = 1,
#        line = list(width = 0),
#        fillcolor =  plotly::toRGB("black", alpha = 0.2)
#    )
#    )

## ----dpi = 72-----------------------------------------------------------------
png("./coveffectsplot2.png",width =9 ,height = 6,units = "in",res=72)
 plotlist<- coveffectsplot::forest_plot(plotdata,
            ref_area = c(0.8, 1/0.8),
            x_facet_text_size = 13,
            y_facet_text_size = 13,
            interval_legend_text = "Median (points)\n90% CI (horizontal lines)",
            ref_legend_text = "Reference\n(vertical line)\n+/- 20% ratios\n(gray area)",
            area_legend_text = "Reference\n(vertical line)\n+/- 20% ratios\n(gray area)",
            xlabel = "Fold Change Relative to Parameter",
            facet_formula = "covname~.",
            facet_switch = "both",
            facet_scales = "free",
            facet_space = "fixed",
            paramname_shape = FALSE,
            table_position = "right",
            table_text_size = 4,
            plot_table_ratio = 4,
            show_table_facet_strip = "none",
            legend_space_x_mult = 0.5,
            ref_area_col = rgb( col2rgb("gray50")[1], col2rgb("gray50")[2],col2rgb("gray50")[3],
             max = 255, alpha = 0.1*255 ) ,
             interval_col = "steelblue",
            return_list = TRUE)
egg::ggarrange(
      plotlist[[1]]+ 
   labs(x= expression(paste("Changes Relative to ",
                                     CL["subscript"]^alpha["beta"], " Reference"),
                               sep=""))+
      theme(strip.text.y =  element_text(colour="blue")),
       plotlist[[2]] ,
      nrow = 1,
      widths = c(4, 1)
    )


 
dev.off()

