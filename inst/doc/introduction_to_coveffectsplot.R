## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- echo=FALSE, results='asis'-----------------------------------------
  set.seed(657687)
  df<- data.frame(
MASS::mvrnorm(n=1000 ,
                mu =c(10,0.75,1.5),
                Sigma=matrix(c(0.2,0.01,0.01,0.01,0.0225,0.01,0.01,0.01,0.0225),3,3,byrow = TRUE) 
))
names(df)<- c("POPCL","dWTdCL","dSEXdCL")
knitr::kable(head(df,5))
dfcov<- data.frame(
MASS::mvrnorm(n=1000 ,
                mu =c(68,85),
                Sigma=matrix(c(15,0.01,0.01,20),2,2,byrow = TRUE) 
))
names(dfcov)<- c("WTWOMAN","WTMAN")
dfcovlong <- tidyr::gather(dfcov)
ggplot2::ggplot(dfcovlong,ggplot2::aes(x=value,fill=key))+
  ggplot2::geom_density(,alpha=0.2)+
  ggplot2::labs(fill="",x="Weight (kg)")+
  ggplot2::theme(legend.position = c(0.65,0.95),legend.background = 
                   ggplot2::element_rect(fill="transparent"))+
  ggplot2::guides(fill=ggplot2::guide_legend(reverse = TRUE))

CLBSVdistribution <- data.frame(CL= 10*exp(rnorm(1000,0,sd=0.30)))
CLBSVdistribution$CLBSV<- CLBSVdistribution$CL/10

dfbsv<- as.data.frame(
  quantile(CLBSVdistribution$CLBSV,probs=c(0,0.05,0.25,0.5,0.75,0.95,1))
)
names(dfbsv)<- "BSVquantilevalue"
dfbsv$quantile<- rownames(dfbsv)


ggplot2::ggplot(CLBSVdistribution,ggplot2::aes(x=CLBSV))+
  ggplot2::geom_density(,alpha=0.2)+
  ggplot2::geom_vline(xintercept = c(0.8179004,1.2271218),size=3,col="blue",alpha=0.6)+
  ggplot2::geom_vline(xintercept = c(0.6073418,1.6259988),size=2,col="blue",alpha=0.3)

knitr::kable(dfbsv,row.names=FALSE)


## ----fig.width= 7--------------------------------------------------------
dfeffects<- df
dfeffects$REF <- dfeffects$POPCL/ median(dfeffects$POPCL)
dfeffects$WT_50 <- dfeffects$REF*(50/70)^dfeffects$dWTdCL
dfeffects$WT_90 <-  dfeffects$REF*(90/70)^dfeffects$dWTdCL
dfeffects$SEX_Male <- dfeffects$dSEXdCL
dfeffects$SEX_Male_WT_90 <- dfeffects$dSEXdCL*dfeffects$REF*(90/70)^dfeffects$dWTdCL

dfeffects$SEX_Male <- dfeffects$dSEXdCL

dfeffects<- dfeffects[,c("WT_50","WT_90","SEX_Male","SEX_Male_WT_90","REF")]
dfeffects$BSV<-  CLBSVdistribution$CLBSV


dflong <- tidyr::gather(dfeffects)
ggplot2::ggplot(dflong,ggplot2::aes(x=value,y=key,fill=factor(..quantile..)))+
ggridges::stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE,
                              quantile_lines = TRUE, rel_min_height = 0.01,
quantiles = c(0.025,0.5, 0.975)) +
  ggplot2::scale_fill_manual(
    name = "Probability", values = c("#FF0000A0", "white","white", "#0000FFA0"),
    labels = c("(0, 0.025]", "(0.025, 0.5]","(0.5, 0.975]", "(0.975, 1]")
  )+
      ggplot2::annotate(
        "rect",
        xmin = 0.8,
        xmax = 1.25,
        ymin = -Inf,
        ymax = Inf,
        fill = "gray",alpha=0.4
      )+
  ggplot2::geom_vline(
      ggplot2::aes(xintercept = 1),
      size = 1
    )+
  ggplot2::theme_bw()+
  ggplot2::labs(x="Effects Relative to parameter reference value",y="")


## ------------------------------------------------------------------------
dfeffects$SEX_Male_WT_90<- NULL
dfeffectslong<- tidyr::gather(dfeffects)
dfeffectslong<- dplyr::group_by(dfeffectslong,key)
dfeffectslongsummaries<- dplyr::summarise(dfeffectslong,mid=quantile(value,0.5),
                                   lower=quantile(value,0.025),
                                   upper=quantile(value,0.975))

dfeffectslongsummaries$paramname <- "CL"
dfeffectslongsummaries$covname <- c("BSV","REF","SEX","Weight","Weight")
dfeffectslongsummaries$label <- c("95% of patients","70 kg/Woman","Man","50 kg", "90 kg")
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


## ----fig.width=7---------------------------------------------------------
ggplot2::ggplot(data = plotdata, ggplot2::aes_string(
      y = "label",
      x = "mid",
      xmin = "lower",
      xmax = "upper"
    )) +
    ggstance::geom_pointrangeh(
      position = ggstance::position_dodgev(height = 0.75),
      ggplot2::aes(color = "95 %CI\nCovariate Effects"),
      size = 1,
      alpha = 1
    )+
  ggplot2::facet_grid(covname~.,scales="free_y",switch="y")+
  ggplot2::labs(y="",x="Effects Relative to Reference Value",
                colour="")

png("coveffectsplot.png",width =9 ,height = 6,units = "in",res=72)
coveffectsplot::forest_plot(plotdata,
            ref_area = c(0.8, 1/0.8),
            x_facet_text_size = 13,
            y_facet_text_size = 13,
            ref_legend_text = "Reference (vertical line)\n+/- 20% ratios (gray area)",
            area_legend_text = "Reference (vertical line)\n+/- 20% ratios (gray area)",
            xlabel = "Fold Change Relative to Parameter",
            facet_formula = "covname~.",
            facet_switch = "both",
            facet_scales = "free",
            facet_space = "fixed",
            paramname_shape = TRUE,
            table_position = "right",
            table_text_size=4,
            plot_table_ratio = 4)
  dev.off()
  

