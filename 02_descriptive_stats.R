#### Libs and data ########################################################################################################
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(devtools)
library(grid)
library(gridExtra)
library(dplyr)
library(tidyr)
library(arsenal)

load("data.Rdata")

#### Tables 4 latex ###############################################################################################

pars = c("node_status","grade", "stage", "tumor_size","survival_time", "event","age","chemo", "hormonal","radiation","ER",
         "PR","HER2","survival_time", "study_id")
data1 = rbind(
  data.ExprImage.fil[,pars],
  data.MB.clin[,pars],
  data.X16.clin[,pars])
data2 = rbind(
  data.X1.clin[,pars],
  data.X2.clin[,pars],
  data.X3.clin[,pars],
  data.X4.clin[,pars])
data3 = rbind(
  data.X5.clin[,pars],
  data.X6A.clin[,pars],
  data.X7.clin[,pars],
  data.X8.clin[,pars])
data4 = rbind(
  data.X9A.clin[,pars],
  data.X10.clin[,pars],
  data.X11.clin[,pars],
  data.X12.clin[,pars])
data5 = rbind(
  data.X13.clin[,pars],
  data.X14.clin[,pars],
  data.X15.clin[,pars])

data3$radiation = "NA"
data4$HER2 = "NA"

mycontrols  <- tableby.control(test=FALSE, total=FALSE,
                               numeric.test="kwt", cat.test="chisq",
                               numeric.stats=c("median", "q1q3"),
                               cat.stats=c("countpct"),
                               stats.labels=list(median='Median', q1q3='Q1,Q3'))
D = tableby(study_id ~ age+node_status+grade+tumor_size+survival_time+event+chemo+ hormonal+radiation+
               ER+PR+HER2+survival_time, data = data3, control = mycontrols)
df = as.data.frame(summary(D, text = T))
rownames(df) = NULL
print(xtable::xtable(df),  include.rownames = FALSE)




##### Plots #########################################################################################################
X = list(     data.X1.clin,
              data.X2.clin,
              data.X3.clin,
              data.X4.clin,
 #             data.X5.clin,
              data.X6A.clin,
              data.X7.clin,
              data.X8.clin,
              data.X9A.clin,
              data.X10.clin,
              data.X11.clin,
              data.X12.clin,
              data.X13.clin,
              data.X14.clin,
              data.X15.clin)
parameters = c("node_status",
              "grade",
              "tumor_size",
              "survival_time", 
              "event",
              "age",
              "chemo", 
              "hormonal",
              "radiation",
#              "stage",
              "ER",
              "PR",
              "HER2",
              "study_id",
              "survival_time")

#subset all dataframes with parameters
subset_df = function(df){return(df[,parameters])}

Xs = lapply(X, subset_df)

#stack all dataframes in one big df
df = bind_rows(Xs)
df = subset(df, df$event!="n/a")
df$PR[df$PR=="n/a"] = NA
df$survival_time = df$survival_time/12.
dim(df)

313 + 1994 + 3273 + 2226


table(df$event)



plot_factor = function(df, parameter, ylim = 5000){
                ggplot(df, aes_string(x=parameter)) + 
                geom_bar(aes(fill = event), color="black",alpha = 0.4) + 
                scale_fill_manual(labels = c("Control", "Event"), 
                                  values = alpha(c("blue", "red", "gray"), .6)) +
                geom_text(aes(label=..count..), stat='count', position=position_dodge(0.9),vjust=-0.2) +
                theme(legend.position="none", axis.title.y=element_blank()) + 
                ylim(0, ylim)}
plot_numeric = function(df, parameter, ylim = 5000, xmin = 0, xmax = 15){
                ggplot(df, aes_string(x=parameter))+
                geom_histogram(data=subset(df,event == 0),color="blue", fill="blue", alpha = 0.2, bins = 25) +
                geom_histogram(data=subset(df,event == 1),color="red", fill="red",  alpha = 0.3, bins = 25) +
                theme(axis.title.y=element_blank())+
                ylim(0, ylim) + xlim(xmin,xmax) + xlab(parameter) }

#GSE family

#factors
study_id_plt =   plot_factor(df, parameters[13], ylim=400)
node_status_plt= plot_factor(df, parameters[1] , ylim=1500)
grade_plt =      plot_factor(df, parameters[2] , ylim=1000)
event_plt =      plot_factor(df, parameters[5] , ylim=1700)
chemo_plt =      plot_factor(df, parameters[7] , ylim=1200)
hormonal_plt =   plot_factor(df, parameters[8] , ylim=1000)
radiation_plt =  plot_factor(df, parameters[9] , ylim=1800)
er_plt =         plot_factor(df, parameters[10], ylim=1600)
pr_plt =         plot_factor(df, parameters[11], ylim=1300)
her2_plt =       plot_factor(df, parameters[12], ylim=1800)
#numericals
age_plt =        plot_numeric(df, parameters[6], ylim = 170, xmin = 20, xmax = 100)
size_plt =       plot_numeric(df, parameters[3], ylim = 350, xmax = 10)
time_plt =       plot_numeric(df, parameters[14], ylim = 200, xmax = 20)


# Move to a new page
grid.newpage()
# Create layout : nrow = 2, ncol = 2
pushViewport(viewport(layout = grid.layout(3, 5)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
# Arrange the plots
print(study_id_plt, vp =    define_region(1, 1:5))
print(event_plt, vp =       define_region(2, 1))
print(age_plt, vp =         define_region(2, 2))
print(size_plt, vp =        define_region(2, 3))
print(node_status_plt, vp = define_region(2, 4))
print(time_plt, vp =        define_region(2, 5))
print(er_plt, vp =          define_region(3, 1))
print(pr_plt, vp =          define_region(3, 2))
print(grade_plt, vp =       define_region(3, 3))
print(hormonal_plt, vp =    define_region(3, 4))
print(chemo_plt, vp =       define_region(3, 5))
#print(radiation_plt, vp = define_region(3, 2))








#swedish data

S = data.X5.clin[,parameters]
S$survival_time = S$survival_time/12.
S$tumor_size = S$tumor_size/10.

#factors
node_status_plt= plot_factor(S, parameters[1] , ylim=2100)
grade_plt =      plot_factor(S, parameters[2] , ylim=1600)
event_plt =      plot_factor(S, parameters[5] , ylim=3100)
chemo_plt =      plot_factor(S, parameters[7] , ylim=2000)
hormonal_plt =   plot_factor(S, parameters[8] , ylim=2600)
radiation_plt =  plot_factor(S, parameters[9] , ylim=1800)
er_plt =         plot_factor(S, parameters[10], ylim=3000)
pr_plt =         plot_factor(S, parameters[11], ylim=2700)
her2_plt =       plot_factor(S, parameters[12], ylim=2800)
#numericalsS
age_plt =        plot_numeric(S, parameters[6], ylim = 370, xmin = 20, xmax = 100)
size_plt =       plot_numeric(S, parameters[3], ylim = 550, xmax = 7)
time_plt =       plot_numeric(S, parameters[14], ylim = 300, xmax = 7)


# Move to a new page
grid.newpage()
# Create layout : nrow = 2, ncol = 2
pushViewport(viewport(layout = grid.layout(2, 5)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
# Arrange the plots
print(event_plt, vp =       define_region(1, 1))
print(age_plt, vp =         define_region(1, 2))
print(size_plt, vp =        define_region(1, 3))
print(node_status_plt, vp = define_region(1, 4))
print(time_plt, vp =        define_region(1, 5))
print(er_plt, vp =          define_region(2, 1))
print(pr_plt, vp =          define_region(2, 2))
print(grade_plt, vp =       define_region(2, 3))
print(hormonal_plt, vp =    define_region(2, 4))
print(chemo_plt, vp =       define_region(2, 5))
#print(radiation_plt, vp = define_region(3, 2))


#MB training and valid 
parameters = c("node_status",
               "grade",
               "tumor_size",
               "survival_time", 
               "event",
               "age",
               "chemo", 
               "hormonal",
               "radiation",
               #              "stage",
               "ER",
               "PR",
               "HER2",
               "survival_time")
M = bind_rows(data.MBtraining.clin[,parameters], data.MBvalid.clin[,parameters])
M$survival_time = M$survival_time/12.
dim(M) #[1] 1483   13


#factors
node_status_plt= plot_factor(M, parameters[1] , ylim=800)
grade_plt =      plot_factor(M, parameters[2] , ylim=800)
event_plt =      plot_factor(M, parameters[5] , ylim=900)
chemo_plt =      plot_factor(M, parameters[7] , ylim=1200)
hormonal_plt =   plot_factor(M, parameters[8] , ylim=900)
radiation_plt =  plot_factor(M, parameters[9] , ylim=1100)
er_plt =         plot_factor(M, parameters[10], ylim=1200)
pr_plt =         plot_factor(M, parameters[11], ylim=800)
her2_plt =       plot_factor(M, parameters[12], ylim=1300)
#numericalsS
age_plt =        plot_numeric(M, parameters[6], ylim = 100, xmin = 20, xmax = 100)
size_plt =       plot_numeric(M, parameters[3], ylim = 150, xmax = 7)
time_plt =       plot_numeric(M, parameters[13], ylim = 120, xmax = 30)


# Move to a new page
grid.newpage()
# Create layout : nrow = 2, ncol = 2
pushViewport(viewport(layout = grid.layout(2, 5)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
# Arrange the plots
print(event_plt, vp =       define_region(1, 1))
print(age_plt, vp =         define_region(1, 2))
print(size_plt, vp =        define_region(1, 3))
print(node_status_plt, vp = define_region(1, 4))
print(time_plt, vp =        define_region(1, 5))
print(er_plt, vp =          define_region(2, 1))
print(pr_plt, vp =          define_region(2, 2))
print(grade_plt, vp =       define_region(2, 3))
print(hormonal_plt, vp =    define_region(2, 4))
print(chemo_plt, vp =       define_region(2, 5))
#print(radiation_plt, vp = define_region(3, 2))





# ExprImage 
parameters = c("node_status",
               "grade",
               "tumor_size",
               "survival_time", 
               "event",
               "age",
               "chemo", 
               "hormonal",
               "radiation",
               #              "stage",
               "ER",
               "PR",
               "HER2",
               "survival_time")
E = bind_rows(data.ExprImage.raw[,parameters])
E$survival_time = E$survival_time/12.
dim(M) #[1] 1483   13


#factors
node_status_plt= plot_factor(E, parameters[1] , ylim=200)
grade_plt =      plot_factor(E, parameters[2] , ylim=200)
event_plt =      plot_factor(E, parameters[5] , ylim=200)
chemo_plt =      plot_factor(E, parameters[7] , ylim=150)
hormonal_plt =   plot_factor(E, parameters[8] , ylim=200)
radiation_plt =  plot_factor(E, parameters[9] , ylim=200)
er_plt =         plot_factor(E, parameters[10], ylim=300)
pr_plt =         plot_factor(E, parameters[11], ylim=300)
her2_plt =       plot_factor(E, parameters[12], ylim=300)
#numericalsS
age_plt =        plot_numeric(E, parameters[6], ylim = 30, xmin = 20, xmax = 100)
size_plt =       plot_numeric(E, parameters[3], ylim = 30, xmax = 7)
time_plt =       plot_numeric(E, parameters[13], ylim = 23, xmax = 16)


# Move to a new page
grid.newpage()
# Create layout : nrow = 2, ncol = 2
pushViewport(viewport(layout = grid.layout(2, 5)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
# Arrange the plots
print(event_plt, vp =       define_region(1, 1))
print(age_plt, vp =         define_region(1, 2))
print(size_plt, vp =        define_region(1, 3))
print(node_status_plt, vp = define_region(1, 4))
print(time_plt, vp =        define_region(1, 5))
print(er_plt, vp =          define_region(2, 1))
print(pr_plt, vp =          define_region(2, 2))
print(grade_plt, vp =       define_region(2, 3))
print(hormonal_plt, vp =    define_region(2, 4))
print(chemo_plt, vp =       define_region(2, 5))
#print(radiation_plt, vp = define_region(3, 2))

997*2

