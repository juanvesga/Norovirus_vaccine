# -------------------------------------------------------------------------
# input -------------------------------------------------------------------
# -------------------------------------------------------------------------
root           <- here::here()
infile         <-file.path(root, "output", "scenario_results.qs2")
outfile        <- file.path(root, "output", "scenario_results.qs2")

#Functions
source(file.path(root, "R", "modify_attach.R"))
source(file.path(root, "R", "utility_functions.R"))

# Packages
modify_attach(qs2, include.only = c("qs_save","qs_read"))

library(ggplot2)
library(dplyr)
library(data.table)


runs <- qs_read(infile)

my_names = c("days", "scenario", "cases_averted")

shortlist = lapply(runs, "[", , my_names)
shortDF <- rbindlist(shortlist)

get_scenario_quantiles<- function(dframe,varname){

  scen_names<-unique(dframe$scenario)  
  n=length(scen_names)
  
  datalist = vector("list", length = n)
  
  for (i in 1:n) {
    # ... make some data
    dat <- dframe[dframe$scenario==scen_names[i],] %>% group_by(days) %>% 
      do(data.frame(t(quantile(.[[varname]], probs = c(0.025,0.5,0.975)))))

    dat$scenario <- scen_names[i]  # maybe you want to keep track of which iteration produced it?
    
    datalist[[i]] <- dat # add it to your list
  }
  
  big_data = do.call(rbind, datalist)
  
}

cases_averted<-get_scenario_quantiles(shortDF,"cases_averted")


p <- ggplot(cases_averted) + geom_line(aes(y=X50., x=days, colour = scenario))+
  geom_ribbon(aes(ymin=X2.5., ymax=X97.5., x=days, fill = scenario), alpha = 0.3)





