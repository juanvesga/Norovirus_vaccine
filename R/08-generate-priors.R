# -------------------------------------------------------------------------
# input -------------------------------------------------------------------
# -------------------------------------------------------------------------
root    <- here::here()
outfile <- file.path(root, "output", "priors.qs2")

# Packages
source(file.path(root, "R", "modify_attach.R"))
modify_attach(qs2, include.only = c("qs_save","qs_read"))
modify_attach(monty, include.only = c("monty_dsl"))

# -------------------------------------------------------------------------
# Load data ---------------------------------------------------------------
# -------------------------------------------------------------------------

# mean
last_best=
  c(
    beta_1      =    0.13410824,
    beta_2      =    0.12708480,
    beta_3      =    0.20300807,
    beta_4      =    0.27545582,
    aduRR       =    0.58903565,
    maternalAB  =  367.16148967,
    imm_yr      =     16.15967444 ,
    imm_fac     =     1.70697752,
    repfac_0    =     64.26091774, 
    repfac_5    =     450.23810289,
    repfac_15   =     345.14095794,  
    repfac_65p  =    59.89184305,
    crossp_GI   =        0.09838119,
    crossp_GII  =      0.01773758
  )

prior <- monty_dsl({
  beta_1 ~ Gamma(shape = 2, scale = 0.15/2)
  beta_2 ~ Gamma(shape = 2, scale = 0.15/2)
  beta_3 ~ Gamma(shape = 3, scale = 0.2/3)
  beta_4 ~ Gamma(shape = 3, scale = 0.2/3)
  aduRR  ~ Uniform(0,1)
  maternalAB ~ Gamma(shape = 3, scale = 150/3)
  imm_yr~ Gamma(shape = 4, scale = 20/4)
  imm_fac~ Uniform(0,5)
  repfac_0 ~ Gamma(shape = 6, scale = 128/6)
  repfac_5 ~ Gamma(shape = 12, scale = 523/12)
  repfac_15~ Gamma(shape = 12, scale = 500/12)
  repfac_65p ~ Gamma(shape = 6, scale = 98/6)
  crossp_GI ~ Uniform(0,1)
  crossp_GII~ Uniform(0,1)
})







# -------------------------------------------------------------------------
# Save data ---------------------------------------------------------------
# -------------------------------------------------------------------------
qs_save(prior, outfile)
