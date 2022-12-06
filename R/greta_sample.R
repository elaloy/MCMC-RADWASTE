#### Load libraries ####

library(tidyverse)
library(here)
library(extraDistr) # rdirichlet
library(abind)
library(matlab)
library(rlist) # list.rbind
library(tictoc)
library(greta)
library(greta.gp) # for the GPs
library(bayesplot)
library(coda)

#### Set global inversion parameters ####
DoSynthetic<-FALSE
SaveModelplot<-FALSE
DoSampling<-TRUE
DoSaveOutput<-TRUE

used_meas<-c("PNCC","3AX")
nseg<-20 # number of segments
ndet<-20 # number of SGS-3AX detector positions
n_pu<-5 # Dirichlet prior for Pu-238, Pu-239, Pu-240, Pu-241 and Pu-242 (5)

sel_nucl<-c("Am-241","Cm-244","Np-237","Pu-238","Pu-239","Pu-240","Pu-241","Pu-242","U-235")
n_nucl<-length(sel_nucl)
sel_nucl_gamma<-c("Am-241","Pu-238","Pu-239","Pu-240","Pu-241")
idx_gamma<-which(sel_nucl %in% sel_nucl_gamma)

if(DoSynthetic){

  DoRefGP<-TRUE  # TRUE: True Pu profile is GP-like, FALSE: true Pu profile is "flat and spiky"

  sampling_case <- paste0("synth_","_refGP_",DoRefGP,"_",paste(used_meas,collapse="_"))
}else{
  sampling_case <- paste0("real_",paste(used_meas,collapse="_"))
}

#### HMC parameters ####
warmup<-10000
n_samples<-5000
n_chains<-8
seed<-1978
set.seed(seed)

#### Set priors ####
prior_for_vec<-"multilevel"
# "flat_dirichlet", "fram_based", "genie2k_based", "multilevel"

prior_for_pu<-"gp"
#"loguniform" or "gp" (log(pu) ~gp) or "gaussian" (log(pu)~N(log(1),log(exp(1)))

LinkCm2Pu<-"2Pu240"

if (LinkCm2Pu %in% c("None","2totPu")){
  prior_for_cm_st<-"loguniform"
  }else{prior_for_cm_st<-"gamma"}


Use1vec<-FALSE # either infer a single vector for all segments (TRUE) or one vector per segment (FALSE)

prior_for_am2pu<-"uniform"
prior_for_u2pu<-"uniform"
prior_for_bckg<-"norm2pois"
prior_for_lam<-"uniform"

if (Use1vec==TRUE){
  vectype<-"single_vec"
}else{
  vectype<-"multi_vec"
}

#### Load data for 3AX-SGS ####

eh_3ax<-read_rds(here("input/sgs/eh_3ax.rds"))
ep_3ax<-read_rds(here("input/sgs/ep_3ax.rds"))
# The eh and ep are of dimensions:
# source location (1 = bottom, 20 = top) * energy (sorted by energy value) * detector location (1=bottom, 20 = top)


#Specific activity [Bq/g]
spec_a<-read_rds(here("input/specific_activity.rds")) %>%
  filter(Nuclide %in% sel_nucl) %>%
  arrange(Nuclide) %>%
  select(Value) %>%
  pull()

spa<-array(rep(spec_a,each=nseg),c(nseg,n_nucl))

#### Load and preprocess count data for 3AX ####
gross_time_3ax<-300 #[s]
bckg_time_3ax<-300
tb_3ax<-bckg_time_3ax


# Get indices of the energy peaks (sorted by ascending values)
# to be used per detector's location (valid measurements) in eh_3ax and ep_3ax
ide3_s<-readRDS(here("input/sgs/idx_valid_meas_per_segment.rds"))

# Get indices of energy peaks (sorted by ascending values), sorted by nuclide
ide3<-readRDS(here("input/sgs/idx_sorted_energies_by_nuclide.rds"))

# Get initial number of peaks per nuclide, sorted by nuclide
idx3<-readRDS(here("input/sgs/peaks_per_nuclide.rds"))

# Get emission probabilities, sorted by nuclide
I3<-readRDS(here("input/sgs/sorted_b_by_nuclide.rds"))
nene<-length(I3) # we have 13 energies but peak number 11 (Pu-239 645.94 keV) is never considered
# because there are no valid measurements for this peak - see index vectors in the ide3_s list

# Get vector of observed net and background counts
# 12 valid peaks * 20 segments = 240 data points, sorted by segment
obs_net <- readRDS(here("input/sgs/NetCounts_3ax.rds"))
obs_bckg <- readRDS(here("input/sgs/BckgCounts_3ax.rds"))


#### Load and process data for PNCC ####
obs_Rexp_m<-read_rds(here("input/pncc/obs_Rexp_m.rds"))
alfa_rexp<-0.05
sigma_Rexp_m<-alfa_rexp*obs_Rexp_m

Qref<-read_rds(here("input/pncc/Qref.rds"))
QiQref<-read_rds(here("input/pncc/QiQref.rds"))
Fs_prime_240Pu<-read_rds(here("input/pncc/Fs_prime_240Pu.rds"))

eh_pncc<-read_rds(here("input/pncc/eh_pncc.rds"))
ep_pncc<-read_rds(here("input/pncc/ep_pncc.rds"))

# eh_pncc and ep_pncc are the Rsim_m variables in the ANIMA paper



qiqr<-QiQref %>%
  arrange(radionuclide) %>%
  filter(radionuclide %in% sel_nucl) %>%
  select(`Q/Qref`) %>%
  replace_na(list(`Q/Qref`=0)) %>%
  as_vector() %>%
  unname()
qiqr<-array(rep(qiqr,each=nseg),c(nseg,n_nucl))

#### Create reference values that are replacing the observations in case of a synthetic experiment ####
if(DoSynthetic){

  ref_bckg<-obs_bckg

  #### Create reference vector(s) - ref_vec [m% Pu]

  if (Use1vec==TRUE){ # a single vector for all segments - not used herein

    ref_vec<-c(0.01,0.73,0.23,0.01,0.02)


    ref_vec<-array(rep(ref_vec,each=nseg),c(nseg,n_pu))

  }else{ # one vector per segment

    if(prior_for_vec=="flat_dirichlet"){
      mu_alph<-c(0.6408,	77.6854,	18.5527,	1.3227,	1.7974) # FRAM/ISOCS results * 100
      ref_vec = rdirichlet(n=20, mu_alph)

    }else if(prior_for_vec=="multilevel"){
      mu_alph<-c(0.6408,	77.6854,	18.5527,	1.3227,	1.7974) # FRAM/ISOCS results * 100
      ref_vec = rdirichlet(n=20, mu_alph)
    }
  }
  # Pu mass [g]
  if(DoRefGP==FALSE){
    low_pu<-0.01
    large_pu<-2.5
    ref_totPu<-numeric(nseg)+low_pu
    ref_totPu[5]<-large_pu
    ref_totPu[16]<-large_pu
  }else{
    # uses a GP-like curve for the true totPu
    low_pu<-0.15
    large_pu<-2.5
    # simulate data
    x <- seq(nseg) %>% as.numeric()
    y <- sin(0.5*x)
    ref_totPu<- low_pu + ((y+1)*0.5)*(large_pu-low_pu)

    DofitGP2RefPu<-FALSE

    if(DofitGP2RefPu){
      # now fit a model to ref_totPu, using the same GP model as for the real data
      gp_len0 <- uniform(0.5, 20) # prior for lengthscale
      gp_var0 <- lognormal(0,1) # prior for variance
      gp_sd0<-sqrt(gp_var0)
      gp_mean0 <- normal(log(1),log(exp(1))) # prior for mean

      obs_sd0 <- lognormal(0, 1) # prior for observation uncertainty

      # kernel & GP
      kernel0 <- rbf(gp_len0, gp_var0)
      f00 <- gp(x, kernel0)

      f0 <- f00 + gp_mean0

      # output data
      y0 <-ref_totPu

      # likelihood
      distribution(y0) <- normal(exp(f0), obs_sd0)

      # prediction
      x0_plot<-x
      f0_plot <- project(f00, x0_plot)

      # fit the model by HMC
      m0 <- model(f0_plot,gp_var0,gp_len0,obs_sd0,gp_mean0)
      z0 <- greta::mcmc(m0,one_by_one = TRUE,warmup=2000,n_samples=5000)

      ref_z0<-colMeans(list.rbind(z0))

      # quick and dirty plot in base R
      plot(y0 ~ x, pch = 16, col = grey(0.4), xlim = c(0, 21), ylim = c(0, 3))
      for (i in 1:200) {
        lines(exp(z0[[1]][i,1:20] + z0[[1]][i,24]) ~ x0_plot, #
              lwd = 2,
              col = rgb(0.7, 0.1, 0.4, 0.1))

      }
    }
  }
  # Am-241
  # Am-241 to totPu ratio
    #FRAM/ISOCS gives 0.0423
    lb_am2pu <- 0.023
    ub_am2pu <-0.063
  if(Use1vec){
    # Am-241/tot-Pu ratio
    ref_am2pu<-rep(0.043,each=nseg)
  }else{
    # Am-241/tot-Pu ratio
    ref_am2pu <- runif(nseg)*(ub_am2pu-lb_am2pu) + lb_am2pu
  }

  # U-235
  lb_u2pu <- 0.05
  ub_u2pu <-0.1
  #FRAM/ISOCS gives 0.0736
  if(Use1vec){
    ref_u2pu<-rep(0.075,each=nseg)
  }else{
    ref_u2pu <- runif(nseg)*(ub_u2pu-lb_u2pu) + lb_u2pu
  }
  # Np-237
  lb_np2pu <- 0.0001
  ub_np2pu <-0.01
  #FRAM/ISOCS gives 0.001682
  if(Use1vec){
    ref_np2pu<-rep(0.002,each=nseg)
  }else{
    ref_np2pu <- runif(nseg)*(ub_np2pu-lb_np2pu) + lb_np2pu
  }


  if(LinkCm2Pu=="2Pu240"){
    ref_pu239<-ref_vec[,2]*ref_totPu
    ref_pu240<-ref_vec[,3]*ref_totPu
    #ref_cm_st=log10((Pu239+Pu240)/Cm244)
    ref_cm_st<-5
    ref_cm244<-(ref_pu240+ref_pu239)/(10^ref_cm_st)
  }

  #### Create reference lambdas ###
  ref_lam<-numeric(20) +1
  ref_lam[16]<-0

  #### Set ref parameters ###
  ref_vec0<-ref_vec
  ref_vec<-as.vector(ref_vec) # ref_vec will come as a vector in the DAG
  ref_par<-list()

  ref_par[[1]]<-ref_vec0
  ref_par[[2]]<-ref_totPu
  ref_par[[3]]<-ref_am2pu
  ref_par[[4]]<-ref_lam
  ref_par[[5]]<-ref_bckg
  ref_par[[6]]<-ref_cm_st
  ref_par[[7]]<-ref_u2pu
  ref_par[[8]]<-ref_np2pu

  ref_am241<-ref_am2pu*ref_totPu

  ref_u235<-ref_u2pu*ref_totPu
  ref_np237<-ref_np2pu*ref_totPu

  ref_mass<-abind(ref_am241,ref_cm244,ref_np237,array(ref_vec,c(nseg,n_pu))*array(rep(ref_totPu,times=n_pu),dim = c(nseg,n_pu)),ref_u235)

  #### activity [Bq] for gamma spec ###
  # only consider the radionuclides detected by gamma spec:
  # spa<-array(rep(spec_a,each=nseg),c(nseg,n_nucl))
  ref_a2<-spa[,idx_gamma]*ref_mass[,idx_gamma]
  ref_a<-array(t(ref_a2)) # flatten the 2D array of activities with C-ordering / row-major order

  #### Create reference net counts - 3AX ###

  # ref_lam is indexed w.r.t. source location and NOT detector location
  # ref_a2 is indexed w.r.t. source location and NOT detector location
  ref_lam3<-array(rep(ref_lam,times=nene*ndet),dim = c(nseg, nene, ndet))
  # ref_lam3 is source's location x energy x detector's location

  ref_ees3 <- ref_lam3*eh_3ax+(1-ref_lam3)*ep_3ax

  ref_am3<-abind(rep(ref_a[1:5],idx3)*I3,rep(ref_a[6:10],idx3)*I3,rep(ref_a[11:15],idx3)*I3,
                 rep(ref_a[16:20],idx3)*I3,rep(ref_a[21:25],idx3)*I3,rep(ref_a[26:30],idx3)*I3,
                 rep(ref_a[31:35],idx3)*I3,rep(ref_a[36:40],idx3)*I3,rep(ref_a[41:45],idx3)*I3,
                 rep(ref_a[46:50],idx3)*I3,rep(ref_a[51:55],idx3)*I3,rep(ref_a[56:60],idx3)*I3,
                 rep(ref_a[61:65],idx3)*I3,rep(ref_a[66:70],idx3)*I3,rep(ref_a[71:75],idx3)*I3,
                 rep(ref_a[76:80],idx3)*I3,rep(ref_a[81:85],idx3)*I3,rep(ref_a[86:90],idx3)*I3,
                 rep(ref_a[91:95],idx3)*I3,rep(ref_a[96:100],idx3)*I3,along=0)

  # cps1 is for detector's location 1, cps2 is for detector's location 2, ...

  cps1<-colSums(ref_am3*ref_ees3[,,1][,ide3])[ide3_s[[1]]]
  cps2<-colSums(ref_am3*ref_ees3[,,2][,ide3])[ide3_s[[2]]]
  cps3<-colSums(ref_am3*ref_ees3[,,3][,ide3])[ide3_s[[3]]]
  cps4<-colSums(ref_am3*ref_ees3[,,4][,ide3])[ide3_s[[4]]]
  cps5<-colSums(ref_am3*ref_ees3[,,5][,ide3])[ide3_s[[5]]]
  cps6<-colSums(ref_am3*ref_ees3[,,6][,ide3])[ide3_s[[6]]]
  cps7<-colSums(ref_am3*ref_ees3[,,7][,ide3])[ide3_s[[7]]]
  cps8<-colSums(ref_am3*ref_ees3[,,8][,ide3])[ide3_s[[8]]]
  cps9<-colSums(ref_am3*ref_ees3[,,9][,ide3])[ide3_s[[9]]]
  cps10<-colSums(ref_am3*ref_ees3[,,10][,ide3])[ide3_s[[10]]]
  cps11<-colSums(ref_am3*ref_ees3[,,11][,ide3])[ide3_s[[11]]]
  cps12<-colSums(ref_am3*ref_ees3[,,12][,ide3])[ide3_s[[12]]]
  cps13<-colSums(ref_am3*ref_ees3[,,13][,ide3])[ide3_s[[13]]]
  cps14<-colSums(ref_am3*ref_ees3[,,14][,ide3])[ide3_s[[14]]]
  cps15<-colSums(ref_am3*ref_ees3[,,15][,ide3])[ide3_s[[15]]]
  cps16<-colSums(ref_am3*ref_ees3[,,16][,ide3])[ide3_s[[16]]]
  cps17<-colSums(ref_am3*ref_ees3[,,17][,ide3])[ide3_s[[17]]]
  cps18<-colSums(ref_am3*ref_ees3[,,18][,ide3])[ide3_s[[18]]]
  cps19<-colSums(ref_am3*ref_ees3[,,19][,ide3])[ide3_s[[19]]]
  cps20<-colSums(ref_am3*ref_ees3[,,20][,ide3])[ide3_s[[20]]]

  ref_cps<-abind(cps1,cps2,cps3,cps4,cps5,cps6,cps7,cps8,cps9,
                 cps10,cps11,cps12,cps13,cps14,cps15,cps16,
                 cps17,cps18,cps19,cps20,along=1)

  ref_net<-ref_cps*tb_3ax

  # replace observed net counts by the true ones and then contaminate them with Poisson noise
  obs_net<-numeric(length(ref_net))
  for (i in 1:length(ref_net)){obs_net[i]<-rpois(1,ref_net[i])}





  #### Create reference Rexp,m - PNCC ###
  # ref_mass is ordered as sel_nucl, that is,
  # [1] "Am-241" "Cm-244" "Pu-238" "Pu-239" "Pu-240" "Pu-241"
  # [7] "Pu-242" [8] "U-235"

  ref_comp_pncc<-ref_mass/array(rep(rowSums(ref_mass),times=n_nucl),c(nseg,n_nucl))

  ref_w_sum_sf_pair_rate<-rowSums(ref_comp_pncc*qiqr)
  ref_m240Pu_eq<-rowSums(ref_mass)*ref_w_sum_sf_pair_rate

  ref_eff_pncc <- ref_lam*eh_pncc+(1-ref_lam)*ep_pncc

  urel_Rsim_m<-0.5/100 # 4.4e-3

  sigma_Rsim_m<-ref_eff_pncc*urel_Rsim_m
  ref_Rexp_m <- Fs_prime_240Pu*base::sum(ref_m240Pu_eq*ref_eff_pncc)
  #ref_Rexp_m is simulated number of reals per second

  #ref_Rexp_m <- Fs_prime_240Pu*base::sum(c*(ref_eff_pn*(rnorm(nseg)*sigma_Rsim_m+1)))

  sigma_ref_Rexp_m<-alfa_rexp*ref_Rexp_m
  sigma_Rexp_m <-sigma_ref_Rexp_m

  obs_Rexp_m<-ref_Rexp_m+rnorm(1)*sigma_ref_Rexp_m
  t_pncc<-550*60 #s
}
#### Convert necessary input data to Greta arrays ####

eh3_g <- eh_3ax %>% as_data()
ep3_g <- ep_3ax %>% as_data()
I3_g <- I3 %>% as_data()
spa_arr_g<-greta_array(rep(spec_a,each=nseg),c(nseg,n_nucl))

ehpn_g <- eh_pncc %>% as_data()
eppn_g  <- ep_pncc %>% as_data()
qiqr_g <- qiqr %>% as_data()
Fs_prime_240Pu_g <- Fs_prime_240Pu %>% as_data()

#### Priors and transformations ####
#### Pu isotopic vector #
dim_puv<-n_pu
if (prior_for_vec=="flat_dirichlet"){
  alph<-numeric(n_pu)+1
  alph<-t(alph)
}else if(prior_for_vec=="fram_based"){ # # FRAM/ISOCS results
  alph<-c(0.6408,	77.6854,	18.5527,	1.3227,	1.7974)
  alph<-t(alph)

}else if(prior_for_vec=="multilevel"){

  mu_alph<-c(0.6408,	77.6854,	18.5527,	1.3227,	1.7974) # FRAM/ISOCS results
  cv<-0.5 #sigma/mu # cv is coefficient of variation, set to 0.5 here
  sd_alph<-mu_alph*cv
  alph_st<-normal(0,1,dim=length(mu_alph),truncation=c(-1/cv,Inf))
  alph <-t(alph_st*sd_alph+mu_alph)


}
if(Use1vec==TRUE){
  puv<-dirichlet(alph,dim=n_pu)

  puv_arr<-greta_array(rep(puv,each=nseg),dim=c(nseg,n_pu))

}else{
  if(prior_for_vec=="flat_dirichlet" | prior_for_vec=="fram_based" | prior_for_vec=="multilevel"){

    dim_full_puv<-n_pu*nseg
    puv1<-dirichlet(alph,dim=n_pu)

    puv2<-dirichlet(alph,dim=n_pu)

    puv3<-dirichlet(alph,dim=n_pu)

    puv4<-dirichlet(alph,dim=n_pu)

    puv5<-dirichlet(alph,dim=n_pu)

    puv6<-dirichlet(alph,dim=n_pu)

    puv7<-dirichlet(alph,dim=n_pu)

    puv8<-dirichlet(alph,dim=n_pu)

    puv9<-dirichlet(alph,dim=n_pu)

    puv10<-dirichlet(alph,dim=n_pu)

    puv11<-dirichlet(alph,dim=n_pu)

    puv12<-dirichlet(alph,dim=n_pu)

    puv13<-dirichlet(alph,dim=n_pu)

    puv14<-dirichlet(alph,dim=n_pu)

    puv15<-dirichlet(alph,dim=n_pu)

    puv16<-dirichlet(alph,dim=n_pu)

    puv17<-dirichlet(alph,dim=n_pu)

    puv18<-dirichlet(alph,dim=n_pu)

    puv19<-dirichlet(alph,dim=n_pu)

    puv20<-dirichlet(alph,dim=n_pu)

  }

  puv_arr<-abind(puv1,puv2,puv3,puv4,puv5,puv6,puv7,puv8,puv9,puv10,puv11,
                 puv12,puv13,puv14,puv15,puv16,puv17,puv18,puv19,puv20,along=1) # correct with along=1?
}


#### Pu mass #
dim_pm<-nseg
if(prior_for_pu=="uniform" | prior_for_pu=="loguniform"){
  lower_pm<-numeric(dim_pm)+1e-2
  upper_pm<-numeric(dim_pm)+100
  lower_logpm<-log(lower_pm)
  upper_logpm<-log(upper_pm)

  lb_pm<-lower_pm %>% as_data()
  ub_pm<-upper_pm %>% as_data()
  lb_logpm<-lower_logpm %>% as_data()
  ub_logpm<-upper_logpm %>% as_data()

  pm_st<-uniform(0,1,dim=dim_pm)

  if(prior_for_pu=="uniform"){
    pm<-lb_pm + pm_st*(ub_pm - lb_pm)
  }else if(prior_for_pu=='loguniform'){
    pm<-lb_logpm + pm_st*(ub_logpm - lb_logpm)
    pm <- exp(pm)
  }
}else if(prior_for_pu=="gp"){

  lengthscale <- uniform(0.5, 20)
  gp_var <- lognormal(0,1)
  gp_sd<-sqrt(gp_var)
  gp_mean <- normal(log(1),log(exp(1)))
  pm_st <- gp(1:20, rbf(lengthscale, gp_var)) + gp_mean
  pm <- exp(pm_st)

}else if(prior_for_pu=="gaussian"){

  pm_st<-normal(0,1,dim=nseg)
  mu_pm<-log(1) %>% as_data()
  sd_pm<-log(exp(1)) %>% as_data()
  pm<-pm_st*sd_pm + mu_pm
  pm<-exp(pm)

}
#### Am-241 mass to Pu mass ratio #
if(Use1vec==TRUE){
  dim_am241r<-1
}else{dim_am241r<-nseg}
if(prior_for_am2pu=="uniform" | prior_for_am2pu=="loguniform"){

  lower_am241r <- numeric(dim_am241r)+0.01
  upper_am241r <-numeric(dim_am241r)+0.15


  lower_logam241r<-log(lower_am241r)
  upper_logam241r<-log(upper_am241r)
  lb_am241r<-lower_am241r %>% as_data()
  ub_am241r<-upper_am241r %>% as_data()
  lb_logam241r<-lower_logam241r %>% as_data()
  ub_logam241r<-upper_logam241r %>% as_data()

  am241r_st<-uniform(0,1,dim=dim_am241r)

  if(prior_for_am2pu=="uniform"){
    am241r<-lb_am241r + am241r_st*(ub_am241r - lb_am241r)
  }else if(prior_for_am2pu=='loguniform'){
    am241r<-lb_logam241r + am241r_st*(ub_logam241r - lb_logam241r)
    am241r <- exp(am241r)
  }

}else if(prior_for_am2pu=="normal"){

  am241r_st<-normal(0,1,dim=nseg,truncation = c(-0.043/0.5, Inf))
  mu_am241r<-c(0.043) %>% as_data()
  sd_am241r<-c(0.5) %>% as_data()
  am241r<-am241r_st*sd_am241r + mu_am241r

}
### U-235 mass to Pu mass ratio #
if(Use1vec==TRUE){dim_u235r<-1}else{dim_u235r<-nseg}
if(prior_for_u2pu=="uniform"){
  lower_u235r <- numeric(dim_u235r)+0.01#c(0.075-0.05)
  upper_u235r <-numeric(dim_u235r)+0.15#c(0.075+0.05)
  lb_u235r<-lower_u235r %>% as_data()
  ub_u235r<-upper_u235r %>% as_data()
  u235r_st<-uniform(0,1,dim=dim_u235r)
  u235r<-lb_u235r + u235r_st*(ub_u235r - lb_u235r)
}else if(prior_for_u2pu=="normal"){
  u235r_st<-normal(0,1,dim=nseg,truncation = c(-0.075/0.5, Inf))
  mu_u235r<-c(0.075) %>% as_data()
  sd_u235r<-c(0.5) %>% as_data()
  u235r<-u235r_st*sd_u235r + mu_u235r
}

### Np-237 #
if(Use1vec==TRUE){dim_np237r<-1}else{dim_np237r<-nseg}
if(prior_for_u2pu=="uniform"){
  lower_np237r <- numeric(dim_np237r)+1e-4
  upper_np237r <-numeric(dim_np237r)+1e-2
  lb_np237r<-lower_np237r %>% as_data()
  ub_np237r<-upper_np237r %>% as_data()
  np237r_st<-uniform(0,1,dim=dim_np237r)
  np237r<-lb_np237r + np237r_st*(ub_np237r - lb_np237r)
}else if(prior_for_u2pu=="normal"){
  np237r_st<-normal(0,1,dim=nseg,truncation = c(-0.0017/0.0006, Inf))
  mu_np237r<-c(0.0017) %>% as_data()
  sd_np237r<-c(0.0006) %>% as_data()
  np237r<-np237r_st*sd_np237r + mu_np237r
}
#### Lambda vector #
dim_lam<-nseg
if(prior_for_lam=="uniform"){
  lam<- uniform(0,1,dim=dim_lam)
}
#### Bckg continuum counts #
# Order of bckg is (1) Spectrum, (2) Nuclide, (3) Energy

dim_b<-length(obs_bckg)
if(prior_for_bckg=="norm2pois"){
  b_st<-normal(0,1,dim=dim_b)
  beta<-1
  mu_b<-c(obs_bckg) %>% as_data()
  sd_b<-c(sqrt(obs_bckg)*beta) %>% as_data()

  sim_bckg<-b_st*sd_b + mu_b
}

#### Cm-244 mass #

if(LinkCm2Pu=="2Pu240"){ # gamma prior for cm_st = log10((Pu239+Pu240)/Cm244)
  #rgamma(1,shape=7.884140,rate=2.513273) for log10((Pu239+Pu240)/Cm244)
  dim_cmr<-1
  cmr_st<-gamma(shape=7.884140,rate=2.513273,dim=dim_cmr)
  if(Use1vec){
    cm244<-((puv[,2]*pm)+(puv[,3]*pm))/(10^cmr_st)
  }else{
    cm244<-((puv_arr[,2]*pm)+(puv_arr[,3]*pm))/(10^cmr_st)
  }

}
dim_cmr_st<-dim_cmr
dim_cm244<-nseg # always equal to nseg

#### Models ####

# mass of each isotope [g]
# mass vector for Pu-related isotopes and Cm-244, in alphabetical order

am241<-am241r*pm
u235<-u235r*pm
np237<-np237r*pm


mass<-abind(am241,cm244,np237,puv_arr*greta_array(rep(pm,times=n_pu),dim = c(nseg,n_pu)),u235)

if ("3AX"%in% used_meas){
  # activity [Bq] for radionuclides related to gamma spec
  a2<-spa_arr_g[,idx_gamma]*mass[,idx_gamma]
  a<-greta_array(t(a2)) # flatten the 2D array of activities with C-ordering / row-major order

  lam3<-greta_array(rep(lam,times=nene*ndet),dim = c(nseg, nene, ndet))
  # lam3_ref is source's location x energy x detector's location

  ees3 <- lam3*eh3_g+(1-lam3)*ep3_g

  am3<-abind(rep(a[1:5],idx3)*I3_g,rep(a[6:10],idx3)*I3_g,rep(a[11:15],idx3)*I3_g,
             rep(a[16:20],idx3)*I3_g,rep(a[21:25],idx3)*I3_g,rep(a[26:30],idx3)*I3_g,
             rep(a[31:35],idx3)*I3_g,rep(a[36:40],idx3)*I3_g,rep(a[41:45],idx3)*I3_g,
             rep(a[46:50],idx3)*I3_g,rep(a[51:55],idx3)*I3_g,rep(a[56:60],idx3)*I3_g,
             rep(a[61:65],idx3)*I3_g,rep(a[66:70],idx3)*I3_g,rep(a[71:75],idx3)*I3_g,
             rep(a[76:80],idx3)*I3_g,rep(a[81:85],idx3)*I3_g,rep(a[86:90],idx3)*I3_g,
             rep(a[91:95],idx3)*I3_g,rep(a[96:100],idx3)*I3_g,along=0)

  #cps1 is for detector's location 1, cps2 is for detector's location 2, ...
  #greta code
  cps1<-colSums(am3*ees3[,,1][,ide3,1])[ide3_s[[1]]]
  cps2<-colSums(am3*ees3[,,2][,ide3,1])[ide3_s[[2]]]
  cps3<-colSums(am3*ees3[,,3][,ide3,1])[ide3_s[[3]]]
  cps4<-colSums(am3*ees3[,,4][,ide3,1])[ide3_s[[4]]]
  cps5<-colSums(am3*ees3[,,5][,ide3,1])[ide3_s[[5]]]
  cps6<-colSums(am3*ees3[,,6][,ide3,1])[ide3_s[[6]]]
  cps7<-colSums(am3*ees3[,,7][,ide3,1])[ide3_s[[7]]]
  cps8<-colSums(am3*ees3[,,8][,ide3,1])[ide3_s[[8]]]
  cps9<-colSums(am3*ees3[,,9][,ide3,1])[ide3_s[[9]]]
  cps10<-colSums(am3*ees3[,,10][,ide3,1])[ide3_s[[10]]]
  cps11<-colSums(am3*ees3[,,11][,ide3,1])[ide3_s[[11]]]
  cps12<-colSums(am3*ees3[,,12][,ide3,1])[ide3_s[[12]]]
  cps13<-colSums(am3*ees3[,,13][,ide3,1])[ide3_s[[13]]]
  cps14<-colSums(am3*ees3[,,14][,ide3,1])[ide3_s[[14]]]
  cps15<-colSums(am3*ees3[,,15][,ide3,1])[ide3_s[[15]]]
  cps16<-colSums(am3*ees3[,,16][,ide3,1])[ide3_s[[16]]]
  cps17<-colSums(am3*ees3[,,17][,ide3,1])[ide3_s[[17]]]
  cps18<-colSums(am3*ees3[,,18][,ide3,1])[ide3_s[[18]]]
  cps19<-colSums(am3*ees3[,,19][,ide3,1])[ide3_s[[19]]]
  cps20<-colSums(am3*ees3[,,20][,ide3,1])[ide3_s[[20]]]

  cps<-abind(cps1,cps2,cps3,cps4,cps5,cps6,cps7,cps8,cps9,
             cps10,cps11,cps12,cps13,cps14,cps15,cps16,
             cps17,cps18,cps19,cps20,along=1)

  sim_net<-cps*tb_3ax

  gross_count_sim<-sim_net + sim_bckg
}

if ("PNCC"%in% used_meas){
  #passive neutron

  mass_pncc<-mass

  comp_pncc<-mass_pncc/greta_array(rep(rowSums(mass_pncc),times=n_nucl),c(nseg,n_nucl))

  w_sum_sf_pair_rate<-rowSums(comp_pncc*qiqr_g)
  m240Pu_eq<-rowSums(mass_pncc)*w_sum_sf_pair_rate

  eff_pncc <- lam*ehpn_g+(1-lam)*eppn_g

  Rexp_m_sim <- Fs_prime_240Pu_g*mean(m240Pu_eq*eff_pncc)*nseg

}
#### Observations and simulated data ####

if("3AX"%in% used_meas){
  gross_count_obs <-c(obs_net + obs_bckg) %>% as_data()
}

if ("PNCC"%in% used_meas){
  Rexp_m_obs<-obs_Rexp_m %>% as_data()
}
#### Likelihood ####
if ("3AX"%in% used_meas){
    distribution(gross_count_obs) <- poisson(gross_count_sim)
}

if ("PNCC"%in% used_meas){
  distribution(Rexp_m_obs) <- normal(Rexp_m_sim,sigma_Rexp_m)
}
#### Compile Greta model ####

if (Use1vec==TRUE){ # a single vector for every segment - not used herein

  if ("3AX"%in% used_meas){ #3AX or 3AX-PNCC
    m <- model(puv,pm,cm244,u235,np237,lam,b_st,am241r,cmr_st,u235r,np237r)
  }else{ # PNCC only
    m <- model(puv,pm,cm244,u235,np237,lam,am241r,cmr_st,u235r,np237r)
  }

}else{ # a different vector is inferred for every segment

  if(prior_for_pu=="gp"){
    if(prior_for_vec!="multilevel"){ # non-multilevel prior for the vector: flat or informative Dirichlet
      if ("3AX"%in% used_meas){ #3AX or 3AX-PNCC
        m <- model(puv1,puv2,puv3,puv4,puv5,puv6,puv7,puv8,puv9,puv10,
                   puv11,puv12,puv13,puv14,puv15,puv16,puv17,puv18,puv19,
                   puv20,pm,cm244,u235,np237,lam,b_st,
                   gp_var,lengthscale,gp_mean,
                   am241r,cmr_st,u235r,np237r)
      }else{ # PNCC only
        m <- model(puv1,puv2,puv3,puv4,puv5,puv6,puv7,puv8,puv9,puv10,
                   puv11,puv12,puv13,puv14,puv15,puv16,puv17,puv18,puv19,
                   puv20,pm,cm244,u235,np237,lam,
                   gp_var,lengthscale,gp_mean,
                   am241r,cmr_st,u235r,np237r)
      }
    }else{ # multilevel prior for the vector: informative Dirichlet with joint inference of the Dirichlet prior distribution parameters
     if ( "3AX"%in% used_meas){
       m <- model(puv1,puv2,puv3,puv4,puv5,puv6,puv7,puv8,puv9,puv10,
                  puv11,puv12,puv13,puv14,puv15,puv16,puv17,puv18,puv19,
                  puv20,pm,cm244,u235,np237,lam,b_st,
                  gp_var,lengthscale,gp_mean,
                  alph,alph_st,
                  am241r,cmr_st,u235r,np237r)
      }else{
        m <- model(puv1,puv2,puv3,puv4,puv5,puv6,puv7,puv8,puv9,puv10,
                   puv11,puv12,puv13,puv14,puv15,puv16,puv17,puv18,puv19,
                   puv20,pm,cm244,u235,np237,lam,
                   gp_var,lengthscale,gp_mean,
                   alph,alph_st,
                   am241r,cmr_st,u235r,np237r)
      }
    }
  }else{# no GP prior for the total Pu mass: (log-)uniform or (log-)Gaussian prior
    if(prior_for_vec!="multilevel"){ # non-multilevel prior for the vector: flat or informative Dirichlet
      if ("3AX"%in% used_meas){ #3AX or 3AX-PNCC
        m <- model(puv1,puv2,puv3,puv4,puv5,puv6,puv7,puv8,puv9,puv10,
                   puv11,puv12,puv13,puv14,puv15,puv16,puv17,puv18,puv19,
                   puv20,pm,cm244,u235,np237,lam,b_st,
                   am241r,cmr_st,u235r,np237r)
      }else{ # PNCC only
        m <- model(puv1,puv2,puv3,puv4,puv5,puv6,puv7,puv8,puv9,puv10,
                   puv11,puv12,puv13,puv14,puv15,puv16,puv17,puv18,puv19,
                   puv20,pm,cm244,u235,np237,lam,
                   am241r,cmr_st,u235r,np237r)
      }
    }else{ # multilevel prior for the vector: informative Dirichlet with joint inference of the Dirichlet prior distribution parameters
      if ( "3AX"%in% used_meas){
        m <- model(puv1,puv2,puv3,puv4,puv5,puv6,puv7,puv8,puv9,puv10,
                   puv11,puv12,puv13,puv14,puv15,puv16,puv17,puv18,puv19,
                   puv20,pm,cm244,u235,np237,lam,b_st,
                   alph,alph_st,
                   am241r,cmr_st,u235r,np237r)
      }else{
        m <- model(puv1,puv2,puv3,puv4,puv5,puv6,puv7,puv8,puv9,puv10,
                   puv11,puv12,puv13,puv14,puv15,puv16,puv17,puv18,puv19,
                   puv20,pm,cm244,u235,np237,lam,
                   alph,alph_st,
                   am241r,cmr_st,u235r,np237r)
      }
    }
  }
}

#### plot model ####
plot(m)
if (SaveModelPlot<-TRUE){
  rsvg::rsvg_png(charToRaw(DiagrammeRsvg::export_svg(plot(m))),
                 file=here(paste0("figures/plot_m_",sampling_case,".png")))
}
#### Sample from the model by HMC ####
if (DoSampling==TRUE){
  if(prior_for_pu=="gp"){

    tic()
    draws <- greta::mcmc(m, warmup = warmup, n_samples = n_samples, pb_update = 100,
                          chains=n_chains,one_by_one = TRUE,sampler=hmc(Lmin = 20, Lmax = 25))
    toc()

  }else{

    tic()
    draws <- greta::mcmc(m, warmup = warmup, n_samples = n_samples, pb_update = 100,
                          chains=n_chains,one_by_one = TRUE,sampler=hmc(Lmin = 20, Lmax = 25))
    toc()

  }

  # Get the corresponding simulated data
  sim<-list()
  tic()
  if ("3AX"%in% used_meas){
    sim[[1]]<-list()
    sim[[1]][[1]]<-calculate(gross_count_sim, values = draws) # deterministic prediction(s) for each posterior parameter set
    sim[[1]][[2]]<-calculate(gross_count_obs, values = draws,nsim=2E3) # posterior predictive model checks
  }else{sim[[1]]<-NA_real_}
  if ("PNCC"%in% used_meas){
    sim[[2]]<-list()
    sim[[2]][[1]]<-calculate(Rexp_m_sim, values = draws) # deterministic posterior predictions
    sim[[2]][[2]]<-calculate(Rexp_m_obs, values = draws,nsim=2E3) # posterior predictive distribution
  }else{sim[[2]]<-NA_real_}
  toc()

  res<-list()
  res[[1]]<-draws
  res[[2]]<-sim
  res[[3]]<-used_meas

  if(DoSynthetic){
    res[[4]]<-ref_par
    res[[5]]<-ref_net
    res[[6]]<-ref_Rexp_m
  }else{
    res[[4]]<-NA_real_
    res[[5]]<-NA_real_
    res[[6]]<-NA_real_
  }
  res[[7]]<-obs_net
  res[[8]]<-obs_bckg
  res[[9]]<-obs_Rexp_m
  res[[10]]<-sigma_Rexp_m
  res[[11]]<-warmup


  if(DoSaveOutput){
  saveRDS(res, file = here(paste0("output/",sampling_case,
                                  "_",vectype,"_",prior_for_vec,
                                  "_prpu_",prior_for_pu,
                                  "_seed",seed,".rds")))


  }
}else{

  res<-readRDS(here(paste0("output/",sampling_case,
                           "_",vectype,"_",prior_for_vec,
                           "_prpu_",prior_for_pu,
                           "_seed",seed,".rds")))

  draws<-res[[1]]
  sim<-res[[2]]
}

#### Quick and dirty trace plots ####
parset<-list.rbind(draws) %>%
  as_tibble()
colnames(draws[[1]])
np<-dim(draws[[1]])[2]
ns<-dim(draws[[1]])[1]
cat(np,ns)
idx<-seq(101,109)
q_start<-1
q_end<-ns
ii<-seq(q_start,q_end,10)
draws_<-draws
for (name in names(draws_)) {
  print(name)
  draws_[[name]]<-draws_[[name]][ii,idx]
}
mcmc_trace(draws_, facet_args = list(ncol = 4))

#### Remove a stucked chain (can happen with the GP prior) ####
RemStuckChain<-TRUE
if(RemStuckChain){ # check on Pu mass for segment 1, col number = 101
  id_bad<-vector()
  teller<-0
  idx<-seq(length(draws))
  for (name in names(draws)) {
    teller<-teller+1
    print(name)
    q<-draws[[name]][,101] #  the check is made based on sampled variable 101 but another one can be chosen
    if(base::sum(diff(q))==0){# chain is stuck
      id_bad<-c(id_bad,teller)
    }
  }
  id_bad
  idx2<- idx[!idx %in% id_bad]
  draws<-draws[idx2]
  n_chains<-length(draws)
  cat(id_bad)
}

#### Convergence of the MCMC ####
rstat<-coda::gelman.diag(draws,autoburnin=F,multivariate=F)
iir<-which(rstat[[1]][,1]>1.1)
cat(length(iir))


