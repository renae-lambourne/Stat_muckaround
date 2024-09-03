x=seq(1, 10, len = 1)
y=40 * 2 + rnorm(10, 0, 5)
plot(x, y)
summary(x)
median(x)
mean(x)

# When collaborating::: commit then pull then push

##Install packages
install.packages("tidyverse")
install.packages("remotes")
install.packages("DHARMa")
install.packages("emmeans")
install.packages("tidybayes")
install.packages("HDInterval")
install.packages("ggeffects")
install.packages("broom.mixed")
install.packages("patchwork")
install.packages("bayestestR")
install.packages("see")
install.packages("easystats")
install.packages("modelsummary")

# install package standist from github as it isn't in cran
remotes::install_github("jmgirard/standist")

#Add libraries
library(tidyverse) #for data wrangling
library(remotes) #to install off github instead of cran
library(cmdstanr) # for cmdstan
library(brms) # for fitting models in STAN
library(standist) #for exploring distributions
library(coda) #for diagnostics
library(bayesplot) #for diagnostics
library(DHARMa) #for residual diagnostics
library(rstan) #for interfacing with stan
library(emmeans) #for marginal means etc
library(broom) #for tidying outputs
library(tidybayes) #for more tidying outputs
library(HDInterval) #for HPD intervals
library(ggeffects) #for partial plots
library(broom.mixed) #for summarising models
library(posterior) #for posterior draws
library(patchwork) #for multi-panel figures
library(bayestestR) #for ROPE
library(see) #for some plots
library(easystats) #framework for stats, modelling and visualisation
library(modelsummary) #for data and model summaries

#set working directory > session> set working directory > this source
# Look at data and check assumptions before fitting model >> Normal, homogeneity of variance, linearity
# ^ to do this need to perform exploratory data analysis

## Exploratory data analysis

boxplot(fert$FERTILIZER)
boxplot(fert$YIELD)

## assessing normality
ggplot(data = fert,aes(x=YIELD))+
  geom_boxplot()
# data ^^ seems normal

## assessing linearity & assessing equal variance
ggplot(data=fert,aes(x=FERTILIZER,y=YIELD))+
  geom_point()+
  geom_smooth(method="lm") #for variance- that the distance between the line and shaded area doesn't increase as the line progresses

#prepare continuous data by centering the data.
fert |> summarise(Median=median(YIELD), MAD=mad(YIELD),MAD_SLOPE=mad(YIELD)/mad(FERTILIZER))
#^ Median and MAD are the names of the columns in the output

#make priors
priors<-prior(normal(162,90),class="Intercept")+
  prior(normal(0,1),class='b')+
  prior(student_t(3,0,90),class='sigma')

#now define the formula of our model
form<-bf(YIELD~FERTILIZER, family=gaussian())

fert.brm<-brm(form,
              data=fert,
              prior=priors,
              sample_prior='only',
              iter=5000,
              warmup=1000,
              chains=3,cores=3,
              thin=5,
              backend='cmdstanr',
              refresh=0)

fert.brm |> 
  conditional_effects() |> 
  plot(points=TRUE)
#^ from above the priors are probably ok so we rerun the model

fert.brm2<-update(fert.brm,sample_prior='yes')

fert.brm2 |> 
  conditional_effects() |> 
  plot(points=TRUE)

fert.brm2 |> SUYR_prior_and_posterior()
#When the coloured point are different to the black points it suggests that the priors have no influence on the data

# diagnostics- check that the chains have been run and it has converged
fert.brm2 |> mcmc_plot(type='trace')
#these plots look good- lots of noise. You don't want to see the plots drifting up or down

# assess autocorrelation- how related are our draws to each other--- we don't want them to be correlated
fert.brm2 |> mcmc_plot(type='acf_bar')
#^ no correlation let--so our thinning was successful!

# measuring convergence-- Rhat (the metric for convergence-- want all Rhat values to be <1.01)
fert.brm2 |> mcmc_plot(type='rhat_hist')
##^ for this all the Rhat values are less than 1.01 so this means that all the chains converged

##sampling efficiency -- ideally you want all the values to be >0.5, the closer to 1 the better
fert.brm2 |> mcmc_plot(type='neff_hist')

#alternate way to look at the above diagnostics
stan_trace(fert.brm2$fit)

#Model Validation
fert.brm2 |> pp_check(type = 'dens_overlay',ndraws=100)

fert.resids<-make_brms_dharma_res(fert.brm2,integerResponse=FALSE)
# there should be no patterns in residuals -- so we need to explore the residuals to see if there are any patterns
testUniformity(fert.resids)
plotResiduals(fert.resids)

fert.brm2 |> summary()
# from ^^ Tail and Bulk ESS you want the values to be greater than 1000



fert.brm2 |> as_draws_df() #returns the full posteriors of every parameter
fert.brm2 |> as_draws_df() |> summarise_draws(median)
fert.brm2 |> as_draws_df() |> summarise_draws(median, HDInterval::hdi)
fert.brm2 |> as_draws_df() |> ggplot(aes(x=b_FERTILIZER))+geom_histogram()

fert.brm2 |> as_draws_df() |> summarise_draws(median,HDInterval::hdi,rhat,ess_bulk,ess_tail,
                                              Pg0=~mean(.x>1)) #hypothesis testing in a bayesian sense

#calculate R2 value for a bayesian model
fert.brm2 |> bayes_R2(summary=FALSE) |> median_hdci()

#Model prediction
fert.brm2 |> predict(newdata=data.frame(FERTILIZER=100))


#Create new dataframe
newdata<-data.frame(FERTILIZER=100)
newdata
fert.brm2 |> emmeans(~FERTILIZER,at=newdata) #prediction
# ^^ so a prediction for fetiliser of 100 would be a yield of 133 

#if want to predict for two values:
newdata<-data.frame(FERTILIZER=c(100,200)) #adds an extra value into the vector
fert.brm2 |> emmeans(~FERTILIZER,at=newdata)

newdata<-with(fert,data.frame(FERTILIZER=seq(min(FERTILIZER),max(FERTILIZER),len=100)))
newdata
fert.brm2 |> emmeans(~FERTILIZER,at=newdata)
fert.pred<-fert.brm2 |> emmeans(~FERTILIZER,at=newdata) |> as.data.frame()
head(fert.pred)

##need the data to be in a dataframe

fig1<-fert.pred |> ggplot(aes(y=emmean,x=FERTILIZER))+
  geom_line()+
  geom_ribbon(aes(ymin=lower.HPD,ymax=upper.HPD),fill="orange",alpha=0.3)+
  theme_classic()+
  scale_y_continuous(expression(Grass~yield~(g.m^-2)),
                     breaks=seq(50,300,by=50))+
  scale_x_continuous(expression(Fertiliser~concentration~(g.m^-2)))+
  geom_point(data=fert,aes(y=YIELD))+
  ggtitle("Relationship between grass yield and fertiliser concentration")
  #use a ~ anywhere you want a space in an expression
fig1
ggsave(filename="figure1.png",fig1,height=5,width=6,dpi=300)


#Methods for ^^
# the relationship between grass yield and fertiliser concentration was explored using linear regression in a Bayesian framework
#specifically, the yield of grass was modelled against fertilizer concentration with a Gaussian family and weakly informative priors
# (equation and priors)
# the model was fit using the brms package in R
# the bayesian model included 3 chains, each of which was run for 5000 iterations, thin to a rate of 5 and a warmup of 1000
# the model was found to be well mixed and converged (all Rhat<1.01) on a stable posterior and was validated via simulated residuals ()
# all statistical models were preformed in the R (4.4.1) Statical and Graphical Environment using the brms package ()







