
# Required packages -----------------
  
library(dplyr)
library(ggplot2)
library(ggarrange)

# JFS_N and JFScut functions are required

# Figure 1--------------------------
# *FIG 1 A -------------------------
retention <- seq(0.5, 1.0, by = 0.001)
samp_size <- ceiling(168/retention)
rate <- (samp_size/36)/2
fig1 <- as.data.frame(cbind(samp_size,rate, retention))
fig1 <- as.data.frame(cbind(samp_size,rate, retention))

add <- c(NA, max(rate), max(retention))
fil_1 <- rbind(fig1, add)

Fig1A <- ggplot(data = fig1, aes(x = retention, y = rate))+
  geom_line() + 
  xlab("Retention") +
  ylab("Recruitment Rate") + theme_bw() +
  geom_polygon(data=fil_1, aes(retention, rate), fill="blue",
               alpha=0.2)

# *FIG 1 B -------------------------

fil_1b <- fil_1[fil_1[,"retention"]>=0.7,]
add2 <- c(NA, max(rate), 0.7)
fil_1b <- rbind(fil_1b, add2)

Fig1B <- ggplot(data = fig1, aes(x = retention, y = rate))+
  geom_line() + 
  xlab("Retention") +
  ylab("") + theme_bw() +
  geom_vline(xintercept = .7, linetype = "dotted")+
  geom_polygon(data=fil_1b, aes(retention, rate), fill="blue", 
               alpha=0.2)

# Combine Fig1A and Fig1B
ggarrange(Fig1A, Fig1B, labels = c("A", "B"),  
          ncol=2, nrow=1, common.legend = TRUE, legend="bottom")


# Fake data creation ---------------

n <- 40
set.seed(123)
lam <- 2.91/4.33
y <- rpois(30,lam) # Simulate weekly recruitment
rate2 <- y[1:22] # Keep weeks needed (22 wks) to get 20 dyads 

set.seed(101)
retention2 <- rbinom(n,1,0.84) # Simulate retention data

set.seed(109)
# Bootstrap
boot1 <- matrix(NA,10000,2)
for (i in 1:10000){
  boot1[i,1] <- mean(sample(rate2,22,replace = TRUE))*4.33 # recruitement rate
  boot1[i,2] <- mean(sample(retention2,n,replace = TRUE))
}

boot1 <- as.data.frame(boot1)


# Figure 2  -------------------------

addcords <- c(NA, 7, 1)
addcords2 <- c(NA, 7, 0.7)

fil_2 <- rbind(fig1, addcords,addcords2)
fil_2 <- fil_2[fil_2[,"retention"]>=0.7,]

ggplot(data = fig1, aes(x = retention, y = rate))+
  geom_line() + 
  xlab("Retention") +
  ylab("Recruitment Rate") + theme_bw() +
  geom_vline(xintercept = .7, linetype = "dotted")+
  geom_polygon(data=fil_2, aes(retention, rate), fill="blue",
               alpha=0.2) +
  geom_jitter(data = boot1, 
              mapping = aes(x = V2, y = V1), size = 0.1) 

# Figure 3--------------------------
# Optimization of cut point via simulation

n <- 40
dyadn <- n/2
# feas_rate <- c(2.91, 2.5)
# feas_ret <- c(0.8, 0.725)

 feas_rate <- c(2.91, 2.19)
 feas_ret <- c(0.8, 0.75)

set.seed(711)
proceedopt <- matrix(NA,20000,2)
for (cond in 1:2){
  feas <- feas_rate[cond]
  ret <- feas_ret[cond]
for (j in 1:10000){

  rate2 <- rpois(200, feas/4.333)
  #Divided by 4.333 to simulate weekly values
  
  # Generates a sample of varying weeks that enrolls a total
  # of 20 participants, this code identified how many weeks were needed
  csum <- cumsum(rate2)
  nl <- which(csum>=dyadn)[1]
  rate2 <- rate2[1:nl]
  diff <- sum(rate2)-dyadn
  rate2[nl] <- rate2[nl]-diff

  # Generate binary retention data
  retention2 <- rbinom(n,1,ret)

  # Bootstrap
  # Retention and rate are defined for Figure 1, which is the
  # same example, this checks for exceeding the JFS boundary
  boot2 <- matrix(NA,1000,3)
  for (i in 1:1000){
    boot2[i,1] <- mean(sample(rate2,nl,replace = TRUE))*4.33 # recruitment rate
    boot2[i,2] <- mean(sample(retention2,n,replace = TRUE))
    boot2[i,3] <- ifelse(boot2[i,2]>=0.7 && boot2[i,1] >= 
             rate[which(round(retention,3) == round(boot2[i,2],3))] ,1,0)
  }
  
  empty_rows <- which(rowSums(is.na(proceedopt)) == ncol(proceedopt))
  empt <- empty_rows[1]
  
  proceedopt[empt,1] <- cond
  proceedopt[empt,2] <- mean(boot2[,3])
} # Simulations loop
} # Condition loop

proceedopt2 <- data.frame(proceedopt)
colnames(proceedopt2) <- c("Combination", "probability")
proceedopt2$Combination <- ifelse(proceedopt2$Combination == 1, "S", "SN")

ggsave("G:\\Papers\\JFS\\Fig32.tiff", units="in", width=6,height=4, dpi=300)

ggplot(proceedopt2,aes(x=probability, fill=Combination)) + geom_density(alpha=0.3) + theme_bw() + ylab("Density") +
  xlab("Probability of being in JFS") +scale_fill_manual(values=c( "#56B4E9","#ffa600"))
dev.off()




# Get cut point:
proceed_S <- proceedopt2[proceedopt2$Combination == "S",]
proceed_SN <- proceedopt2[proceedopt2$Combination == "SN",]

picker <- matrix(NA,1001,4)
picker[,1] <- seq(0,1,0.001)
for (l in 1:1001){
  picker[l,2] <- length(proceed_S[proceed_S[,2]>=picker[l,1],2])/length(proceed_S[,2])
  picker[l,3] <- length(proceed_SN[proceed_SN[,2]>=picker[l,1],2])/length(proceed_SN[,2])
}



# Gives cut point to hit P_G(S) = 0.80
cutpoint <- picker[tail(which(picker[,2]>=0.80), n = 1),1]

length(proceed_S[proceed_S[,2]>= cutpoint,2])/length(proceed_S[,2])
length(proceed_SN[proceed_SN[,2]>= cutpoint,2])/length(proceed_SN[,2])

# Cut points for largest difference P_G(S) and P_G(SN)
picker[,4] <- picker[,2]-picker[,3]
picker[which(picker[,4] == max(picker[,4])),c(1,2,3)]


mean(proceed_S[,2])
sd(proceed_S[,2])

mean(proceed_SN[,2])
sd(proceed_SN[,2])

quantile(proceed_SN[,2])
quantile(proceed_S[,2])


# Simulation study 1 endpoint ---------------
# Simulation to compare, Threshold, Tolerance, & Hyp. Test.

n4sim1 <- 40
threshold <- 0.8
effect <- 0.8
tolerance <- 0.05

# **Threshold theoretical values----------------
sum(dbinom(32:40,40,.8)) #P_FA(S) = 0.59
sum(dbinom(32:40,40,.7)) # P_FA(SN) = 0.11
# For N = 119
sum(dbinom(96:119,119,.8)) #P_FA(S) = 0.48
sum(dbinom(96:119,119,.7)) # P_FA(SN) = 0.006

# Optimization of methods before simulation:
# **Optimize Tolerance ----------------
# For N = 40
sum(dbinom(30:40,40,.8))
# Need to observe 30/40 success to have P_FA(S) >= 0.80 (0.839)
# 30/40 = 0.75, so tolerance = 0.05
# We would expect P_FA(SN) = 0.31
sum(dbinom(30:40,40,.7))
 
 # For N = 119
 sum(dbinom(92:119,119,.8))
 # Need to observe 92/119 success to have P_FA(S) >= 0.80 (0.803)
 # 92/119 = 0.7731092, so tolerance = 0.0268908
 # We would expect P_FA(SN) = 0.047
 sum(dbinom(92:119,119,.7))


 # **Optimize Hyp. Test ----------------
 # Hyp. test procedure, 2 options:
 # 1) Simple: E >= GLL -> Go, otherwise Don't (not considering amend proceed)
 # this will have same operating char. as threshold
 # 2) Expanded: (E >= GLL) -> Go or (E >= RUL & pval < 0.05) -> Go
 
 
 threshold <- 0.8
 effect <- 0.8
 set.seed(229)
 checkHT1 <- matrix(NA,5000,4)
 for (f in 1:5000){
   samp80 <- rbinom(40,1,0.8)
   samp70 <- rbinom(40,1,0.7)
   
   # Simple
   checkHT1[f,1] <- ifelse(mean(samp80)>=threshold,1,0)
   checkHT1[f,2] <- ifelse(mean(samp70)>=threshold,1,0)
   
   # Expanded
   # This calculates if either the estimate is in the green zone,
   # or the p-value is < 0.05, while the estimate is not in RED zone
   NotinRED <- ifelse(mean(samp80) > 0.70,1,0)
   p80 <- ifelse(binom.test(sum(samp80), 40, p = 0.70, alternative = "greater")$p.value < 0.05,1,0)
    Amber_minor <- ifelse(NotinRED + p80 >1,1,0)
    GO80 <- ifelse((checkHT1[f,1]+Amber_minor)>0 ,1,0)
   checkHT1[f,3] <- GO80
   
   NotinRED70 <- ifelse(mean(samp70) > 0.70,1,0)
   p70 <- ifelse(binom.test(sum(samp70), 40, p = 0.70, 
                            alternative = "greater")$p.value < 0.05,1,0)
   Amber_minor70 <- ifelse(NotinRED70 + p70 >1,1,0)
   GO70 <- ifelse((checkHT1[f,2]+Amber_minor70)>0 ,1,0)
   checkHT1[f,4] <- GO70
  
 }
 colMeans(checkHT1)
 # We get the exact same operating characteristics as the threshold approach,
 # because the p-value for the binomial exact test cannot be < 0.05, unless 
 # we are already in the green zone, so the amended approach reverts to the 
 # simple approach
 
 # We'd like to get fair comparison between FA approach and Hyp. Test so we set the 
 # Type I error rate for the Hyp. Test to 0.31, this gives as close to P_FA(S) = 0.80 
 # as possible,
 # it also happens to correspond exactly with the Tolerance approach, since we get a
 # p-value < 0.31
 # iff we observe at least 30 successes (30/40 = 0.75), simulation below:
 
 #For N = 40
 
 set.seed(239)
 threshold <- 0.8
 effect <- 0.8
 checkHT1o <- matrix(NA,5000,4)
 for (f in 1:5000){
   samp80 <- rbinom(40,1,0.8)
   samp70 <- rbinom(40,1,0.7)
   
   # Simple
   checkHT1o[f,1] <- ifelse(mean(samp80)>=threshold,1,0)
   checkHT1o[f,2] <- ifelse(mean(samp70)>=threshold,1,0)
   
   # Expanded
   # This calcualtes if either the estimate is in the green zone,
   # or the p-value is < 0.05, while the estimate is not in RED zone
   NotinRED <- ifelse(mean(samp80) > 0.70,1,0)
   p80 <- ifelse(binom.test(sum(samp80), 40, p = 0.70, 
                            alternative = "greater")$p.value < 0.31,1,0)
   Amber_minor <- ifelse(NotinRED + p80 >1,1,0)
   GO80 <- ifelse((checkHT1o[f,1]+Amber_minor)>0 ,1,0)
   checkHT1o[f,3] <- GO80
   
   NotinRED70 <- ifelse(mean(samp70) > 0.70,1,0)
   p70 <- ifelse(binom.test(sum(samp70), 40, p = 0.70, 
                            alternative = "greater")$p.value < 0.31,1,0)
   Amber_minor70 <- ifelse(NotinRED70 + p70 >1,1,0)
   GO70 <- ifelse((checkHT1o[f,2]+Amber_minor70)>0 ,1,0)
   checkHT1o[f,4] <- GO70
  
 }
 colMeans(checkHT1o)
 
 #For N = 119
 
 set.seed(239)
 threshold <- 0.8
 effect <- 0.8
 checkHT1o2 <- matrix(NA,5000,4)
 for (f in 1:5000){
   samp80 <- rbinom(119,1,0.8)
   samp70 <- rbinom(119,1,0.7)
   
   # Simple
   checkHT1o2[f,1] <- ifelse(mean(samp80)>=threshold,1,0)
   checkHT1o2[f,2] <- ifelse(mean(samp70)>=threshold,1,0)
   
   # Expanded
   # This calcualtes if either the estimate is in the green zone,
   # or the p-value is < 0.05, while the estimate is not in RED zone
   NotinRED <- ifelse(mean(samp80) > 0.70,1,0)
   p80 <- ifelse(binom.test(sum(samp80), 119, p = 0.70, 
                            alternative = "greater")$p.value < 0.05,1,0)
   Amber_minor <- ifelse(NotinRED + p80 >1,1,0)
   GO80 <- ifelse((checkHT1o2[f,1]+Amber_minor)>0 ,1,0)
   checkHT1o2[f,3] <- GO80
   
   NotinRED70 <- ifelse(mean(samp70) > 0.70,1,0)
   p70 <- ifelse(binom.test(sum(samp70), 119, p = 0.70, 
                            alternative = "greater")$p.value < 0.05,1,0)
   Amber_minor70 <- ifelse(NotinRED70 + p70 >1,1,0)
   GO70 <- ifelse((checkHT1o2[f,2]+Amber_minor70)>0 ,1,0)
   checkHT1o2[f,4] <- GO70
   
 }
 colMeans(checkHT1o2)
 # For simple, we get the threshold approach as expected, for Expanded we get 
 # the tolerance
 # approach values (within simulation error) of 0.8346 and 0.3168 
 
 # **Optimize JFS  ---------------------------
 # Need to optimize the cut point probability
 # (Note: This code takes ~ 6 minutes to run)
 n4sim1 <- 40 # Change to 119 to run for N = 119
 #a <- Sys.time()
 feas_ret <- c(0.70, 0.80)
 set.seed(249)
 FAsimopt <- matrix(NA,10000,2)
 for (cond in 1:2){
   ret <- feas_ret[cond]
   for (j in 1:5000){
     # Generate binary retention data
     retentionE <- rbinom(n4sim1,1,ret)
     
     # Bootstrap
     boots1o <- matrix(NA,1000,1)
     for (i in 1:1000){
       boots1o[i,1] <- ifelse(mean(sample(retentionE,n4sim1,replace = TRUE))>=0.8,1,0) # recruitment rate
     }
     
     empty_rows <- which(rowSums(is.na(FAsimopt)) == ncol(FAsimopt))
     empt <- empty_rows[1]
     
     FAsimopt[empt,1] <- cond
     FAsimopt[empt,2] <- mean(boots1o)
     
    # print(empt)
   } # Simulations loop
 } # Condition loop
 #Sys.time()-a
 
 # Cut point to discriminate
 
 breaklist <- c(seq(0,1, 0.001))
 FAsimopt_S <- FAsimopt[FAsimopt[,1] == 2,]
 FAsimopt_SN <- FAsimopt[FAsimopt[,1] == 1,]
 
 # Check:
 picker <- matrix(NA,1001,3)
 picker[,1] <- seq(0,1,0.001)
 for (l in 1:1001){
   picker[l,2] <- length(FAsimopt_S[FAsimopt_S[,2]>=picker[l,1],2])/length(FAsimopt_S[,2])
   picker[l,3] <- length(FAsimopt_SN[FAsimopt_SN[,2]>=picker[l,1],2])/length(FAsimopt_SN[,2])
 }

 cutpoint <- picker[tail(which(picker[,2]>=0.80), n = 1),1]
 # Gives cut point of 0.294
 # P_FA(S) = 0.8022 , P_FA(SN) = 0.2704
 # We'll round down to 0.29, due to possible simulation error and the desire to 
 # exceed 80%
 # For comparison to tolerance and HT, a cut point of 0.22 
 # gives P_FA(SN) = 0.3084
 
 # Re-run with N = 119 (Table 2 in HT paper) with same seed,
 # we get a cut point of 0.19
 # P_FA(S) = 0.8000 , P_FA(SN) = 0.0524
 # For comparison to HT, 0.207

# **Simulation 1 code ---------------
#a <- Sys.time()
 nsim1vect <- c(40, 119)
 effectv <- c(0.7,0.8)

 out_sim1 <- matrix(NA,20000,8)
 for (nv in 1:2){
   nsim1 <- nsim1vect[nv]
   
   for (ev in 1:2){
     ret_effect <- effectv[ev]
     
     for (k in 1:5000){
       
       empty_rows <- which(rowSums(is.na(out_sim1)) == ncol(out_sim1))
       empt <- empty_rows[1]
       out_sim1[empt,1] <- nsim1
       out_sim1[empt,2] <- ret_effect
       
       sampls <- rbinom(nsim1,1,ret_effect)
       
       # Threshold
       out_sim1[empt,3] <- ifelse(mean(sampls)>=0.8,1,0)
       
       # Tolerance
       thresh <- ifelse(nsim1 == 40, 0.05,0.0268908 )
       out_sim1[empt,4] <- ifelse(mean(sampls)>=0.80-thresh,1,0)
       
       # Hypothesis test simple (same as threshold)
       out_sim1[empt,5] <- ifelse(mean(sampls)>=0.8,1,0)

       pcut <- ifelse(nsim1 == 40, 0.31, 0.05)
       NotinRED <- ifelse(mean(sampls) > 0.70,1,0)
       pval <- ifelse(binom.test(sum(sampls), nsim1, p = 0.70, 
                                alternative = "greater")$p.value < pcut,1,0)
       Amber_minor <- ifelse(NotinRED + pval >1,1,0)
       
       # HT Expanded
       out_sim1[empt,6] <- ifelse((out_sim1[empt,3]+Amber_minor) > 0 ,1,0)
       
       # JFS
       bootsim1 <- matrix(NA,1000,1)
       for (i in 1:1000){
         bootsim1[i,1] <- ifelse(mean(sample(sampls,nsim1,replace = TRUE))>=0.8,1,0) 
       }
       
      # Cut point to achieve P_{G}(S) >= 0.80
       disc <- ifelse(nsim1 == 40,0.29, 0.19)
       
       out_sim1[empt,7] <- ifelse(mean(bootsim1) >= disc,1,0)
       
       # Specified cut point to match with HT
       directcompare <- ifelse(nsim1 == 40,0.22, 0.207)
       out_sim1[empt,8] <- ifelse(mean(bootsim1) >= directcompare,1,0)
       
    print(empt)
     }
   }
 }
 totaltime <- Sys.time()-a
 out_sim1work <- data.frame(out_sim1)
 colnames(out_sim1work) <- c("N", "effect", "Threshold", "Tolerance", "HTS", 
                             "HTE", "FA", "FAfair")
 
 out_sim1work %>% dplyr::group_by(N, effect) %>% 
   dplyr::summarise(across(everything(), list(mean))) -> tablesim1
 

 # Simulation study 2 endpoints -----------------------
 
 # All optimizations below assume independence
 
 # **Threshold theoretical values------------------
 #n = 40
 sum(dbinom(32:40,40,.8))^2
 sum(dbinom(32:40,40,.7))^2
 # Both hitting threshold
 # P_{G}(S) = 0.35
 # P_{G}(SN) = 0.012
 #n = 119
 sum(dbinom(96:119,119,.8))^2
 sum(dbinom(96:119,119,.7))^2
 # Both hitting threshold
 # P_{G}(S) = 0.23
 # P_{G}(SN) = < 0.0001
 
 # **Optimize tolerance ------------------
 # For N = 40
 sum(dbinom(29:40,40,.8))^2
 # Need to observe 29/40 success  in both endpoints to have probability of both
 # being in P_{G}(S) >= 0.80 (0.8326467)
 # 29/40 = 0.725, so tolerance = 0.075
 # We would expect P_{G}(SN) = 0.19
 sum(dbinom(29:40,40,.7))^2
 
 # For N = 119
 sum(dbinom(90:119,119,.8))^2 # 0.8135
 # Need to observe 90/119 success in both endpoints
 # 90/119 = 0.7563025, so tolerance = 0.0436975
 # We would expect P_{G}(SN) = 0.01
 sum(dbinom(90:119,119,.7))^2
 
 
 # **Optimize Hypothesis test
 
 # **Optimize Hyp. Test ----------------
 # Hyp. test procedure, 2 options:
 # 1) Simple: E >= GLL -> Go, otherwise Don't (not considering amend proceed)
 # this will have same operating char. as threshold
 # 2) Expanded: (E >= GLL) -> Go or (E >= RUL & pval < alpha threshold) -> Go
 
 # N = 40
 set.seed(429)
 pcut <- 0.45 # Specified alpha level to achieve P_{G}(S) >= 0.80
 checkHT2 <- matrix(NA,10000,4)
 for (f in 1:10000){
   samp801 <- rbinom(40,1,0.8)
   samp802 <- rbinom(40,1,0.8)
   
   samp701 <- rbinom(40,1,0.7)
   samp702 <- rbinom(40,1,0.7)
   
   # Simple (all within green zone)
   checkHT2[f,1] <- ifelse(mean(samp801)>=0.80 & mean(samp802) >= 0.80,1,0)
   checkHT2[f,2] <- ifelse(mean(samp701)>=0.80 & mean(samp702) >= 0.80,1,0)
   
   # Expanded
   # This calculates if either the estimate is in the green zone,
   # or the p-value is < alpha level, while the estimate is not in RED zone
   # Both endpoints need to be Go for us to proceed
   NotinREDs1 <- ifelse(mean(samp801) > 0.70,1,0)
   p801 <- ifelse(binom.test(sum(samp801), 40, p = 0.70, 
                             alternative = "greater")$p.value < pcut,1,0)
   Amber_minor1 <- ifelse(NotinREDs1 + p801 >1,1,0)
   GO801 <- ifelse((ifelse(mean(samp801)>=0.80,1,0)+Amber_minor1)>0 ,1,0)
   
   NotinREDs2 <- ifelse(mean(samp802) > 0.70,1,0)
   p802 <- ifelse(binom.test(sum(samp802), 40, p = 0.70, 
                             alternative = "greater")$p.value < pcut,1,0)
   Amber_minor2 <- ifelse(NotinREDs2 + p802 >1,1,0)
   GO802 <- ifelse((ifelse(mean(samp802)>=0.80,1,0)+Amber_minor2)>0 ,1,0)
   
   checkHT2[f,3] <- ifelse(GO801 + GO802==2,1,0)   
   
   NotinRED701 <- ifelse(mean(samp701) > 0.70,1,0)
   p701 <- ifelse(binom.test(sum(samp701), 40, p = 0.70, 
                             alternative = "greater")$p.value < pcut,1,0)
   Amber_minor701 <- ifelse(NotinRED701 + p701 >1,1,0)
   GO701 <- ifelse((ifelse(mean(samp701)>=0.80,1,0)+Amber_minor701)>0 ,1,0)
   
   NotinRED702 <- ifelse(mean(samp702) > 0.70,1,0)
   p702 <- ifelse(binom.test(sum(samp702), 40, p = 0.70, 
                             alternative = "greater")$p.value < pcut,1,0)
   Amber_minor702 <- ifelse(NotinRED702 + p702 >1,1,0)
   GO702 <- ifelse((ifelse(mean(samp702)>=0.80,1,0)+Amber_minor702)>0 ,1,0)
   
   checkHT2[f,4] <- ifelse(GO701 + GO702==2,1,0)   
   
 }
 colMeans(checkHT2)
 
 # p-value threshold of 0.45
 # For N = 40, 2 endpoints, need p-value < 0.45 in Amber
 # Get probabilities of 0.8338 and 0.19 (same as tolerance)

 # For N = 119
 set.seed(478)
 pcut <- 0.11 # Alpha level to achieve P_{G}(S) >= 0.80
 checkHT2 <- matrix(NA,10000,4)
 for (f in 1:10000){
   samp801 <- rbinom(119,1,0.8)
   samp802 <- rbinom(119,1,0.8)
   
   samp701 <- rbinom(119,1,0.7)
   samp702 <- rbinom(119,1,0.7)
   
   # Simple (all within green zone)
   checkHT2[f,1] <- ifelse(mean(samp801)>=0.80 & mean(samp802) >= 0.80,1,0)
   checkHT2[f,2] <- ifelse(mean(samp701)>=0.80 & mean(samp702) >= 0.80,1,0)
   
   # Expanded
   # This calculates if either the estimate is in the green zone,
   # or the p-value is < alpha level, while the estimate is not in RED zone
   # Both endpoints need to be Go for us to proceed
   
   NotinREDs1 <- ifelse(mean(samp801) > 0.70,1,0)
   p801 <- ifelse(binom.test(sum(samp801), 119, p = 0.70, 
                             alternative = "greater")$p.value < pcut,1,0)
   Amber_minor1 <- ifelse(NotinREDs1 + p801 >1,1,0)
   GO801 <- ifelse((ifelse(mean(samp801)>=0.80,1,0)+Amber_minor1)>0 ,1,0)
   
   NotinREDs2 <- ifelse(mean(samp802) > 0.70,1,0)
   p802 <- ifelse(binom.test(sum(samp802), 119, p = 0.70, 
                             alternative = "greater")$p.value < pcut,1,0)
   Amber_minor2 <- ifelse(NotinREDs2 + p802 >1,1,0)
   GO802 <- ifelse((ifelse(mean(samp802)>=0.80,1,0)+Amber_minor2)>0 ,1,0)
   
   checkHT2[f,3] <- ifelse(GO801 + GO802==2,1,0)   
   
   NotinRED701 <- ifelse(mean(samp701) > 0.70,1,0)
   p701 <- ifelse(binom.test(sum(samp701), 119, p = 0.70, 
                             alternative = "greater")$p.value < pcut,1,0)
   Amber_minor701 <- ifelse(NotinRED701 + p701 >1,1,0)
   GO701 <- ifelse((ifelse(mean(samp701)>=0.80,1,0)+Amber_minor701)>0 ,1,0)
   
   NotinRED702 <- ifelse(mean(samp702) > 0.70,1,0)
   p702 <- ifelse(binom.test(sum(samp702), 119, p = 0.70, 
                             alternative = "greater")$p.value < pcut,1,0)
   Amber_minor702 <- ifelse(NotinRED702 + p702 >1,1,0)
   GO702 <- ifelse((ifelse(mean(samp702)>=0.80,1,0)+Amber_minor702)>0 ,1,0)
   
   checkHT2[f,4] <- ifelse(GO701 + GO702==2,1,0)   
   
 }
 colMeans(checkHT2)
 
 # p-value threshold of 0.11
 # For simple, we get the threshold approach as expected, for
 # Expanded we get the tolerance
 # approach values (within simulation error) of 0.8178  and 0.01 
 
 # **Optimize Joint Feasibility Space ---------------------------
 # Need to optimize the cut point probability
 # (Note: This code takes ~ 15 minutes to run)
 #a <- Sys.time()
 feas_ret <- c(0.70, 0.80)
 set.seed(632)
 FAsimopt2 <- matrix(NA,10000,2)
 n <- 119
 for (cond in 1:2){
   ret <- feas_ret[cond]
   for (j in 1:5000){
     # Generate binary retention data
     retentionE <- rbinom(n,1,ret)
     retentionE2 <- rbinom(n,1,ret)
     
     # Bootstrap
     # Retention and rate are defined for Figure 1, which is the
     # same example, this checks for exceeding the JFS boundary
     boots2o <- matrix(NA,1000,3)
     for (i in 1:1000){
       resamprows <-sample(1:n,n,replace = TRUE) 
       boots2o[i,1] <- ifelse(mean(retentionE[resamprows])>=0.8,1,0) # recruitment rate
       boots2o[i,2] <- ifelse(mean(retentionE2[resamprows])>=0.8,1,0) # recruitment rate
       boots2o[i,3] <- ifelse(boots2o[i,1] + boots2o[i,2] ==2,1,0) # recruitment rate
     }
     
     empty_rows <- which(rowSums(is.na(FAsimopt2)) == ncol(FAsimopt2))
     empt <- empty_rows[1]
     
     FAsimopt2[empt,1] <- cond
     FAsimopt2[empt,2] <- mean(boots2o[,3])
     
    # print(empt)
   } # Simulations loop
 } # Condition loop
 #Sys.time()-a
 
 # Cut point to discriminate
 breaklist <- c(seq(0,1, 0.001))
 FAsimopt_S <- FAsimopt2[FAsimopt2[,1] == 2,]
 FAsimopt_SN <- FAsimopt2[FAsimopt2[,1] == 1,]
 
 # Check:
 picker <- matrix(NA,1001,3)
 picker[,1] <- seq(0,1,0.001)
 for (l in 1:1001){
   picker[l,2] <- length(FAsimopt_S[FAsimopt_S[,2]>=picker[l,1],2])/length(FAsimopt_S[,2])
   picker[l,3] <- length(FAsimopt_SN[FAsimopt_SN[,2]>=picker[l,1],2])/length(FAsimopt_SN[,2])
 }
 
 cutpoint <- picker[tail(which(picker[,2]>=0.80), n = 1),1]
 picker[tail(which(picker[,2]>=0.80), n = 1),2]
 picker[tail(which(picker[,2]>=0.80), n = 1),3]
 # Gives cut point of 0.096
 # P_FA(S) = 0.8016 , P_FA(SN) = 0.129

 # Re-run with N = 119 (Table 2 in HT paper) with same seed,
 # we get a cut point of 0.049
 # P_FA(S) = 0.8006 , P_FA(SN) = 0.0034
 
 
 
 # **Simulation 2 code ------------------
 # Note before running: This code takes ~22 minutes to run on a single core
 
 set.seed(5473)
 nsim1vect <- c(40, 119)
 effectv <- c(0.7,0.8)
 rhos <- c(0)
 
 out_sim2 <- matrix(NA,(2*2*3*5000),9)

 at <- Sys.time()
 for (nv in 1:2){
   nsim1 <- nsim1vect[nv]
   
   for (ev in 1:2){
     ret_effect <- effectv[ev]
     
     for (rv in 1:1){
       rho <- rhos[rv]
       # Calculate probabilities for correlated binomial data
       a <- function(rho, p , q ) {
         rho * sqrt(p*q*(1-p)*(1-q)) + (1-p)*(1-q)
       }
       
       # same assumed effect in both endpoints
       p <- ret_effect
       q <- ret_effect
       a.0 <- a(rho, ret_effect, ret_effect)
       
       prob <- c(`(0,0)`=a.0, `(1,0)`=1-q-a.0, `(0,1)`=1-p-a.0, `(1,1)`=a.0+p+q-1)
       if (min(prob) < 0) {
         print(prob)
         stop("Error: a probability is negative.")
       }
       

       for (p in 1:5000){
         
         # Generate correlated binomial data
         samps <- sample.int(4, nsim1, replace=TRUE, prob=prob)
         #table(samps)/nsim1
         sams2 <- matrix(NA,nsim1,2)
         
         for (jj in 1:length(samps)){
           if( samps[jj] == 1){
             sams2[jj,1] <- 0
             sams2[jj,2] <- 0
             
           } else if (samps[jj]==2){
             sams2[jj,1] <- 1
             sams2[jj,2] <- 0
           } else if (samps[jj]==3){
             sams2[jj,1] <- 0
             sams2[jj,2] <- 1
           }else if (samps[jj] == 4){
             sams2[jj,1] <- 1
             sams2[jj,2] <- 1
           }
         }
         
         empty_rows <- which(rowSums(is.na(out_sim2)) == ncol(out_sim2))
         empt <- empty_rows[1]
         
         out_sim2[empt,1] <- nsim1
         out_sim2[empt,2] <- ret_effect
         out_sim2[empt,3] <- rho
         
         # Threshold
         e1 <- ifelse(mean(sams2[,1])>=0.8,1,0)
         e2 <- ifelse(mean(sams2[,2])>=0.8,1,0)
         out_sim2[empt,4] <- ifelse(e1+e2 == 2,1,0)
         
         # Tolerance
         tol <- ifelse(nsim1 == 40, 0.075,0.0436975)
         rounder <- ifelse(nsim1 == 40, 3,9) # to avoid floating point issues
         e1_tol <- ifelse(round(mean(sams2[,1]),rounder)>=round((0.8-tol),rounder),1,0)
         e2_tol <- ifelse(round(mean(sams2[,2]),rounder)>=round((0.8-tol),rounder),1,0)
         out_sim2[empt,5] <- ifelse(e1_tol+e2_tol==2,1,0) 
         
         # Hypothesis Test
         # Simple (same as threshold)
         out_sim2[empt,6] <- ifelse(e1+e2 == 2,1,0)

         # Expanded
         # This calculates if either the estimate is in the green zone,
         # or the p-value is < 0.05, while the estimate is not in RED zone
         # Both endpoints need to be Go for us to proceed
         
         pcut <- ifelse(nsim1 == 40, 0.45, 0.11 )
         NotinREDs1 <- ifelse(mean(sams2[,1]) > 0.70,1,0)
         p801 <- ifelse(binom.test(sum(sams2[,1]), nsim1, p = 0.70, 
                                   alternative = "greater")$p.value < pcut,1,0)
         Amber_minor1 <- ifelse(NotinREDs1 + p801 ==2,1,0)
         GO801 <- ifelse((ifelse(mean(sams2[,1])>=0.80,1,0)+Amber_minor1)>0 ,1,0)
         
         NotinREDs2 <- ifelse(mean(sams2[,2]) > 0.70,1,0)
         p802 <- ifelse(binom.test(sum(sams2[,2]), nsim1, p = 0.70, 
                                   alternative = "greater")$p.value < pcut,1,0)
         Amber_minor2 <- ifelse(NotinREDs2 + p802 ==2,1,0)
         GO802 <- ifelse((ifelse(mean(sams2[,2])>=0.80,1,0)+Amber_minor2)>0 ,1,0)
         
         # Expanded
         out_sim2[empt,7] <- ifelse(GO801 + GO802==2,1,0)   
         
         
         # JFS 
         boot2 <- matrix(NA,1000,1)
         for (l in 1:1000){
           rowsamp <- sample(1:nsim1, nsim1, replace = TRUE)
           boot2[l,1] <- ifelse(mean(sams2[rowsamp,1]) >= 0.8 & 
                                mean(sams2[rowsamp,2]) >= 0.8,1,0)
         }
         
         #Straight
         disc1 <- ifelse(nsim1 == 40,0.096, 0.049)
         #Fair
         disc2 <- ifelse(nsim1 == 40, 0.062, 0.028)
         # JFS cut point to get P_{G}(S) >= 0.80
         out_sim2[empt,8] <- ifelse( mean(boot2) >= disc1,1,0)   
         
         # JFS cut point used in paper (this is matched to HT)
         out_sim2[empt,9] <- ifelse( mean(boot2) >= disc2,1,0)   
         
       } # p loop
     } # rho loop
   } # effect loop
 } #n loop
#Sys.time()-at
 

out_sim2work <- data.frame(out_sim2)
colnames(out_sim2work) <- c("N", "effect", "rho", "Threshold", "Tolerance", "HTS", 
                            "HTE", "FA", "FAfair")

out_sim2work %>% dplyr::group_by(N, effect, rho) %>% 
  dplyr::summarise(across(everything(), list(mean))) ->tablesim2


# Simulation study 3 endpoints -------------------
# All optimizations below assume independence

# **Threshold theoretical values------------------
#n = 40
sum(dbinom(32:40,40,.8))^3
sum(dbinom(32:40,40,.7))^3
# All hitting threshold
# P_{G}(S) = 0.21
# P_{G}(SN) = 0.001

#n = 119
sum(dbinom(96:119,119,.8))^3
sum(dbinom(96:119,119,.7))^3
# All hitting threshold
# P_{G}(S) = 0.11
# P_{G}(SN) = < 0.0001


# **Optimize tolerance ------------------
# This tolerance is chosen to achieve P_G(S) = 0.76 to match the HT
# For N = 40
sum(dbinom(29:40,40,.8))^3
# Need to observe 29/40 success in all endpoints to have probability of all
# being in P_{G}(S) ~ 0.76 (0.7597857)
# 29/40 = 0.725, so tolerance = 0.075
# We would expect P_{G}(SN) = 0.085
sum(dbinom(29:40,40,.7))^3

# For N = 119
sum(dbinom(89:119,119,.8))^3 # 0.815
# Need to observe 89/119 success in all endpoints
# 89/119 = 0.7478992, so tolerance = 0.0521008
# We would expect P_{G}(SN) = 0.003
sum(dbinom(89:119,119,.7))^3

# With N = 40 and 3 endpoints the Hypothesis test is not able to achieve
# P_{G}(S) >= 0.80 because the probability of 3 binomials with p = 0.8 all 
# having mean > 0.7 (i.e., not in the red zone is only ~0.75, so the
# p-value for the amber zone is meaningless)
# For the simulation, we'll let the p-value threshold be = 1
tol <- 0.0521008
checktolerance3 <- matrix(NA,5000,2)
for (f in 1:5000){
  samp801 <- rbinom(119,1,0.8)
  samp802 <- rbinom(119,1,0.8)
  samp803 <- rbinom(119,1,0.8)
  
  samp701 <- rbinom(119,1,0.7)
  samp702 <- rbinom(119,1,0.7)
  samp703 <- rbinom(119,1,0.7)

  
  checktolerance2[f,1] <- ifelse(round(mean(samp801),4)>=round((0.8-tol),4) & 
                                   round(mean(samp802),4)>=round((0.8-tol),4) &
                                   round(mean(samp803),4)>=round((0.8-tol),4),1,0)
  checktolerance2[f,2] <- ifelse(round(mean(samp701),4)>=round((0.8-tol),4) &
                                   round(mean(samp702),4)>=round((0.8-tol),4) &
                                   round(mean(samp703),4)>=round((0.8-tol),4) ,1,0)
  
}
colMeans(checktolerance2)


# **Optimize Hypothesis test

# **Optimize Hyp. Test ----------------
# Hyp. test procedure, 2 options:
# 1) Simple: E >= GLL -> Go, otherwise Don't (not considering amend proceed)
# this will have same operating char. as threshold
# 2) Expanded: (E >= GLL) -> Go or (E >= RUL & pval < 0.05) -> Go
# N = 40
set.seed(429)

pcut <- 1.0
checkHT2 <- matrix(NA,2000,4)
for (f in 1:2000){
  samp801 <- rbinom(40,1,0.8)
  samp802 <- rbinom(40,1,0.8)
  samp803 <- rbinom(40,1,0.8)
  
  samp701 <- rbinom(40,1,0.7)
  samp702 <- rbinom(40,1,0.7)
  samp703 <- rbinom(40,1,0.7)
  
  # Simple (all within green zone)
  checkHT2[f,1] <- ifelse(mean(samp801)>=0.80 & mean(samp802) >= 0.80  & mean(samp803) >= 0.80,1,0)
  checkHT2[f,2] <- ifelse(mean(samp701)>=0.80 & mean(samp702) >= 0.80  & mean(samp703) >= 0.80,1,0)
  
  # Expanded
  # This calculates if either the estimate is in the green zone,
  # or the p-value is < alpha level, while the estimate is not in RED zone
  # Both endpoints need to be Go for us to proceed
  NotinREDs1 <- ifelse(mean(samp801) > 0.70,1,0)
  p801 <- ifelse(binom.test(sum(samp801), 40, p = 0.70, 
                            alternative = "greater")$p.value < pcut,1,0)
  Amber_minor1 <- ifelse(NotinREDs1 + p801 >1,1,0)
  GO801 <- ifelse((ifelse(mean(samp801)>=0.80,1,0)+Amber_minor1)>0 ,1,0)
  
  NotinREDs2 <- ifelse(mean(samp802) > 0.70,1,0)
  p802 <- ifelse(binom.test(sum(samp802), 40, p = 0.70, 
                            alternative = "greater")$p.value < pcut,1,0)
  Amber_minor2 <- ifelse(NotinREDs2 + p802 >1,1,0)
  GO802 <- ifelse((ifelse(mean(samp802)>=0.80,1,0)+Amber_minor2)>0 ,1,0)
  
  NotinREDs3 <- ifelse(mean(samp803) > 0.70,1,0)
  p803 <- ifelse(binom.test(sum(samp803), 40, p = 0.70, 
                            alternative = "greater")$p.value < pcut,1,0)
  Amber_minor3 <- ifelse(NotinREDs3 + p803 >1,1,0)
  GO803 <- ifelse((ifelse(mean(samp803)>=0.80,1,0)+Amber_minor3)>0 ,1,0)
  
  checkHT2[f,3] <- ifelse(GO801 + GO802 + GO803==3,1,0)   
  
  NotinRED701 <- ifelse(mean(samp701) > 0.70,1,0)
  p701 <- ifelse(binom.test(sum(samp701), 40, p = 0.70, 
                            alternative = "greater")$p.value < pcut,1,0)
  Amber_minor701 <- ifelse(NotinRED701 + p701 >1,1,0)
  GO701 <- ifelse((ifelse(mean(samp701) >= 0.80,1,0)+Amber_minor701) > 0 ,1,0)
  
  NotinRED702 <- ifelse(mean(samp702) > 0.70,1,0)
  p702 <- ifelse(binom.test(sum(samp702), 40, p = 0.70, 
                            alternative = "greater")$p.value < pcut,1,0)
  Amber_minor702 <- ifelse(NotinRED702 + p702 >1,1,0)
  GO702 <- ifelse((ifelse(mean(samp702)>=0.80,1,0)+Amber_minor702)>0 ,1,0)
  
  NotinRED703 <- ifelse(mean(samp703) > 0.70,1,0)
  p703 <- ifelse(binom.test(sum(samp703), 40, p = 0.70, 
                            alternative = "greater")$p.value < pcut,1,0)
  Amber_minor703 <- ifelse(NotinRED703 + p703 >1,1,0)
  GO703 <- ifelse((ifelse(mean(samp703)>=0.80,1,0)+Amber_minor703)>0 ,1,0)
  
  checkHT2[f,4] <- ifelse(GO701 + GO702 + GO703==3,1,0)   
  
}
colMeans(checkHT2)

# p-value threshold of 1.0
# For N = 40, 3 endpoints the probability of all 3 endpoints
# not being in the red zone is roughly 0.76 as discussed in the paper.


set.seed(478)
pcut <- 0.15
checkHT2 <- matrix(NA,5000,4)
for (f in 1:5000){
  samp801 <- rbinom(119,1,0.8)
  samp802 <- rbinom(119,1,0.8)
  samp803 <- rbinom(119,1,0.8)
  
  samp701 <- rbinom(119,1,0.7)
  samp702 <- rbinom(119,1,0.7)
  samp703 <- rbinom(119,1,0.7)
  
  # Simple (all within green zone)
  checkHT2[f,1] <- ifelse(mean(samp801)>=0.80 & mean(samp802) >= 0.80 & mean(samp803) >= 0.80,1,0)
  checkHT2[f,2] <- ifelse(mean(samp701)>=0.80 & mean(samp702) >= 0.80 & mean(samp703) >= 0.80,1,0)
  
  # Expanded
  # This calculates if either the estimate is in the green zone,
  # or the p-value is < alpha level, while the estimate is not in RED zone
  # Both endpoints need to be Go for us to proceed
  
  NotinREDs1 <- ifelse(mean(samp801) > 0.70,1,0)
  p801 <- ifelse(binom.test(sum(samp801), 119, p = 0.70, 
                            alternative = "greater")$p.value < pcut,1,0)
  Amber_minor1 <- ifelse(NotinREDs1 + p801 >1,1,0)
  GO801 <- ifelse((ifelse(mean(samp801)>=0.80,1,0)+Amber_minor1)>0 ,1,0)
  
  NotinREDs2 <- ifelse(mean(samp802) > 0.70,1,0)
  p802 <- ifelse(binom.test(sum(samp802), 119, p = 0.70, 
                            alternative = "greater")$p.value < pcut,1,0)
  Amber_minor2 <- ifelse(NotinREDs2 + p802 >1,1,0)
  GO802 <- ifelse((ifelse(mean(samp802)>=0.80,1,0)+Amber_minor2)>0 ,1,0)
  
  NotinREDs3 <- ifelse(mean(samp803) > 0.70,1,0)
  p803 <- ifelse(binom.test(sum(samp803), 119, p = 0.70, 
                            alternative = "greater")$p.value < pcut,1,0)
  Amber_minor3 <- ifelse(NotinREDs3 + p803 >1,1,0)
  GO803 <- ifelse((ifelse(mean(samp803)>=0.80,1,0)+Amber_minor3)>0 ,1,0)
  
  checkHT2[f,3] <- ifelse(GO801 + GO802 + GO803==3,1,0)   
  
  NotinRED701 <- ifelse(mean(samp701) > 0.70,1,0)
  p701 <- ifelse(binom.test(sum(samp701), 119, p = 0.70, 
                            alternative = "greater")$p.value < pcut,1,0)
  Amber_minor701 <- ifelse(NotinRED701 + p701 >1,1,0)
  GO701 <- ifelse((ifelse(mean(samp701)>=0.80,1,0)+Amber_minor701)>0 ,1,0)
  
  NotinRED702 <- ifelse(mean(samp702) > 0.70,1,0)
  p702 <- ifelse(binom.test(sum(samp702), 119, p = 0.70, 
                            alternative = "greater")$p.value < pcut,1,0)
  Amber_minor702 <- ifelse(NotinRED702 + p702 >1,1,0)
  GO702 <- ifelse((ifelse(mean(samp702)>=0.80,1,0)+Amber_minor702)>0 ,1,0)
  
  NotinRED703 <- ifelse(mean(samp703) > 0.70,1,0)
  p703 <- ifelse(binom.test(sum(samp703), 119, p = 0.70, 
                            alternative = "greater")$p.value < pcut,1,0)
  Amber_minor703 <- ifelse(NotinRED703 + p703 >1,1,0)
  GO703 <- ifelse((ifelse(mean(samp703)>=0.80,1,0)+Amber_minor703)>0 ,1,0)
  
  checkHT2[f,4] <- ifelse(GO701 + GO702 + GO703==3,1,0)   
  
}
colMeans(checkHT2)

# p-value threshold of 0.15
# For simple, we get the threshold approach as expected, for
# Expanded we get the tolerance
# approach values (within simulation error) of 0.8178  and 0.01 



# **Optimize Joint Feasibility Space ---------------------------
# Need to optimize the cut point probability
# (Note: This code takes ~ 15 minutes to run)
#a <- Sys.time()
feas_ret <- c(0.70, 0.80)
set.seed(632)
FAsimopt3 <- matrix(NA,10000,2)
n <- 40
for (cond in 1:2){
  ret <- feas_ret[cond]
  for (j in 1:5000){
    # Generate binary retention data
    retentionE <- rbinom(n,1,ret)
    retentionE2 <- rbinom(n,1,ret)
    retentionE3 <- rbinom(n,1,ret)
    
    # Bootstrap
    # Retention and rate are defined for Figure 1, which is the
    # same example, this checks for exceeding the JFS boundary
    boots2o <- matrix(NA,1000,4)
    for (i in 1:1000){
      
      resamprows <-sample(1:n,n,replace = TRUE) 
      boots2o[i,1] <- ifelse(mean(retentionE[resamprows])>=0.8,1,0) # recruitment rate
      boots2o[i,2] <- ifelse(mean(retentionE2[resamprows])>=0.8,1,0) # recruitment rate
      boots2o[i,3] <- ifelse(mean(retentionE3[resamprows])>=0.8,1,0) # recruitment rate
      boots2o[i,4] <- ifelse(boots2o[i,1] + boots2o[i,2] + boots2o[i,3] == 3,1,0) # recruitment rate
      
    }
    
    
    empty_rows <- which(rowSums(is.na(FAsimopt3)) == ncol(FAsimopt3))
    empt <- empty_rows[1]
    
    FAsimopt3[empt,1] <- cond
    FAsimopt3[empt,2] <- mean(boots2o[,4])
    
    
    print(empt)
  } # Simulations loop
} # Condition loop
#Sys.time()-a

# Cut point to discriminate
breaklist <- c(seq(0,1, 0.001))
FAsimopt_S <- FAsimopt3[FAsimopt3[,1] == 2,]
FAsimopt_SN <- FAsimopt3[FAsimopt3[,1] == 1,]

# Check:
picker <- matrix(NA,1001,3)
picker[,1] <- seq(0,1,0.001)
for (l in 1:1001){
  picker[l,2] <- length(FAsimopt_S[FAsimopt_S[,2]>=picker[l,1],2])/length(FAsimopt_S[,2])
  picker[l,3] <- length(FAsimopt_SN[FAsimopt_SN[,2]>=picker[l,1],2])/length(FAsimopt_SN[,2])
  
  
}

cutpoint <- picker[tail(which(picker[,2]>=0.80), n = 1),1]
picker[tail(which(picker[,2]>=0.80), n = 1),2]
picker[tail(which(picker[,2]>=0.80), n = 1),3]

# Choose such that we have P_{G}(SN) = 0.09 to match the HT
picker[20:25,] # 0.022


# Again check and see which cutpoint are used in the table

# For N = 40, gives 
# Gives cut point of 0.035
# P_FA(S) = 0.8018 , P_FA(SN) = 0.0584
# For comparison to tolerance and HT, we want a cut point to give 
# P_{G}(SN) = 0.09, which is 0.022

# Re-run with N = 119 (Table 2 in HT paper) with same seed,
# we get a cut point of 0.012


# Simulation 3 code -----------------------------------
set.seed(473)
nsim1vect <- c(40, 119)
effectv <- c(0.7,0.8)
rhos <- c(0)

out_sim3 <- matrix(NA,(2*2*5000),9)
nv <- ev <- rv <-1

at <- Sys.time()
for (nv in 1:2){
  nsim1 <- nsim1vect[nv]
  
  for (ev in 1:2){
    ret_effect <- effectv[ev]
    
      for (p in 1:5000){
        empty_rows <- which(rowSums(is.na(out_sim3)) == ncol(out_sim3))
        empt <- empty_rows[1]
        # Generate correlated binomial data
        #table(samps)/nsim1
        sams2 <- matrix(NA,nsim1,3)
        sams2[,1] <- rbinom(nsim1,1, ret_effect  )
        sams2[,2] <- rbinom(nsim1,1, ret_effect  )
        sams2[,3] <- rbinom(nsim1,1, ret_effect  )
        
      
        out_sim3[empt,1] <- nsim1
        out_sim3[empt,2] <- ret_effect
        out_sim3[empt,3] <- 0
        
        # Threshold
        e1 <- ifelse(mean(sams2[,1])>=0.8,1,0)
        e2 <- ifelse(mean(sams2[,2])>=0.8,1,0)
        e3 <- ifelse(mean(sams2[,3])>=0.8,1,0)
        out_sim3[empt,4] <- ifelse(e1+e2+e3 == 3,1,0)
        
        
        # Tolerance
        tol <- ifelse(nsim1 == 40, 0.075,0.0521008)
        rounder <- ifelse(nsim1 == 40, 3,7) # to avoid floating point issues
        e1_tol <- ifelse(round(mean(sams2[,1]),rounder)>=round((0.8-tol),rounder),1,0)
        e2_tol <- ifelse(round(mean(sams2[,2]),rounder)>=round((0.8-tol),rounder),1,0)
        e3_tol <- ifelse(round(mean(sams2[,3]),rounder)>=round((0.8-tol),rounder),1,0)
        
        out_sim3[empt,5] <- ifelse(e1_tol+e2_tol+e3_tol==3,1,0) 
        
        # Hypothesis Test
        # Simple (same as threshold)
        out_sim3[empt,6] <- ifelse(e1+e2 +e3 == 3,1,0)
        
        # Expanded
        
        pcut <- ifelse(nsim1 == 40, 1.01, 0.15 )
        NotinREDs1 <- ifelse(mean(sams2[,1]) > 0.70,1,0)
        p801 <- ifelse(binom.test(sum(sams2[,1]), nsim1, p = 0.70, 
                                  alternative = "greater")$p.value < pcut,1,0)
        Amber_minor1 <- ifelse(NotinREDs1 + p801 ==2,1,0)
        GO801 <- ifelse((ifelse(mean(sams2[,1])>=0.80,1,0)+Amber_minor1)>0 ,1,0)
        
        NotinREDs2 <- ifelse(mean(sams2[,2]) > 0.70,1,0)
        p802 <- ifelse(binom.test(sum(sams2[,2]), nsim1, p = 0.70, 
                                  alternative = "greater")$p.value < pcut,1,0)
        Amber_minor2 <- ifelse(NotinREDs2 + p802 ==2,1,0)
        GO802 <- ifelse((ifelse(mean(sams2[,2])>=0.80,1,0)+Amber_minor2)>0 ,1,0)
        
        NotinREDs3 <- ifelse(mean(sams2[,3]) > 0.70,1,0)
        p803 <- ifelse(binom.test(sum(sams2[,3]), nsim1, p = 0.70, 
                                  alternative = "greater")$p.value < pcut,1,0)
        Amber_minor3 <- ifelse(NotinREDs3 + p803 ==2,1,0)
        GO803 <- ifelse((ifelse(mean(sams2[,3])>=0.80,1,0)+Amber_minor3)>0 ,1,0)
        
        # Expanded
        out_sim3[empt,7] <- ifelse(GO801 + GO802 + GO803 == 3,1,0)   
        
        
        # JFS 
        
        boot3 <- matrix(NA,1000,1)
        for (l in 1:1000){
          rowsamp <- sample(1:nsim1, nsim1, replace = TRUE)
          boot3[l,1] <- ifelse(mean(sams2[rowsamp,1]) >= 0.8 & 
                                 mean(sams2[rowsamp,2]) >= 0.8 & 
                                 mean(sams2[rowsamp,3]) >= 0.8,1,0)
        }
        
        # This is the cut point for the comparison to the HT
        disc1 <- ifelse(nsim1 == 40,0.022, 0.012)
        # A different cut point
        disc2 <- ifelse(nsim1 == 40, 0.047, 0.004)
        
        # JFS 80
        out_sim3[empt,8] <- ifelse( mean(boot3) >= disc1,1,0)   
        
        # JFS matched to HT
        out_sim3[empt,9] <- ifelse( mean(boot3) >= disc2,1,0)   
        
      } # p loop
  } # effect loop
} #n loop
Sys.time()-at


out_sim3work <- data.frame(out_sim3)
colnames(out_sim3work) <- c("N", "effect", "rho", "Threshold", "Tolerance", "HTS", 
                            "HTE", "FA", "FAfair")

out_sim3work %>% dplyr::group_by(N, effect) %>% 
  dplyr::summarise(across(everything(), list(mean))) -> tablesim3

# Efficiency comparison simulations -------------
# The code takes ~ 2 days to run for each # of endpoints, in total ~ 6 days.

# *One Endpoint Estimates -----------------------
# Code uses JFS_N to estimate sample size for each S and SN combination
# in Table 3 and then checks with JFScut

rm(list = setdiff(ls(), lsf.str())) # Remove everything except functions
library(dplyr)

dat <- read.csv("G:\\Papers\\JFS\\Submissions\\Submitted code\\Exact_gpower_done.csv", header = TRUE)
HT_1 <- matrix(NA,48,6)
for (r in 1:48){
  
  HT_1[r,1] <- dat[r,1] # SN
  HT_1[r,2] <- dat[r,2] # S
  HT_1[r,3] <- dat[r,3] # N for HT, column 3 for 1, column 5 for 2, column 7 for 3
  
  jfs <-  JFS_N(n_max = 200, 
                PpS = 0.8, 
                PpSN = 0.05,
                S = rep(dat[r,2]*0.01,1),
                SN = rep(dat[r,1]*0.01,1),
                maxit = 5,
                paral = TRUE, 
                ncor = 4, 
                nsimdata = 10000,
                bootrep = 5000)
  if (length(jfs)>1){
    HT_1[r,4] <- as.numeric(jfs$N)
    HT_1[r,5] <- as.numeric(jfs$PpS)
    HT_1[r,6] <- as.numeric(jfs$PpSN)
    
  } else{
    HT_1[r,4] <- NA
    HT_1[r,5] <- NA
    HT_1[r,6] <- NA
  }
  
  print(r)
}

# Save the estimated data set
# write.csv(HT_1,"G:\\Papers\\JFS\\Submissions\\Submitted code\\Estimate_one.csv" )

# *Checking one endpoint -------------------------

endpoints_1 <- read.csv("G:\\Papers\\JFS\\Submissions\\Submitted code\\Estimate_one.csv", header = TRUE)

end1_check <- matrix(NA,48,3)
for (l in 1:48){
  SNv <- 0.01*endpoints_1[l,2]
  Sv <- 0.01*endpoints_1[l,3]
  Nv <- endpoints_1[l,5]
  
  jfs <- JFScut(S = Sv ,
                SN = SNv ,
                nsims = 10000,
                N = Nv, 
                nboot = 10000, 
                PpS = 0.80, 
                parallel = TRUE, ncores = 4)
  
  end1_check[l,1] <- Nv
  end1_check[l,2] <- jfs[[1]]
  end1_check[l,3] <- jfs[[2]]
  print(l)
}

# Any values that did not have P(P|SN) < 0.05 were re-ran iteratively (+1 N)
# until P(P|SN) <= 0.05

# Save the checked data set
# write.csv(end1_check,"G:\\Papers\\JFS\\Submissions\\Submitted code\\endpoint1_check.csv") 


# *Two endpoint estimates --------------------------

rm(list = setdiff(ls(), lsf.str())) # Remove everything except functions
library(dplyr)

dat <- read.csv("G:\\Papers\\JFS\\Submissions\\Submitted code\\Exact_gpower_done.csv", header = TRUE)
a <- Sys.time()
HT_2 <- matrix(NA,48,7)
for (r in 1:48){
  
  HT_2[r,1] <- dat[r,1] # SN
  HT_2[r,2] <- dat[r,2] # S
  HT_2[r,3] <- dat[r,5] # N for HT, column 3 for 1, column 5 for 2, column 7 for 3
  
  jfs <-  JFS_N(n_max = 300, 
                PpS = 0.8, 
                PpSN = 0.05,
                S = rep(dat[r,2]*0.01,2),
                SN = rep(dat[r,1]*0.01,2),
                maxit = 5,
                paral = TRUE, 
                ncor = 4, 
                nsimdata = 10000,
                bootrep = 5000)
  if (length(jfs)>1){
    HT_2[r,4] <- as.numeric(jfs$N)
    HT_2[r,5] <- as.numeric(jfs$PpS)
    HT_2[r,6] <- as.numeric(jfs$PpSN)
    
  } else{
    HT_2[r,4] <- NA
    HT_2[r,5] <- NA
    HT_2[r,6] <- NA
  }
  
  print(r)
}
Sys.time()-a

# Save the estimated data set
#write.csv(HT_2,"G:\\Papers\\JFS\\Submissions\\Submitted code\\Estimate_two.csv" )


# *Checking two endpoints ------------------
end2_check <- matrix(NA,48,3)
endpoints_2 <- read.csv("G:\\Papers\\JFS\\Submissions\\Submitted code\\Estimate_two.csv", header = TRUE)
for (m in 1:48){
  
  SN <- rep(endpoints_2[m,2]*0.01,2)
  S <- rep(endpoints_2[m,3]*0.01,2)
  N <- ceiling(endpoints_2[m,5]*1.0)
  out <- JFScut(S,
                SN,
                N,
                PpS = 0.80,
                parallel = TRUE,
                ncores = 4, 
                nsims = 10000,
                nboot = 10000)
  
  end2_check[m,1] <- N
  end2_check[m,2] <- out[[2]]
  end2_check[m,3] <- out[[3]]
  print(m)
}

# Save the checked data set
#write.csv(end2_check,"G:\\Papers\\JFS\\Submissions\\Submitted code\\endpoint2_check.csv") 

# Any values  that did not have P(P|SN) < 0.05 were re-ran iteratively (+1 N)
# until P(P|SN) <= 0.05


# *Three endpoint estimates ------------

rm(list = setdiff(ls(), lsf.str())) # Remove everything except functions 
library(dplyr)
dat <- read.csv("G:\\Papers\\JFS\\JFS\\Exact_gpower_done.csv", header = TRUE)
a <- Sys.time()
HT_3 <- matrix(NA,48,7)
for (r in 15:48){
  
  HT_3[r,1] <- dat[r,1] # SN
  HT_3[r,2] <- dat[r,2] # S
  HT_3[r,3] <- dat[r,7] # N for HT, column 3 for 1 endpoint, column 5 for 2, column 7 for 3
  
  jfs <-  JFS_N(n_max = HT_3[r,3]*2, 
                PpS = 0.8, 
                PpSN = 0.05,
                S = rep(dat[r,2]*0.01,3),
                SN = rep(dat[r,1]*0.01,3),
                maxit = 5,
                paral = TRUE, 
                ncor = 4, 
                nsimdata = 10000,
                bootrep = 5000)
  
  if (length(jfs)>1){
    HT_3[r,4] <- as.numeric(jfs[[1]])
    HT_3[r,5] <- as.numeric(jfs[[2]])
    HT_3[r,6] <- as.numeric(jfs[[3]])
    
  } else{
    HT_3[r,4] <- NA
    HT_3[r,5] <- NA
    HT_3[r,6] <- NA
  }
  
  print(r)
}
Sys.time()-a

# Save the estimated data set
write.csv(HT_3,"G:\\Papers\\JFS\\Submissions\\Submitted code\\Estimate_three.csv" )



# *Checking three endpoints ------------------
end3_check <- matrix(NA,48,3)
endpoints_3 <- read.csv("G:\\Papers\\JFS\\Submissions\\Submitted code\\Estimate_three.csv", header = TRUE)
for (m in 1:48){
  
  SN <- rep(endpoints_3[m,2]*0.01,3)
  S <- rep(endpoints_3[m,3]*0.01,3)
  N <- ceiling(endpoints_3[m,5]*1.0)
  out <- JFScut(S, 
                SN, 
                N,
                PpS = 0.80,
                parallel = TRUE,
                ncores = 4, 
                nsims = 10000, 
                nboot = 10000)
  
  end3_check[m,1] <- N
  end3_check[m,2] <- out[[2]]
  end3_check[m,3] <- out[[3]]
  print(m)
}

write.csv(end3_check,"G:\\Papers\\JFS\\Submissions\\Submitted code\\endpoint3_check.csv") 

# Any values that did not have P(P|SN) < 0.05 were re-ran iteratively (+1 N)
# until P(P|SN) <= 0.05

