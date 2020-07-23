require(wiqid)
require(car)
require(lmtest)
require(ggpubr)
require(lme4)
require(lmerTest)
require(blme)
require(MuMIn)

data <- read.table("data.csv",sep=";", header = T)

################ Latitude and Elevation: MAYFLIES ################

data_mayf <- subset(data, order=="mayflies") # Subset mayflies data

data_mayf$residual_maxvel <- residuals(lm(log(max.vel+1) ~ log(body_length_mm),na.action=na.exclude, data_mayf))
data_mayf$residual_avgvel <- residuals(lm(log(avg.vel+1) ~ log(body_length_mm),na.action=na.exclude, data_mayf))
plot(log(max.vel+1) ~ log(body_length_mm), data_mayf, pch=20, ylab="log(velocity)")
points(log(avg.vel+1) ~ log(body_length_mm), data_mayf)
abline(lm(log(max.vel+1) ~ log(body_length_mm), data_mayf))
abline(lm(log(avg.vel+1) ~ log(body_length_mm), data_mayf), lty=3)
data_mayf$residual_maxvel <- scale(data_mayf$residual_maxvel)
data_mayf$residual_avgvel <- scale(data_mayf$residual_avgvel)
data_mayf$temp_scale <- scale(data_mayf$temp_c)

### Low elevation (mayflies)
mod_1_Low <- lmer(residual_maxvel ~ temp_scale + lat + temp_scale:lat + (1|ID), data_mayf[data_mayf$elev=="Low",]) # linear model
mod_2_Low <- lmer(residual_maxvel ~ temp_scale + I(temp_scale^2) + lat + temp_scale:lat + (1|ID), data_mayf[data_mayf$elev=="Low",]) # quadratic model
AICc(mod_1_Low, mod_2_Low)
anova(mod_2_Low)

# Split temperate and tropical data (only low elevation)
data_sub_low <- subset(data_mayf, elev == "Low")
mod_1_N_low <- lmer(residual_maxvel ~ temp_scale + (1|ID), data_sub_low[data_sub_low$lat=="N",]) # linear
mod_2_N_low <- lmer(residual_maxvel ~ temp_scale + I(temp_scale^2) + (1|ID), data_sub_low[data_sub_low$lat=="N",]) #quadratic
AICc(mod_1_N_low, mod_2_N_low) 
summary(mod_2_N_low)

mod_1_S_low <- lmer(residual_maxvel ~ temp_scale + (1|ID), data_sub_low[data_sub_low$lat=="S",]) # Note: Singularity occurs because the estmated variance of individual ID is zero. Compare these results with those obtained with blmer (which computes variance ID using Bayesian approach) -> results are consistent. 
mod_2_S_low <- lmer(residual_maxvel ~ temp_scale + I(temp_scale^2) + (1|ID), data_sub_low[data_sub_low$lat=="S",])
AICc(mod_1_S_low, mod_2_S_low)
summary(mod_2_S_low)

### Mid elevation (mayflies)
mod_1_Mid <- lmer(residual_maxvel ~ temp_scale + lat + temp_scale:lat + (1|ID), data_mayf[data_mayf$elev=="Mid",])
mod_2_Mid <- lmer(residual_maxvel ~ temp_scale + I(temp_scale^2) + lat + temp_scale:lat + (1|ID) , data_mayf[data_mayf$elev=="Mid",])
AICc(mod_1_Mid, mod_2_Mid)
summary(mod_1_Mid)

# Split temperate and tropical data (only mid elevation)
data_sub_mid <- subset(data_mayf, elev == "Mid")
mod_1_N_mid <- lmer(residual_maxvel ~ temp_scale + (1|ID), data_sub_mid[data_sub_mid$lat=="N",]) # linear
mod_2_N_mid <- lmer(residual_maxvel ~ temp_scale + I(temp_scale^2) + (1|ID), data_sub_mid[data_sub_mid$lat=="N",]) #quadratic
AICc(mod_1_N_mid, mod_2_N_mid) 
summary(mod_1_N_mid)

mod_1_S_mid <- lmer(residual_maxvel ~ temp_scale + (1|ID), data_sub_mid[data_sub_mid$lat=="S",])
mod_2_S_mid <- lmer(residual_maxvel ~ temp_scale + I(temp_scale^2) + (1|ID), data_sub_mid[data_sub_mid$lat=="S",])
AICc(mod_1_S_mid, mod_2_S_mid)
summary(mod_1_S_mid) 

###  High elevation (mayflies)
mod_1_High <- lmer(residual_maxvel ~ temp_scale + lat + temp_scale:lat + (1|ID), data_mayf[data_mayf$elev=="High",])
mod_2_High <- lmer(residual_maxvel ~ temp_scale + I(temp_scale^2) + lat + temp_scale:lat + (1|ID), data_mayf[data_mayf$elev=="High",])
AICc(mod_1_High, mod_2_High) 
anova(mod_2_High) 

# Split temperate and tropical data (only high elevation)
data_sub_high <- subset(data_mayf, elev == "High")
mod_1_N_high <- lmer(residual_maxvel ~ temp_scale + (1|ID), data_sub_high[data_sub_high$lat=="N",]) # linear
mod_2_N_high <- lmer(residual_maxvel ~ temp_scale + I(temp_scale^2) + (1|ID), data_sub_high[data_sub_high$lat=="N",]) #quadratic
AICc(mod_1_N_high, mod_2_N_high) 
summary(mod_2_N_high)

mod_1_S_high <- lmer(residual_maxvel ~ temp_scale + (1|ID), data_sub_high[data_sub_high$lat=="S",])
mod_2_S_high <- lmer(residual_maxvel ~ temp_scale + I(temp_scale^2) + (1|ID), data_sub_high[data_sub_high$lat=="S",])
AICc(mod_1_S_high, mod_2_S_high)
summary(mod_2_S_high)

################ Latitude and Elevation: STONEFLIES ################ 

data_stonef <- subset(data, order=="stoneflies") # Subset stoneflies

data_stonef$residual_maxvel <- residuals(lm(log(max.vel+1) ~ log(body_length_mm),na.action=na.exclude, data_stonef))
data_stonef$residual_avgvel <- residuals(lm(log(avg.vel+1) ~ log(body_length_mm),na.action=na.exclude, data_stonef))
plot(log(max.vel+1) ~ log(body_length_mm), data_stonef, pch=20)
points(log(avg.vel+1) ~ log(body_length_mm), data_stonef)
abline(lm(log(max.vel+1) ~ log(body_length_mm), data_stonef))
abline(lm(log(avg.vel+1) ~ log(body_length_mm), data_stonef), lty=3)
data_stonef$residual_maxvel <- scale(data_stonef$residual_maxvel)
data_stonef$residual_avgvel <- scale(data_stonef$residual_avgvel)
data_stonef$temp_scale <- scale(data_stonef$temp_c)

### Low elevation (stoneflies)
mod_1_Low <- lmer(residual_maxvel ~ temp_scale + lat + temp_scale:lat + (1|ID), data_stonef[data_stonef$elev=="Low",]) # linear model
mod_2_Low <- lmer(residual_maxvel ~ temp_scale + I(temp_scale^2) + lat + temp_scale:lat + (1|ID), data_stonef[data_stonef$elev=="Low",]) # quadratic model
AICc(mod_1_Low, mod_2_Low)
summary(mod_1_Low) 

## Split N and S data (only low elevation)
data_sub_low <- subset(data_stonef, elev == "Low")
mod_1_N <- lmer(residual_maxvel ~ temp_c + (1|ID), data_sub_low[data_sub_low$lat=="N",])
mod_2_N <- lmer(residual_maxvel ~ temp_c + I(temp_c^2) + (1|ID), data_sub_low[data_sub_low$lat=="N",])
AICc(mod_1_N, mod_2_N)
summary(mod_1_N)

mod_1_S <- lmer(residual_maxvel ~ temp_c + (1|ID), data_sub_low[data_sub_low$lat=="S",])
mod_2_S <- lmer(residual_maxvel ~ temp_c + I(temp_c^2) + (1|ID), data_sub_low[data_sub_low$lat=="S",])
AICc(mod_1_S, mod_2_S) 
summary(mod_1_S) 

### Mid elevation (stoneflies)
mod_1_Mid <- lmer(residual_maxvel ~ temp_scale + lat + temp_scale:lat + (1|ID), data_stonef[data_stonef$elev=="Mid",])
mod_2_Mid <- lmer(residual_maxvel ~ temp_scale + I(temp_scale^2) + lat + temp_scale:lat + (1|ID), data_stonef[data_stonef$elev=="Mid",])
AICc(mod_1_Mid, mod_2_Mid)
summary(mod_1_Mid) 

## Split N and S data (only mid elevation)
data_sub_mid <- subset(data_stonef, elev == "Mid")
mod_1_N <- lmer(residual_maxvel ~ temp_c + (1|ID), data_sub_mid[data_sub_mid$lat=="N",])
mod_2_N <- lmer(residual_maxvel ~ temp_c + I(temp_c^2) + (1|ID), data_sub_mid[data_sub_mid$lat=="N",])
AICc(mod_1_N, mod_2_N)
summary(mod_1_N)

mod_1_S <- lmer(residual_maxvel ~ temp_c + (1|ID), data_sub_mid[data_sub_mid$lat=="S",])
mod_2_S <- lmer(residual_maxvel ~ temp_c + I(temp_c^2) + (1|ID), data_sub_mid[data_sub_mid$lat=="S",])
AICc(mod_1_S, mod_2_S) 
summary(mod_1_S) 

################ MAYFLIES vs STONEFLIES ################ 

data_mayf <- subset(data, order=="mayflies") # Subset mayflies data
data_mayf$residual_maxvel <- residuals(lm(log(max.vel+1) ~ log(body_length_mm), na.action=na.exclude, data_mayf))
data_mayf$residual_maxvel <- scale(data_mayf$residual_maxvel)
data_mayf$temp_scale <- scale(data_mayf$temp_c)

data_stonef <- subset(data, order=="stoneflies") # Subset stoneflies
data_stonef$residual_maxvel <- residuals(lm(log(max.vel+1) ~ log(body_length_mm),na.action=na.exclude, data_stonef))
data_stonef$residual_maxvel <- scale(data_stonef$residual_maxvel)
data_stonef$temp_scale <- scale(data_stonef$temp_c)

data_total <- rbind(data_mayf, data_stonef)

data_mid <- subset(data_total, elev == "Mid")
data_mid_N <- subset(data_mid, lat == "N")
data_mid_S <- subset(data_mid, lat == "S")

data_low <- subset(data_total, elev == "Low")
data_low_N <- subset(data_low, lat == "N")
data_low_S <- subset(data_low, lat == "S")

# Temperate, mid elevations
mod_mid_N_1 <- lmer(residual_maxvel ~ temp_scale + order + temp_scale:order + (1|ID), data_mid_N)
mod_mid_N_2 <- lmer(residual_maxvel ~ temp_scale + I(temp_scale^2) + order + temp_scale:order + (1|ID), data_mid_N)
AICc(mod_mid_N_1, mod_mid_N_2)
anova(mod_mid_N_1)

# Tropical, mid elevations
mod_mid_S_1 <- lmer(residual_maxvel ~ temp_scale + order + temp_scale:order + (1|ID), data_mid_S)
mod_mid_S_2 <- lmer(residual_maxvel ~ temp_scale + I(temp_scale^2) + order + temp_scale:order + (1|ID), data_mid_S)
AICc(mod_mid_S_1, mod_mid_S_2)
anova(mod_mid_S_1)

# Temperate, low elevations
mod_low_N_1 <- lmer(residual_maxvel ~ temp_scale + order + temp_scale:order + (1|ID), data_low_N)
mod_low_N_2 <- lmer(residual_maxvel ~ temp_scale + I(temp_scale^2) + order + temp_scale:order + (1|ID), data_low_N)
AICc(mod_low_N_1, mod_low_N_2)
anova(mod_low_N_2)

# Tropical, low elevations
mod_low_S_1 <- lmer(residual_maxvel ~ temp_scale + order + temp_scale:order + (1|ID), data_low_S)
mod_low_S_2 <- lmer(residual_maxvel ~ temp_scale + I(temp_scale^2) + order + temp_scale:order + (1|ID), data_low_S)
AICc(mod_low_S_1, mod_low_S_2)
anova(mod_low_S_2)













