# Libraries for R script
library(tidyverse)
library(GGally)

# load the data
d <- read.csv("Drug_Trial.csv", 
                      #na.strings=c(""," ","NA"), 
                      stringsAsFactors = FALSE, 
                      encoding = "UTF-8", header = TRUE)

count(d) # 1307 rows of data
sum(is.na(d)) # count for missing values

# check for duplicate data, 0 in this instance
d %>%
  count(PATIENT_ID) %>%
  filter(n > 1)

# We have 1307 rows (patients) and 9 columns
dim(d)

# To see column names
dimnames(d)[[2]]

# See the data types of each column
sapply(d, class)
# or
str(d)

# changing the data types as R assumes that some data types are continuous
# and/or not categorical
gender <- as.factor(d[, "GENDER"])

condition_severity <- as.factor(d[, "CONDITION_SEVERITY"])

drug <- as.factor(d[, "DRUG"])

succ_treat <- as.factor(d[, "SUCCESSFUL_TREATMENT"])

# how many males and females
t <- table(gender)
addmargins(t)
round(prop.table(t), digits = 2)

# range of the values
summary(d)

# see the Age distribution
hist(d$AGE)

d$SUCCESSFUL_TREATMENT = factor(d$SUCCESSFUL_TREATMENT)
d$CONDITION_SEVERITY = factor(d$CONDITION_SEVERITY)
d$GENDER = factor(d$GENDER)
d$DRUG = factor(d$DRUG)
d$AGE = as.numeric(d$AGE)
#d$PATIENT_ID = character(d$PATIENT_ID)

ggcorr(d,
       nbreaks = 6,
       label = TRUE,
       label_size = 3,
       color = "grey50")

# See how the successful treatment looks like for our patients
library(ggplot2)
ggplot(d, aes(x = d$SUCCESSFUL_TREATMENT)) +
  geom_bar(width=0.5, fill = "coral") +
  geom_text(stat='count', aes(label=stat(count)), vjust=-0.5) +
  theme_classic()

# women are 50/50, men not very successful
ggplot(d, aes(x = d$GENDER, fill=d$SUCCESSFUL_TREATMENT)) +
  geom_bar(position = position_dodge()) +
  geom_text(stat='count', 
            aes(label=stat(count)), 
            position = position_dodge(width=1), vjust=-0.5)+
  theme_classic()

# how successful a drug is, in terms of treatment (maybe better in excel)
ggplot(d, aes(x = d$DRUG, fill=d$SUCCESSFUL_TREATMENT)) +
  geom_bar(position = position_dodge()) +
  geom_text(stat='count', 
            aes(label=stat(count)), 
            position = position_dodge(width=1), 
            vjust=-0.5)+
  theme_classic()

# how successful a drug is, in terms of treatment (maybe better in excel)
ggplot(d, aes(x = d$CONDITION_SEVERITY, fill=d$SUCCESSFUL_TREATMENT)) +
  geom_bar(position = position_dodge()) +
  geom_text(stat='count', 
            aes(label=stat(count)), 
            position = position_dodge(width=1), 
            vjust=-0.5)+
  theme_classic()

# Discretize age to plot survival
d$Discretized.age = cut(d$AGE, c(0,10,20,30,40,50,60,70,80,100))
# Plot discretized age
ggplot(d, aes(x = d$Discretized.age, fill=d$SUCCESSFUL_TREATMENT)) +
  geom_bar(position = position_dodge()) +
  geom_text(stat='count', aes(label=stat(count)), position = position_dodge(width=1), vjust=-0.5)+
  theme_classic()
d$Discretized.age = NULL

# age density plot
ggplot(d, aes(x = d$AGE)) +
  geom_density(fill='coral')

# condition severity plots (done in excel for speed)


# cross-tabulations
age <- d[,'AGE']
missed_doses <- d[, 'MISSED_DOSES']
dosage <- d[, 'DOSAGE']
missed_checks <- d[, 'MISSED_CHECKS']

age_grouped <- ifelse(age < 45, 'under 45', age_grouped)
age_grouped <- ifelse(age >= 45 & age < 65, '45 - 64', age_grouped)
age_grouped <- ifelse(age >= 65 & age < 75, '65 - 74', age_grouped)
age_grouped <- ifelse(age >= 75, '75 or over', age_grouped)


table(age_grouped, exclude = NULL)

age_group_by_gender <- table(age_grouped, gender, exclude = NULL)

age_group_by_gender

round(100 * prop.table(age_group_by_gender, margin = 2), digits = 1)

head(cbind(age_grouped, age))

### filter out outlier patients, do a boxplot of their age
boxplot(d$AGE~d$GENDER,data=d, main="Age vs Gender of Patients",
        xlab="Gender", ylab="AGE")

boxplot(d$AGE~d$CONDITION_SEVERITY,data=d, main="Age vs Condition Severity of Patients",
        xlab="Condition Severity", ylab="AGE")


### Chi squared tests (out of scope)


### Logistic regression 
# outcome variable
# predictors
# distribution for the outcome variable, binomial
# link function logit
# Null model with outcome variable and 1 intercept
m <- glm(succ_treat ~ 1, family = binomial (link = logit))
summary(m)
table(m$y) # 0: 966, 1: 341 35%
exp(-1.04128) # and then divided by 1 plus the odds, x100 you get the odds of having a successful treatment!

# Model with 1 independent/predictor variable
m <- glm(succ_treat ~ gender, family = binomial (link = logit))
summary(m)

contrasts(gender)
levels(gender)

m$coefficients
exp(m$coefficients)
# logistic regression for age
m <- glm(succ_treat ~ age, family = binomial (link = logit))
summary(m) # log odds of having a successful treatment decreases with age but age not significant

# successful treatment and age cross-tabulation
st_by_age <- table(age, succ_treat)

# output the frequencies of st by age
freq_table <- prop.table(st_by_age, margin = 1)

# calculate the odds of getting successful treatment
odds <- freq_table[, "1"]/freq_table[, "0"]

# calculate the log odds
logodds <- log(odds)

# plot the ages found in the sample against the log odds of getting successful treatment
plot(rownames(freq_table), logodds)

#### Multiple logistic regression

# Correlation test to see if missed_doses correlates with missed_checks
cor.test(x=d$MISSED_DOSES, y=d$MISSED_CHECKS, method='pearson')
cor.test(x=d$MISSED_DOSES, y=d$DOSAGE, method='pearson')
cor.test(x=d$DOSAGE, y=d$MISSED_CHECKS, method='pearson')

# model 
m <- glm(succ_treat ~ age + gender + condition_severity + drug + missed_checks + missed_doses + dosage,
         family=binomial (link=logit))
summary(m)

exp(cbind(coef(m)))
exp(confint(m))

# model again without drug variable
m <- glm(succ_treat ~ age + gender + condition_severity + missed_checks + missed_doses + dosage,
         family=binomial (link=logit))
summary(m)
