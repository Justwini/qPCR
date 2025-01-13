####
# Is the expression of this tumour suppressor affected by this drug? - Part 2 - Q4
####
## read the data
# make sure you specify the full path to it or that CSV file is in your working directory
d = read.csv("IV.csv")

## calcylate dCTs
dCt_treated = d$treated_p27 - d$treated_GAPDH
dCt_untreated = d$untreated_p27 - d$untreated_GAPDH
mean_dCt_untreated = mean(dCt_untreated)
dCt_negativeControl = d$negativeControl_p27 - d$negativeControl_GAPDH
dCt_positiveControl = d$positiveControl_p27 - d$positiveControl_GAPDH

## calculate ddCTs
ddCt_treated = dCt_treated - mean_dCt_untreated
ddCt_untreated = dCt_untreated - mean_dCt_untreated
ddCt_negativeControl = dCt_negativeControl - mean_dCt_untreated
ddCt_positiveControl = dCt_positiveControl - mean_dCt_untreated

## calculate expressions
treated_p27_expression = 2**-ddCt_treated
untreated_p27_expression = 2**-ddCt_untreated
negativeControl_p27_expression = 2**-ddCt_negativeControl
positiveControl_p27_expression = 2**-ddCt_positiveControl

## visualize the data
# plot barplots
img_obj=barplot(c(mean(untreated_p27_expression),
                  mean(negativeControl_p27_expression),
                  mean(treated_p27_expression),
                  mean(positiveControl_p27_expression)),
                names.arg = c("UT","NC","Treated","PC"),
                main = "p27 expression")
# plot the error bars
error.bar = function(x, y, upper, lower=upper, length=0.1,...) { 
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper)) stop("vectors must be same length")     
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)}

error.bar(img_obj, 
          c(mean(untreated_p27_expression),mean(negativeControl_p27_expression),mean(treated_p27_expression),mean(positiveControl_p27_expression)),
          c(sd(untreated_p27_expression),sd(negativeControl_p27_expression),sd(treated_p27_expression),sd(positiveControl_p27_expression)))

## run the statistical analysis
# prepare the data for ANOVA and pairwise t tests
expression = c(treated_p27_expression,
               untreated_p27_expression,
               negativeControl_p27_expression,
               positiveControl_p27_expression)
condition = c(rep("treated_p27_expression",10),
              rep("untreated_p27_expression",10),
              rep("negativeControl_p27_expression",10),
              rep("positiveControl_p27_expression",10))

data=as.data.frame(t(rbind(condition,expression)))
data$expression = as.numeric(data$expression)
# run ANOVA
summary(aov(formula = expression ~ condition, data = data))
# run pairwise t tests and note the P value adjustment method
see=pairwise.t.test(x=expression,g=condition,p.adjust.method = "fdr")
see$p.value

#Which of the following best describes the outcome from this experiments assuming that by "difference" we mean statistically significant alteration in expression?

#Although it appears that the experiment was overall successful as judged by enhanced
#expression of p27 in positive control and no difference between un-treated cells and
#negative control, as well as degree of up-regulation of p27 in treated cells 
#realtive to un-treated cells andnegative control, it is not possible to say that p27 is up-regulated as  a result of drug treatment.
#This is  becauyse there is no statistical significant difference in expression between negative control and treated cells