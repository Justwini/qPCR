####
# Is the expression of this tumour suppressor affected by this drug? - Part 2 - Q3
####
## read the data
# make sure you specify the full path to it or that CSV file is in your working directory
d = read.csv("I.csv")

## calcylate dCTs
dCt_treated = d$treated_p53 - d$treated_GAPDH
dCt_untreated = d$untreated_p53 - d$untreated_GAPDH
mean_dCt_untreated = mean(dCt_untreated)

## calcylate ddCTs
ddCt_treated = dCt_treated - mean_dCt_untreated
ddCt_untreated = dCt_untreated - mean_dCt_untreated

## calculate expressions
treated_p53_expression = 2**-ddCt_treated
untreated_p53_expression = 2**-ddCt_untreated

## visualize the data
# plot barplots
img_obj=barplot(c(mean(untreated_p53_expression),
                  mean(treated_p53_expression)),
                names.arg = c("Cnt","Trt"),
                main = "Expression of p53")
# plot the error bars
error.bar = function(x, y, upper, lower=upper, length=0.1,...) { 
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper)) stop("vectors must be same length")     
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)}

error.bar(img_obj, 
          c(mean(untreated_p53_expression),mean(treated_p53_expression)),
          c(sd(untreated_p53_expression),sd(treated_p53_expression)))
## run the statistical analysis

#You were excited to see difference in the expression of p53 between treated cells and un-treated control.
#After repeating the experiment, is there still enough evidence to suggest the existence of such a difference?
 
t.test(treated_p53_expression,untreated_p53_expression)

#There is no sufficient evidence to reject null hypothesis stating
#that  expression between un-treated and treated cell population is the same


