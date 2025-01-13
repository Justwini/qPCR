####
# Is the expression of this tumour suppressor affected by this drug? - Part 2 - Q3
####
## read the data
# make sure you specify the full path to it or that CSV file is in your working directory
d = read.csv("II.csv")

## calcylate dCTs
dCt_treated = d$treated_p27 - d$treated_GAPDH
dCt_untreated = d$untreated_p27 - d$untreated_GAPDH
mean_dCt_untreated = mean(dCt_untreated)

## calcylate ddCTs
ddCt_treated = dCt_treated - mean_dCt_untreated
ddCt_untreated = dCt_untreated - mean_dCt_untreated

## calculate expressions
treated_p27_expression = 2**-ddCt_treated
untreated_p27_expression = 2**-ddCt_untreated

## visualize the data
# plot barplots
img_obj=barplot(c(mean(untreated_p27_expression),
                  mean(treated_p27_expression)),
                names.arg = c("Cnt","Trt"),
                main = "Expression of p27")
# plot the error bars
error.bar(img_obj, 
          c(mean(untreated_p27_expression),mean(treated_p27_expression)),
          c(sd(untreated_p27_expression),sd(treated_p27_expression)))
## run the statistical analysis

t.test(treated_p27_expression,untreated_p27_expression)

#Is the expression of p27 different between un-treated control and treated cells?
#Yes p27 appears to be up-regulated with probability of null hypothesis
#that there is no difference in expression being below the 0.05 thresjold