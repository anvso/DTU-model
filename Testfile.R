pdf("font_plot.pdf", family="Impact", width=13.2, height=4)
plot(mtcars$mpg, mtcars$wt, 
     main = "Fuel Efficiency of 32 Cars",
     xlab = "Weight (x1000 lb)",
     ylab = "Miles per Gallon")
dev.off()

?pdf