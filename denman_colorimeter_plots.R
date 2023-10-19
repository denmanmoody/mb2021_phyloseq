# Test script to learn to read in data and run simple plot

# Set working directory (where you pull files from)
setwd("~/Desktop")

# How to open info on a command (argument) = ? 
?read.delim()

# Open data file into environment
data.df <- read.delim(file = "ssm_491_data_collection_2023.txt", header = TRUE, sep = "\t")

# First 6 rows of data (head), Last 6 rows of data (tail), structure of data (str), 
head(data.df)
tail(data.df)
str(data.df)

# How many from the "family" column are unique, What is the length of the family column
unique(data.df$family)
length(data.df$family)

# plot percent fouling by variables
plot(data.df$percent.fouling, data.df$L)
plot(data.df$percent.fouling, data.df$L, col = as.factor(data.df$family))
plot(data.df$percent.fouling, data.df$a)
plot(data.df$percent.fouling, data.df$b)

?plot()

# subsetting variables

# boxplot

# adding legend

# arguments within plot
# xlab
# ylab
# main

# Plot with title (main), x axis title (xlab), y axis title (ylab)
plot(data.df$percent.fouling, data.df$L, 
     main = "Relationship Between Luminosity and Percent Fouling in Yesso Scallops", 
     xlab = "Percent Fouling (%)", ylab = "Luminosity")

# Break title (main) into two lines by using \n (new line)
plot(data.df$percent.fouling, data.df$L, 
     main = "Relationship Between Luminosity and\nPercent Fouling in Yesso Scallops", 
     xlab = "Percent Fouling (%)", ylab = "Luminosity")

# Plot of L vs Percent Fouling with titles, coloured by family
plot(data.df$percent.fouling, data.df$L, 
     main = "Relationship Between Luminosity and\nPercent Fouling in Yesso Scallops", 
     xlab = "Percent Fouling (%)", ylab = "Luminosity", col = as.factor(data.df$family))

# Plot of L (y-axis) vs Percent Fouling (x-axis) with titles, coloured by fouling
plot(data.df$percent.fouling, data.df$L, 
     main = "Relationship Between Luminosity and\nPercent Fouling in Yesso Scallops", 
     xlab = "Percent Fouling (%)", ylab = "Luminosity", col = as.factor(data.df$percent.fouling))

# Plot of L (x-axis) vs Percent Fouling (y-axis) with titles, coloured by fouling, with legend separated into bins of 10% fouling

par(xpd = TRUE)
par(mar = c(5, 4.5, 4, 6) + 0.1)
x <- c(1, 2, 3, 4)
y <- c(10, 15, 12, 8)
percent_fouling <- data.df$percent.fouling
breakpoints <- seq(0, 100, 10)
percent_fouling_bins <- cut(percent_fouling, breakpoints, include.lowest = TRUE, labels = FALSE)

plot(data.df$L, data.df$percent.fouling,
     main = "Relationship Between Luminosity and\nPercent Fouling in Yesso Scallops", 
     ylab = "Percent Fouling (%)", xlab = "Luminosity", col = percent_fouling_bins)

legend_labels <- paste0(breakpoints[-length(breakpoints)], "-", breakpoints[-1], "%")

x_legend <- max(data.df$L) + 2.0
y_legend <- max(data.df$percent.fouling) - 1.9

legend(x_legend, y_legend, legend = legend_labels, col = 1:length(legend_labels), pch = 1,
       title = "Percent Fouling", cex = 0.8)

par(xpd = FALSE)


# Plot of Green to Red (a) (x-axis) vs Percent Fouling (y-axis) with titles, coloured by fouling, with legend separated into bins of 10% fouling
par(xpd = TRUE)
par(mar = c(5, 4.5, 4, 6) + 0.1)
x <- c(1, 2, 3, 4)
y <- c(10, 15, 12, 8)
percent_fouling <- data.df$percent.fouling
breakpoints <- seq(0, 100, 10)
percent_fouling_bins <- cut(percent_fouling, breakpoints, include.lowest = TRUE, labels = FALSE)

plot(data.df$a, data.df$percent.fouling,
     main = "Relationship Between Green to Red (a) and\nPercent Fouling in Yesso Scallops", 
     ylab = "Percent Fouling (%)", xlab = "Green to Red (a)", col = percent_fouling_bins)

legend_labels <- paste0(breakpoints[-length(breakpoints)], "-", breakpoints[-1], "%")

x_legend <- max(data.df$a) + 0.5
y_legend <- max(data.df$percent.fouling) - 1.9

legend(x_legend, y_legend, legend = legend_labels, col = 1:length(legend_labels), pch = 1,
       title = "Percent Fouling", cex = 0.8)

par(xpd = FALSE)







# Plot of Green to Red (a) (x-axis) vs Percent Fouling (y-axis) with titles, coloured by fouling, with legend separated into bins of 10% fouling
### CURRENTLY EDITING
par(xpd = TRUE)
par(mar = c(5, 4.5, 4, 6) + 0.1)
x <- c(1, 2, 3, 4)
y <- c(10, 15, 12, 8)
percent_fouling <- data.df$percent.fouling
breakpoints <- seq(0, 100, 10)
percent_fouling_bins <- cut(percent_fouling, breakpoints, include.lowest = TRUE, labels = FALSE)

plot(data.df$a, data.df$b,
     main = "Relationship Between Green to Red (a) and\nBlue to Yellow (b) in Yesso Scallops", 
     ylab = "Blue to Yellow (b)", xlab = "Green to Red (a)", col = percent_fouling_bins)

legend_labels <- paste0(breakpoints[-length(breakpoints)], "-", breakpoints[-1], "%")

x_legend <- max(data.df$a) + 0.5
y_legend <- max(data.df$b) - 0.5

legend(x_legend, y_legend, legend = legend_labels, col = 1:length(legend_labels), pch = 1,
       title = "Percent Fouling", cex = 0.8)

par(xpd = FALSE)










# Test boxplot
boxplot(data.df$percent.fouling, data.df$L, 
     main = "Relationship Between Luminosity and\nPercent Fouling in Yesso Scallops", 
     xlab = "Percent Fouling (%)", ylab = "Luminosity", col = as.factor(data.df$family))



# Percent fouling of all scallops box plot
boxplot(data.df$percent.fouling,
        main = "Percent Fouling of Yesso Scallops",
        ylab = "Percent Fouling")

# Percent fouling across families box plot
boxplot(percent.fouling ~ family, data = data.df,
        main = "Comparison of Percent Fouling across Families",
        xlab = "Family", ylab = "Percent Fouling",
        las = 2, cex.axis = 0.8)


# Set the bottom margin to create space for the x-axis label
par(mar = c(8, 4, 4, 4) + 0.1)

# Create the box plot with adjusted x-axis label position
boxplot(percent.fouling ~ family, data = data.df,
        main = "Comparison of Percent Fouling across Families",
        xlab = "", ylab = "Percent Fouling",
        las = 2, cex.axis = 0.8)

# Add the x-axis label manually with adjusted position
mtext("Family", side = 1, line = 6)


