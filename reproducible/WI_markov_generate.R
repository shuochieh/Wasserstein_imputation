library(ggplot2)
library(gridExtra)

q = 3
q. = 4

MAR = T
p1 = 1
p2 = 1
p3 = 1
p4 = 1

p1. = 1
p2. = 1
p3. = 1
p4. = 1

a1_star = 0.6 * p1 / (0.6 * p1 + 0.15 * p2)
a2_star = 0.6 * p1. / (0.6 * p1. + 0.15 * p2.)
b1_star = 0.15 * p3 / (0.15 * p3 + 0.1 * p4)
b2_star = 0.15 * p3. / (0.15 * p3. + 0.1 * p4.)

slope1 = ((0.6 * p1 + 0.15 * p2) * q.) / ((0.6 * p1. + 0.15 * p2.) * q)
itcp1 = - (0.6 * (p1 / q - p1. / q.)) / ((0.6 * p1. + 0.15 * p2.) / q.)
slope2 = ((0.15 * p3 + 0.1 * p4) * q.) / ((0.15 * p3. + 0.1 * p4.) * q)
itcp2 = -(0.15 * (p3 / q - p3. / q.)) / ((0.15 * p3. + 0.1 * p4.) / q.)

x_grid = seq(from = 0, to = 1, length.out = 500)
y_grid = slope1 * x_grid + itcp1

dta = cbind(x_grid, y_grid)
dta[which(((y_grid < 0) | (y_grid > 1))),] = NA

df = as.data.frame(dta)
colnames(df) = c("a1", "a2")

if (MAR) {
  plot1 = ggplot(df, aes(x = a1, y = a2)) +
    geom_line(size = 1.5, color = "blue") +
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), linetype = "dashed", color = "orange", size = 1) +
    geom_point(aes(x = a1_star, y = a2_star), size = 3, color = "red") +
    labs(title = "",
         x = "a1",
         y = "a2") +
    ylim(0, 1) + 
    xlim(0, 1)
} else {
  plot1 = ggplot(df, aes(x = a1, y = a2)) +
    geom_line(size = 1.5, color = "blue") +
    geom_point(aes(x = a1_star, y = a2_star), size = 3, color = "red") +
    labs(title = "",
         x = "a1",
         y = "a2") +
    ylim(0, 1) + 
    xlim(0, 1)
}

x_grid = seq(from = 0, to = 1, length.out = 500)
y_grid = slope2 * x_grid + itcp2

dta2 = cbind(x_grid, y_grid)
dta2[which(((y_grid < 0) | (y_grid > 1))),] = NA

df2 = as.data.frame(dta2)
colnames(df2) = c("b1", "b2")

if (MAR) {
  plot2 = ggplot(df2, aes(x = b1, y = b2)) +
    geom_line(aes(color = "TWI solutions"), size = 1.5) +
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), linetype = "dashed", color = "orange", size = 1) +
    geom_point(aes(x = b1_star, y = b2_star, color = "Correct imputation"), size = 3) +
    labs(title = "",
         x = "b1",
         y = "b2", 
         color = "") +
    ylim(0, 1) + 
    xlim(0, 1) +
    scale_color_manual(values = c("TWI solutions" = "blue", "Correct imputation" = "red")) +
    guides(color = guide_legend(override.aes = list(
      linetype = c("solid", "blank"),
      shape = c(NA, 16))))
} else {
  plot2 = ggplot(df2, aes(x = b1, y = b2)) +
    geom_line(aes(color = "TWI solutions"), size = 1.5) +
    geom_point(aes(x = b1_star, y = b2_star, color = "Correct imputation"), size = 3) +
    labs(title = "",
         x = "b1",
         y = "b2", 
         color = "") +
    ylim(0, 1) + 
    xlim(0, 1) +
    scale_color_manual(values = c("TWI solutions" = "blue", "Correct imputation" = "red")) +
    guides(color = guide_legend(override.aes = list(
      linetype = c("solid", "blank"), 
      shape = c(NA, 16))))
}

grid.arrange(plot1, plot2, ncol = 2, widths = c(1, 1.3))