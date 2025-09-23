library(ggplot2)
library(dplyr)

find_steady_point <- function(y0, a1, t1, x_range = c(1000, 100000), threshold = 0.001) {
  # Define ExpDesc1
  f <- function(x) {
    y0 + a1 *  exp(-x/t1)
  }
  
  # Define prime
  f_prime <- function(x) {
    -(a1/t1) * exp(-x/t1)
  }
  
  # Set x steps
  x_seq <- seq(x_range[1], x_range[2], length.out = 1000000)
  derivatives <- f_prime(x_seq)
  
  # Find steady index
  steady_index <- which(abs(derivatives) < threshold)[1]
  
  if (!is.na(steady_index)) {
    steady_point <- x_seq[steady_index]
    steady_value <- f(steady_point)
    
    return(list(
      steady_point = steady_point,
      steady_value = steady_value,
      derivative = derivatives[steady_index],
      convergence = TRUE
    ))
  } else {
    min_deriv_index <- which.min(abs(derivatives))
    return(list(
      steady_point = x_seq[min_deriv_index],
      steady_value = f(x_seq[min_deriv_index]),
      derivative = derivatives[min_deriv_index],
      convergence = FALSE
    ))
  }
}

# Set y0, a1, t1 from origin software
y0_origin =45.31447 
a1_origin = -35.09572
t1_origin = 13387.37122

#Calculate steady point
result <- find_steady_point(y0 =y0_origin , a1 = a1_origin, t1 = t1_origin, threshold = 0.000001)
cat("Steady x ≈", round(result$steady_point), "\n")
cat("Steady y值:", round(result$steady_value, 2), "\n")
cat("Derivative:", result$derivative, "\n")

#Prepare spots for ExpDesc1 funtion 
my_function <- function(x, y0, a1, t1) {
  y0 + a1 * exp(-x/t1) 
}
plot_data <- data.frame(
  x = seq(1000, 200000, length.out = 1000)
) %>%
  mutate(y = my_function(x, y0 =y0_origin , a1 = a1_origin, t1 = t1_origin))

#Calculate 0.7, 0.8, 0.9 Power points
steady_point0 <- result$steady_point
steady_value0 <- result$steady_value
steady_point <- -t1_origin*log((0.8*y0_origin-y0_origin)/a1_origin)
steady_point2 <- -t1_origin*log((0.7*y0_origin-y0_origin)/a1_origin)
steady_point3 <- -t1_origin*log((0.9*y0_origin-y0_origin)/a1_origin)
steady_value <- result$steady_value*0.8
steady_value2 <- result$steady_value*0.7
steady_value3 <- result$steady_value*0.9


#Plot by Simulation Power y/y0_origin
p <- ggplot(plot_data, aes(x = x, y = y/y0_origin)) +
  geom_line(color = "steelblue", linewidth = 1.2) +
  # Add 80% power point
  geom_vline(xintercept = steady_point, linetype = "dashed", color = "orange", linewidth = 1) +
  geom_hline(yintercept = steady_value/y0_origin, linetype = "dashed", color = "orange", linewidth = 1) +
  geom_point(data = data.frame(x = steady_point, y = steady_value/y0_origin), 
             aes(x = x, y = y), color = "red", size = 3) +
  # Add 70% power point
  geom_vline(xintercept = steady_point2, linetype = "dashed", color = "orange", linewidth = 1) +
  geom_hline(yintercept = steady_value2/y0_origin, linetype = "dashed", color = "orange", linewidth = 1) +
  geom_point(data = data.frame(x = steady_point2, y = steady_value2/y0_origin), 
             aes(x = x, y = y), color = "red", size = 3) +
  # Add 90% power point
  geom_vline(xintercept = steady_point3, linetype = "dashed", color = "orange", linewidth = 1) +
  geom_hline(yintercept = steady_value3/y0_origin, linetype = "dashed", color = "orange", linewidth = 1) +
  geom_point(data = data.frame(x = steady_point3, y = steady_value3/y0_origin), 
             aes(x = x, y = y), color = "red", size = 3) +
  labs(
    title = "Origin simulation line",
    x = "Sample Size (x)",
    y = "Simulation Power (y)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.line = element_line(linewidth = 0.5, arrow = arrow(type = "closed", length = unit(0.2, "cm"))),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
    panel.grid.minor = element_blank()
  ) +
  scale_x_continuous(labels = scales::comma) +
  # Add coordinates
  annotate("text", x = steady_point * 2, y = steady_value * 0.990/y0_origin, 
           label = paste("(", round(steady_point), ",", round(steady_value/y0_origin, 1), ")"), 
           color = "red") +
  annotate("text", x = steady_point2 * 2.5, y = steady_value2 * 0.990/y0_origin, 
           label = paste("(", round(steady_point2), ",", round(steady_value2/y0_origin, 1), ")"), 
           color = "red") +
  annotate("text", x = steady_point3 * 1.8, y = steady_value3 * 0.990/y0_origin, 
           label = paste("(", round(steady_point3), ",", round(steady_value3/y0_origin, 1), ")"), 
           color = "red")

print(p) 
ggsave('D:\\生信\\ctDRTF\\Revision\\Origin_simulation_line.pdf',p)










