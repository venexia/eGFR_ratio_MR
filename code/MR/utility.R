# Function to make display numbers ---------------------------------------------

display <- function(x) {
  ifelse(
    abs(x) >= 100,
    sprintf("%.0f", x),
    ifelse(
      abs(x) >= 10,
      sprintf("%.1f", x),
      ifelse(abs(x) >= 0.01, sprintf("%.2f", x), sprintf("%.3f", x))
    )
  )
}
