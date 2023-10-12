## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(clinfun)

## -----------------------------------------------------------------------------
# Specify the parameters and constraints
trial = ph2simon(0.05, 0.2, 0.1, 0.1)

# Print
trial

## -----------------------------------------------------------------------------
# Plot
plot(trial)

## -----------------------------------------------------------------------------
# Obtain admissible design
twostage.admissible(trial)

## -----------------------------------------------------------------------------
# Specify the parameters and constraints
trial = ph2simon(0.25, 0.45, 0.05, 0.1)

# Obtain admissible design
twostage.admissible(trial)

## -----------------------------------------------------------------------------
# Specify the parameters and constraints & print
ph2single(0.4, 0.6, 0.1, 0.1)

## -----------------------------------------------------------------------------
fe.ssize(p1 = 0.2, p2 = 0.3, power = 0.8)

## -----------------------------------------------------------------------------
power.prop.test(p1 = 0.2, p2 = 0.3, power = 0.8)

