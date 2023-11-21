# E(PFS) helpers:
survPFS <- function(transition, t) {
  UseMethod("survPFS")
}

survPFS.ExponentialTransition <- function(transition, t) {
  ExpSurvPFS(t = t, h01 = transition$hazards$h01, h02 = transition$hazards$h02)
}

survPFS.WeibullTransition <- function(transition, t) {
  WeibSurvPFS(t = t, h01 = transition$hazards$h01, h02 = transition$hazards$h02,
              p01 = transition$weibull_rates$p01, p02 = transition$weibull_rates$p02)
}

survPFS.PWCTransition <- function(transition, t) {
  PWCsurvPFS(t, h01 = transition$hazards$h01, h02 = transition$hazards$h02,
             pw01 = transition$intervals$pw01, pw02 = transition$intervals$pw02)
}

# E(OS) helpers:
survOS <- function(transition, t) {
  UseMethod("survOS")
}

survOS.ExponentialTransition <- function(transition, t) {
  ExpSurvOS(t = t, h01 = transition$hazards$h01, h02 = transition$hazards$h02,
            h12 = transition$hazards$h12)
}

survOS.WeibullTransition <- function(transition, t) {
  WeibSurvOS(t = t, h01 = transition$hazards$h01, h02 = transition$hazards$h02,
             h1 = transition$hazards$h1, p01 = transition$weibull_rates$p01,
             p02 = transition$weibull_rates$p02, p12 = transition$weibull_rates$p12)
}

survOS.PWCTransition <- function(transition, t) {
  PWCsurvOS(t = t, h01 = transition$hazards$h01, h02 = transition$hazards$h02,
            h12 = transition$hazards$h12, pw01 = transition$intervals$pw01,
            pw02 = transition$intervals$pw02, pw12 = transition$intervals$pw12)
}

# Var(PFS)/Var(OS) helpers:
expvalPFSInteg <- function(transition, t) {
  t * survPFS(transition, t)
}

expvalOSInteg <- function(transition, t) {
  t * survOS(transition, t)
}

# E(PFS*OS) helpers:

#------- P_11 (see Meller et al, 4274)

p11Integ <- function(x, transition) {
  haz(transition = transition, t = x, trans = 3)
}

p11 <- function(transition, s, t) {
  intval <- mapply(function(s, t) integrate(p11Integ, lower = s, upper = t, transition)$value, s, t)
  exp(-intval)
}

#------- Surv(PFS*OS)

PFSOSInteg <- function(u, t, transition) {
  p11(transition, u, t/u) * survPFS(transition, u) * haz(transition, u, 1)
}

survPFSOS <- function(t, transition) {
  sapply(t, function(x) {
    intval <- integrate(PFSOSInteg, lower = 0, upper = sqrt(x), x, transition)$value
    survPFS(transition, sqrt(x)) + intval
  })
}

# Main function:
corPFSOS <- function(transition) {
  # E(PFS) & E(OS)
  expvalPFS <- integrate(survPFS, lower = 0, upper = Inf,
                         transition = transition)$value

  expvalOS <- integrate(survOS, lower = 0, upper = Inf,
                        transition = transition)$value

  # Var(PFS) & Var(OS)
  expvalPFS2 <- 2 * integrate(expvalPFSInteg, lower = 0, upper = Inf,
                          transition = transition)$value

  expvalOS2 <- 2 * integrate(expvalOSInteg, lower = 0, upper = Inf,
                         transition = transition)$value

  varPFS <- expvalPFS2 - expvalPFS^2

  varOS <- expvalOS2 - expvalOS^2

  # E(PFS*OS)
  expvalPFSOS <- integrate(survPFSOS, lower = 0, upper = Inf,
                           transition)$value

  # Cor(PFS, OS)
  (expvalPFSOS - expvalPFS*expvalOS) / (varPFS * varOS)^0.5
}
