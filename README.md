# Dynamic-network-effects
Non-parametric estimation of reciprocity and triadic effects in relational event networks

This repository contains R code for modelling the dynamic structure of reciprocal and triadic effects in relational event networks via stratified baseline hazards. This avoids having to define arbitrary temporal windows for what can be considered reciprocity or triadic closure. A two-step estimation framework is used, based on the stratified Cox proportional hazards model, followed by non-parametric estimation of the stratified cumulative hazard functions. The proposed approach is illustrated using cyclic closure, but it can also be used for other triadic closure types such as transitive closure, receiving balance closure, etc.

**The directory structure is as follows:**\
simulation code - R source code for network simulation and modelling\
data preparation - function that incorporates censored events\
strata - function that assigns events to corresponding strata\
likelihood - function for the likelihood estimation\
modelling - functions for the penalised spline smoothness selection and estimation of the baseline hazard

**Contact:** Rūta Juozaitienė (Užupytė) - r.uzupyte@gmail.com
