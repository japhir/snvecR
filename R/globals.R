VER <- "snvecR VERSION: 3.7.5 2023-04-25"

AU <- 1.49597870700e11 # m
GM <- 1.32712440041e20 # m3/s2
OM <- 7.292115e-5 # 1/s EarthRot
R0 <- 3.8440e8 # m Moon R0
GK <- 0.9925194 # Kinoshita75,77
ED0 <- 0.0032738134 # DynEll (C-A)/C
FGP <- 0.99961908
AU3 <- AU * AU * AU
R03 <- R0 * R0 * R0
# no need to define pi in R
R2D <- 180. / pi # radians to deg

D2S <- 3600. * 24.
Y2D <- 365.25
KY2D <- 1.e3 * Y2D

# set default Ed, Td
ED <- 1.0000 # set factor 1.0
TD <- 0.0000 # set factor 0.0

# mass ratios
MSEL <- 328900.5596 # MS/(ME+ML)
MEL <- 81.300568 # ME/ML
MLS <- 1. / (MSEL * (1 + MEL)) # ML/MS
# K0, beta0 for torques
K0 <- (3. / 2.) * GM * ED0 * ED / (OM * AU3)
K0D <- K0 * D2S # 1/s => 1/d */
BET0 <- GK * MLS * AU3 / R03
K0B0 <- K0D * (1. + BET0)
# Moon mean motion
N0 <- sqrt(GM / MSEL / R03)
NW0 <- (N0 / OM) # ratio (n/om)_0
# Tidal dissipation Quinn91 Eqs. (3, 11)
# NDN = (dndt/n)_0, WDW = (domdt/om)_0
NDN <- (-4.6e-18 * D2S * TD) # 1/s => 1/d
WDW <- (51. * NDN * NW0) # Lambeck80
# tidal effect on obliquity
UEPSDOT <- -4.17e-19

# SunRot Angles (Transform to HCI)
OMT <- 75.5940
INCT <- 7.155
EP0 <- 23.439291111111110 # Obliquity t0
