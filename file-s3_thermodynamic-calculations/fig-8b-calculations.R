# Copyright (C) 2024 Derek Leung <derekdvleung@gmail.com>
# This program is free software: you can redistribute it and/or modify it under the terms 
# of the GNU General Public License version 3 as published by the Free Software Foundation. 
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
# See the GNU General Public License for more details. You should have received a copy of 
# the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

# This program calculates the reaction path assemblages involved in the hydrothermal alteration and
# carbonation of ultramafic rocks (see Leung et al. 2024), but considers fluid-mobile parameters in 
# the quartz-saturated proximal alteration sequence.  

#adds the CHNOSZ library for thermodynamic modelling
library (CHNOSZ)

#resets the workspace
reset()
getwd()

###################################################################
# set T, P, and ranges for log fCO2 and log aSiO2

T <- 300
#325 for Bohlke (1989); 300 for Kishida and Kerrich (1987)
P <- 2600
#2000 for Bohlke (1989); 2000-3200 for Hagemann and Brown (1996), Channer and Spooner (1991)

###################################################################
# pre-define fluid-mobile constraints

Na <- -1
K <- -2
CO2 <- c(1,3, 600)
pH <- c(3,7)

###################################################################
# display settings - use names = FALSE for a diagram that is labeled post-export for publication
names = FALSE

###################################################################
# import minerals from thermoddem database

#switch to kelvin for T
T.units("K")

dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
#setwd("..")
add.OBIGT("thermoddem-mineral-data.csv")

#swich temperature back to degrees celsius
T.units("C")

###################################################################
# Main space for calculating activity diagrams

# define the minerals in MSHC compositional space  
Mg.cr <- c("Magnesite(Synth)", "Antigorite", "Lizardite", "Chrysotile", "Talc", "Brucite",  "Cummingtonite")

# define the minerals in NKMASHC compositional space  
Al.cr <- c("Clinochlore",
           "Amesite", 
           "Hydrotalcite(CO3)","Spinel", 
           "Muscovite(ordered)", 
           "Albite(low)", 
           "Diaspore",
           "Clinochlore95",
           "Clinochlore90",
           "Clinochlore85",
           "Clinochlore80",
           "Clinochlore75",
           "Clinochlore70",
           "Clinochlore65",
           "Clinochlore60",	
           "Clinochlore55",
           "Clinochlore50",
           "Celadonite", 
           "Kaolinite",
           "Pyrophyllite",
           "Corundum(alpha)",
           "Corundum(gamma)",
           "Phlogopite"
)

# colours to fill the NKMASHC fields
fillColors <- c(
  "#48506130",
  "#A9A66960",
  "#FF000010",
  "#FF00FF20",
  "#99FFE330",
  "#FFF5CE60",
  NA,
  "#51586144",
  "#5B616259",
  "#6569636E",
  "#6E726482",
  "#787B6597",
  "#828365AC",
  "#8B8C66C0",
  "#959467D5",
  "#9F9D68EA",
  "#A9A669FF",
  NA,  
  NA,
  NA,
  NA,
  NA,
  NA,
  NA,
  NA,
  NA,
  NA,
  NA,  
  NA,
  NA,
  NA,
  NA,
  NA,
  NA,
  NA,
  NA,
  NA
)

# define the minerals in CMSHC compositional space    
Ca.cr <- c("Dolomite(ordered)","Tremolite", "Diopside", "Calcite")

# define the minerals in NKFASHC compositional space (experimental)    
#Fe.cr <- c("Magnetite", "Siderite", "Greenalite", "Minnesotaite", "Fe(OH)2")
#"Hematite"

# define the minerals in CMASHC compositional space (experimental)   
#CaAl.cr <- c("Clinozoisite", "Edenite(alpha)", "Pyroxene(CaAl)")


#sets the compositional space to NCKMASHC
#Mg is based on conservation of Mg in magnesite
#Al is based on conservation of Al in muscovite
#Ca is based on conservation of Ca in dolomite
#Si is fixed by the activity of quartz

basis (c("Magnesite(Synth)", "Quartz(alpha)", "H2O", "CO2", "O2", "Muscovite(ordered)", "H+", "Dolomite(ordered)", "Na+","K+"))
basis ("CO2", "gas")
basis ("Na+", Na)
basis ("K+", K)
#basis ("Ca+2", Ca)
#basis ("O2", -37)

###########################
# plot the species in MSHC and overlay NKMASHC over the base MSHC

species (c(Al.cr))

# correct for the activity of spinel because pure spinel is metastable with respect to diaspore
# (Mg0.52Fe2+0.44Mn0.03Zn0.01)Σ1.00(Cr1.13Al0.79Fe3+0.06V0.01)Σ1.99O4; Menzel et al. (2018)
species ("Spinel", -1.09)

mAl <- mosaic (Mg.cr, "CO2" = CO2, "pH" = pH, T = T, P = P)
diagram(mAl$A.bases, add = FALSE, col = "#b4b4b4", names = names, col.names = "#b4b4b4", lty =1, italic = TRUE, lwd = 2, fill = NA)
dAl <- diagram(mAl$A.species, add = TRUE, lwd = 3, fill = fillColors, names = names)

# Manually remove the y-axis
axis(side = 2, labels = FALSE, tick = FALSE)

###########################
# overlay CMSHC over MSHC      
#species (c(Ca.cr))
#mCa <- mosaic (Mg.cr, "CO2" = CO2, "pH" = pH, T = T, P = P)
#dCa <- diagram(mCa$A.species, add = TRUE, col = "#aa0000", names = names, col.names = "#aa0000", lwd = 2)

#species (c(CaAl.cr))
#mCaAl <- mosaic (Mg.cr, "CO2" = CO2, "pH" = pH, T = T, P = P)
#dCaAl <- diagram(mCaAl$A.species, add = TRUE, col = 2, col.names = 2, bold = TRUE)
#a11 <- mix (dCa, dAl, dCaAl, c(1,1))

#species (c(Fe.cr))
#species ("Siderite", -0.7)
#mFe <- affinity ("CO2" = CO2, "pH" = pH, T = T, P = P)
#dFe <- diagram (mFe, add = TRUE)

# Add legend and title

TP <- describe.property(c("T", "P", "pH"), c(T, P, 6))
#SP <- describe.basis(c(8,9))
#legend1 <- lex(TP, SP)
legend1 <- lex(TP)
legend("bottomleft", legend1, bg = "white", cex = 1.0,y.intersp = 1.5, text.width = 0.28)

