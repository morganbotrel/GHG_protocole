#Auteur: Richard LaBrie

#DATA est le fichier de donn?es qui contient les ID, ECD, TCD
#Temp?rature de l'eau, Temp?rature d'analyse, Salinit? et pression atmosph?rique

#Headspace est la quantit? d'espace dans la wheaton pour faire l'?quilibre

#VolumeWheaton est le volume de la Wheaton, g?n?ralement 0,06L

#N2Oair et CO2air est la quantit? de gaz pour faire le headspace. 0 ici car on utilise du N2 pur

#UnitePression est l'unit? de pression lors de l'?chantillonnage. Par d?faut ? kPa, mais on peut mettre en mmHg

#N2Ocal=c(0.3,1,3) et CO2cal = c(200,1500,10000) sont les concentrations en N2O et CO2 de la calibration

#N2Oppmair CO2ppmair sont les valeurs de N2O et CO2 dans l'air sur le site d'?chantillonnage

#Methode est un param?tre logique qui indique quelle m?thode a ?t? utilis?e, Wheaton ou Shaky. Par d?faut, il est ? Shaky

#ATTENTION, IL FAUT ROULER LE SCRIPT DEUX FOIS. UNE FOIS POUR AVOIR LA CONCENTRATION DANS L'AIR
#ET UNE DEUXI?ME FOIS POUR AVOIR LES VALEURS CORRIG?ES DANS L'EAU

#La fonction sort un tableau avec les ID et les valeurs de N2O et CO2 ? saturation, observ?
#d?viation de la saturation (ppm et nM) et supersaturation (%)

GasWheaton <- function(DATA, Headspace = 0.012, VolumeWheaton = 0.06, N2Oair = 0, CO2air = 0, methode="Shaky",  UnitePression = "kPa", N2Ocal=c(0.3,1,3), CO2cal = c(200,1500,10000), N2Oppmair=1, CO2ppmair=1,air=F)
{
  #attach(DATA)
  if(UnitePression== "mmHg")
  {
    Pression = Pression / 7.500615613026438    #COnversion mmHg en pKa, source: www.convertunits.com    
  }
  
  #Graphique et r?gression pour la calibration du N2O
  
  plot(as.numeric(c(N2Ocal,N2Ocal)) ~ ECD[c(1,2,3,length(ECD)-5,length(ECD)-4,length(ECD)-3)],
       xlab = "Valeur ECD", ylab = "N2O ppm", main="Calibration N2O")
  N2Ocalib = lm(c(N2Ocal,N2Ocal) ~ ECD[c(1,2,3, length(ECD)-5,length(ECD)-4,length(ECD)-3)])
  abline(N2Ocalib)
  N2Opente = N2Ocalib$coefficients[2]
  N2Ointercept = N2Ocalib$coefficients[1]
  N2ORcarree = summary(N2Ocalib)$r.squared
  legend(x="topleft", paste("N2O = ", round(N2Opente,3), "* ECD + ", round(N2Ointercept,3),"R? = ", round(N2ORcarree,4)))
  
  #Graphique et r?gression quadratique pour la calibration du CO2
  
  plot(as.numeric(c(CO2cal,CO2cal)) ~ TCD[c(4,5,6,length(TCD)-2,length(TCD)-1,length(TCD))],
       xlab = "Valeur TCD", ylab = "CO2 ppm", main="Calibration CO2",xlim=c(0,TCD[6]))
  CO2calib = lm(c(CO2cal,CO2cal) ~ poly(TCD[c(4,5,6, length(TCD)-2,length(TCD)-1,length(TCD))],2,raw=T))
  CO2carree = CO2calib$coefficients[3]
  CO2pente = CO2calib$coefficients[2]
  CO2intercept = CO2calib$coefficients[1]
  par(new=T)
  curve(CO2calib$coefficients[3]*x^2 + CO2calib$coefficients[2]*x + CO2calib$coefficients[1], xlim=c(0,TCD[6]),axes=F,xlab="",ylab="")
  CO2Rcarree = summary(CO2calib)$r.squared
  legend(x="topleft", paste("CO2 = ", round(CO2carree,3), "* TCD? + ", round(CO2pente,3), "* TCD + ", round(CO2intercept,3),"R? = ", round(CO2Rcarree,4)))
  
  #Calcul des concentration en N2O et CO2 ? partir des calibrations
  N2Oppm = ECD * N2Opente + N2Ointercept
  CO2ppm = TCD^2 * CO2carree + TCD * CO2pente + CO2intercept
  if(air==T)
  {
    output = as.data.frame(cbind(as.character(ID),N2Oppm,CO2ppm))
    return(output)
  }
  
  #Variables transitoires pour le calcul des concentrations corrigees
  Ratio = Headspace / (VolumeWheaton - Headspace)
  VolumeMolaire = 0.082057 * (Tb + 273.15) * 101.325/Pression
  SolubiliteN2Ow = exp(-62.7062 + (97.3066*100/(Tw+273.15)) + 24.1406*log((Tw + 273.15)/100) +
                       Salinite * (-0.05842 + 0.033193*(Tw+273.15)/100 - 0.0051313 * ((Tw + 
                       273.15)/100)^2)) #Solubilit? du N2O dans l'eau

  SolubiliteN2Ob = exp(-62.7062 + (97.3066*100/(Tb+273.15)) + 24.1406*log((Tb + 273.15)/100) + 
                       Salinite * (-0.05842 + 0.033193*(Tb+273.15)/100 - 0.0051313 * ((Tb + 
                       273.15)/100)^2)) #Solubilit? du N2O dans la bouteille
  
  SolubiliteCO2w = 10^-(-((9345.17/(Tw + 273.15))-60.2409+23.3585*log((273.15+Tw)/100)+
                    Salinite*(0.023517-0.023656*((273.15+Tw)/100)+0.0047036*
                    ((273.15+Tw)/100)^2))/log(10)) #Solubilit? du CO2 dans l'eau 
  
  SolubiliteCO2b = 10^-(-((9345.17/(Tb + 273.15))-60.2409+23.3585*log((273.15+Tb)/100)+
                   Salinite*(0.023517-0.023656*((273.15+Tb)/100)+0.0047036*
                   ((273.15+Tb)/100)^2))/log(10)) #Solubilit? du CO2 dans la bouteille
  
  #Calcul des concentrations en C2O et CO2 corrig?es
  if(methode == "Shaky")
  {
    N2Oair = N2Oppmair
    CO2air = CO2ppmair
  }
  N2Oppmcor = (N2Oppm * SolubiliteN2Ob + Ratio * (N2Oppm - N2Oair)/VolumeMolaire) / SolubiliteN2Ow
  CO2ppmcor = ((CO2ppm - CO2air)*Ratio/VolumeMolaire + CO2ppm * SolubiliteCO2b) / SolubiliteCO2w
  
  #Calcul de la saturation en N2O et CO2
  N2OsatnM = SolubiliteN2Ow * N2Oppmair * 1000 #Saturation du N2O si ?quilibre avec air
  N2OwnM = SolubiliteN2Ow * N2Oppmcor * 1000
  DeltaN2OnM = N2OwnM - N2OsatnM
  DeltaN2Oppm = DeltaN2OnM / (SolubiliteN2Ow * 1000)
  N2Osupersaturation = DeltaN2OnM / N2OsatnM * 100
  
  CO2satnM = SolubiliteCO2w * CO2ppmair * 1000 #Saturation du CO2 si ?quilibre avec air
  CO2wnM = SolubiliteCO2w * CO2ppmcor * 1000
  DeltaCO2nM = CO2wnM - CO2satnM
  DeltaCO2ppm = DeltaCO2nM / (SolubiliteCO2w * 1000)
  CO2supersaturation = DeltaCO2nM / CO2satnM * 100
  
  Gazconcentration = as.data.frame(cbind(as.character(ID),N2OsatnM,N2OwnM,DeltaN2OnM,DeltaN2Oppm,N2Osupersaturation,CO2satnM,CO2wnM,DeltaCO2nM,DeltaCO2ppm,CO2supersaturation))
  return(Gazconcentration)
}