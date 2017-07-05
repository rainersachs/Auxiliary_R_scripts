#Cost-free, open-source software under the (lenient) GNU GPLv3 license. Comes with no warraty. Concerns radiogenic mouse HG tumorigenesis
#Written by Mark Ebert, Edward Huang, and Rainer Sachs, summer 2017.
#Relevant references and their abbreviations in commenting this script are the following.
#".93Alp"=Alpen et al. "Tumorigenic potential of high-Z, high-LET charged-particle radiations." Rad Res 136:382-391 (1993)
#".94Alp"=Alpen et al. "Fluence-based relative biological effectiveness for charged particle carcinogenesis in mouse Harderian gland." Adv Space Res 14(10): 573-581. (1994).  
#"16Chang"=Chang et al. "Harderian Gland Tumorigenesis: Low-Dose and LET Response." Radiat Res 185(5): 449-460. (2016).  
#"16Srn"=Siranart et al."Mixed Beam Murine Harderian Gland Tumorigenesis: Predicted Dose-Effect Relationships if neither Synergism nor Antagonism Occurs." Radiat Res 186(6): 577-591 (2016).  
#"17Cuc"=Cucinotta & Cacao. "Non-Targeted Effects Models Predict Significantly Higher Mars Mission Cancer Risk than Targeted Effects Models." Sci Rep 7(1): 1832. (2017). PMC5431989
rm(list=ls())
ddfr=data.frame( #Includes all HZE data used in 16Chang  
  Dose=c(0.2,0.4,0.6,1.2,2.4,3.2,5.1,7,0.05,0.1,0.15,0.2,0.4,0.8,1.6,0.05,0.1,0.2,0.4,0,0.1,0.2,0.4,0.8,1.6,0.4,0.8,1.6,3.2,0.05,0.1,0.2,0.4,0.8,0.1,0.2,0.4,0.8,0.1,0.2,0.4,0.04,0.08,0.16,0.32,0.033,0.066,0.13,0.26,0.52,.2, .4, .6),
  HG=c(0.091,0.045,0.101,0.169,0.347,0.431,0.667,0.623,0.156,0.215,0.232,0.307,0.325,0.554,0.649,0.123,0.145,0.207,0.31,0.026,0.083,0.25,0.39,0.438,0.424,0.093,0.195,0.302,0.292,0.109,0.054,0.066,0.128,0.286,0.183,0.167,0.396,0.536,0.192,0.234,0.317,0.092,0.131,0.124,0.297,0.082,0.088,0.146,0.236,0.371,.154,.132,.333),#HG prevalence as defined in 16Chang
  NWeight=c(520,2048,1145,584,313,232,293,221,1162,877,455,409,374,223,320,742,661,347,131,6081,1091,251,244,191,131,645,255,199,111,649,378,973,833,201,468,381,197,109,496,257,185,1902,1063,884,350,1767,1408,874,299,261,322,206,67),#nominal weight for weighted least squares regression; see .93Alp. The Lanthanum entries were obtained by measuring the main graph in 17Cuc 
  index=c(rep(1,8),rep(0,17), rep(1,4),  rep(0,24)),#index=0 for Z>3 ions, index=1 for proton p (=1H1) and 2HE4 ions
  L=c(rep(1.6,8),rep(193,7),rep(250,4),rep(195,6),rep(0.4,4),rep(25,5),rep(464,4),rep(193,3),rep(70,4),rep(100,5),rep(953,3)), #L=LET=LET_infinity=stopping power
  Z=c(rep(2,8),rep(26,17),rep(1,4),rep(10,5),rep(43,4),rep(26,3),rep(14,4),rep(22,5),rep(57,3)),#proton #, e.g. 2 for 2He4
  Zeff=c(rep("TBD",53)),# effective ion charge according to the formula of W.H Barkas. Zeff <= Z. Calculated below
  beta=c(rep("TBD",53)),# ion speed, relative to speed of light, calculated below
  MeVperu=c(rep(228,8),rep(600,7),rep(300,4),rep(600,6),rep(250,4),rep(670,5),rep(600,4),rep(600,3),rep(260,4),rep(1000,5),rep(593,3)),#Kinetic energy in MeV, divided by atomic mass, e.g. divided by 4u=4x931.5 MeV/c^2 for 2He4
  Katz=c(rep("TBD",53)),#for fully ionized nuclei, Katz's Z^2/beta^2, Calculated below. It is part of the Bethe Barkas Bloch equation for stopping power. Our calculations don't use Katz, but various similar calculations do.
  ion=c(rep("He4",8),rep("Fe56",17),rep("p",4),rep("Ne20",5),rep("Nb93",4),rep("Fe56",3),rep("Si28",4),rep("Ti48",5),rep("La139",3)),
  comments=c(".93AlpLooksOK",rep("",7),".93AlplooksOK",rep('',11),".93Alp.no.iso", "not in 17Cuc (or 16Chang?)",rep("",3),"16Chang all OK?",rep('',24),".94Alp","From graphs",'e.g. in 17Cuc')) 
Y=0.001*ddfr[,"MeVperu"]# convert to GeV/u for convenience in a calculation
ddfr[,"Katz"]=round((ddfr[,"Z"])^2*(2.57*Y^2+4.781*Y+2.233)/(2.57*Y^2+4.781*Y),2)#special relativistic calculation of Z^2/beta^2. The numerics include conversion from GeV to joules and from u to kg.
ddfr[,"beta"]=round(ddfr[,"Z"]*sqrt(1/ddfr[,"Katz"]),3)#i.e. Z*sqrt(beta^2/Z^2)
ddfr[,"Zeff"]=round(ddfr[,"Z"]*(1-exp(-125*ddfr[,"Z"]^(-2.0/3))),2)#Barkas formula for Zeff; for us Zeff is almost Z

### Data for HG induced by photons from Cs-137 or Co-60 beta decay; from 16Chang; needed
ddd=data.frame(
  Dose=c(0, 0.4, 0.8, 1.6, 3.2, 7, 0, .4, .8, .12, 1.6),
  HG=c(.026, .048, .093, .137, .322, .462, .0497, .054, .067, .128, .202),
  NWeight=c(6081.2, 4989.5, 1896.8, 981.1, 522.2, 205.2, 7474.1, 2877.6, 1423.7, 689.9, 514.9),
  Nucleus=c(rep("Cobalt-60", 6), rep("Cesium-137", 5)),
  Comments=c(rep("TBD", 11))
  )
LL=ddd[,"Dose"]
HG=ddd[,"HG"]
WT=ddd[,"NWeight"]
QQ=LL^2
LQ=lm(HG~LL+QQ,weights=WT)
summary(LQ, correlation=T)

### Chunk to enable visual tests of ddfr
aa=20;bb=25# checking ddfr for individual ions against graphs in .93Alp, 16Chang, and 17Cuc
plot(ddfr[aa:bb,"Dose"],ddfr[aa:bb,"HG"], ann='F') #example for checking ddfr, in this case Fe56 at 600 MeV/u data
####
ddfra=ddfr[c(1:19,26:53),] ##removes the zero dose case and the no isograft data

######New model. When calibrating from data could replace (1-exp(-150*phi*Dose/L) by 1 at every observed dose point !=0
phi=3e3#controls how fast NTE build up from zero; not really needed since 150*phi*Dose/L>>1 at every observed Dose !=0  
IDERm=nls(HG~1-exp(-0.01*(2.75+(1-index)*aa1*L*Dose*exp(-aa2*L) +index*bb*Dose^2*exp(-ll0*Dose)+
                            (1-exp(-150*phi*Dose/L))*(1-index)*kk1)),data=ddfra, weights=NWeight,
  start=list(aa1=.9,aa2=.01,bb=4.5,ll0=.2, kk1=0.048))#calibrating parameters; only need 5 parameters for very good results
summary(IDERm,correlation=T)# parameter values & accuracy. All 5 parameters have p-values <<10^(-3)!
cc=coef(IDERm)#calibrated central values of the 5 parameters
HH=function(Dose,L,index){#The Hazard function; Model will be 1-exp(-HH); see 17Cuc
  0.01*(2.75+(1-index)*cc[1]*L*Dose*exp(-cc[2]*L) +index*cc[3]*Dose^2*exp(-cc[4]*Dose)+
           (1-exp(-150*phi*Dose/L))*(1-index)*cc[5])
}
HGc=function(Dose,L,index) 1-exp(-HH(Dose,L,index))##this is the calibrated model. 

#######Next: visual checks to make sure our calibration is consistent with 16Chang, .93Alp, .94Alp and 17Cuc
## Put various values in our calibrated new model to check with numbers and graphs in these references, e.g.
L=1.6; Dose=ddfr[1:8,"Dose"];HGe=ddfr[1:8,"HG"];index=1#He in .93Alp. HGe=experimental HG.
#L=.4; Dose=ddfr[26:29,"Dose"];HGe=ddfr[26:29,"HG"];index=1#same for protons. 
#L=193;Dose=ddfr[9:15,"Dose"]; HGe=ddfr[9:15,"HG"];index=0#same for Fe
plot(c(0,7),c(0,1),col='red', ann='F')# biggest Dose-effect range ever considered
lines(Dose,HGc(Dose,L,index))
points(Dose,HGe)#looks great for Helium, OK for protons; very good for Iron
