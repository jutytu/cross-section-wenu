import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math

# loading the data from the CERN website
url = "http://opendata.cern.ch/record/545/files/Wenu.csv"
wenu = pd.read_csv(url)

# adding iso variables
wenu["iTpT"] = wenu["isoTrack"]/wenu["pt"]
wenu["iEpT"] = wenu["isoEcal"]/wenu["pt"]
wenu["iHpT"] = wenu["isoHcal"]/wenu["pt"]
wenu["mt"] = np.sqrt(2 * wenu["pt"]* wenu["MET"]*(1 - np.cos( wenu['phi'] - wenu['phiMET'])))
print(wenu.describe().to_string())
print(wenu.head().to_string())
# listing the columns (variables)
print("Wenu variables: " + str(list(wenu.columns)))
print("Number of events: " + str(wenu.shape[0]))
print("\n")
# number of accelerator runs from which the data was collected
print(wenu["Run"].value_counts().to_frame())



# histograms

# using cuts on the eta and pt variables to differentiate between regions of the detector
# EB = end barrel, EE = end cap
EB = wenu[(abs(wenu.eta) < 1.44) & ( wenu.pt > 25)]
EE = wenu[(abs(wenu.eta) > 1.57) & (abs(wenu.eta) < 2.5) & (wenu.pt > 25)]

print(EB.describe().to_string())
print(EE.describe().to_string())

figEB = plt.figure(figsize = (10,15))

plt.subplot(3, 2, 2)
plt.hist(EB.eta, bins = 50, range = [-1.5, 1.5])
plt.xlabel("eta", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.ylim(0, 2000)

plt.subplot(3, 2, 1)
plt.hist(EB.pt, bins = 50, range = [20, 100])
plt.xlabel("electron pt [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.ylim(0, 5000)

plt.subplot(3, 2, 3)
plt.hist(EB.Q, bins = 3, range = [-1, 1])
plt.xlabel("charge", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.ylim(0, 30000)
plt.xticks([-1, 0, 1])

plt.subplot(3, 2, 4)
plt.hist(EB.MET, bins = 50, range = [0, 100])
plt.xlabel("MET [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.ylim(0, 4000)
plt.savefig('Fig1.png')
plt.show()


figEE = plt.figure(figsize = (10,15))

plt.subplot(3, 2, 1)
plt.hist(EE.pt, bins = 50, range = [20, 100])
plt.xlabel("electron pt [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.ylim(0, 5000)

plt.subplot(3, 2, 2)
plt.hist(EE.eta, bins = 50, range = [-3, 3])
plt.xlabel("eta", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.ylim(0, 4000)

plt.subplot(3, 2, 3)
plt.hist(EE.Q, bins = 3, range = [-1, 1])
plt.xlabel("charge", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.ylim(0, 30000)
plt.xticks([-1, 0, 1])

plt.subplot(3, 2, 4)
plt.hist(EE.MET, bins = 50, range = [0, 100])
plt.xlabel("MET [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.ylim(0, 4000)
plt.savefig('Fig2.png')
plt.show()

# identification variables: cuts on multiple variables to reduce background (mistakenly recorded events
# from other decays)

cutEB1 = EB[(EB.HoverE < 0.04) & (EB.iTpT < 0.09) & (EB.iEpT < 0.07) & (EB.iHpT < 0.1)]
cutEB2 = EB[(EB.sigmaEtaEta < 0.01) & (EB.iTpT < 0.09) & (EB.iEpT < 0.07) & (EB.iHpT < 0.1)]
cutEB3 = EB[(EB.sigmaEtaEta < 0.01) & (EB.HoverE < 0.04) & (EB.iEpT < 0.07) & (EB.iHpT < 0.1)]
cutEB4 = EB[(EB.sigmaEtaEta < 0.01) & (EB.HoverE < 0.04) & (EB.iTpT < 0.09) & (EB.iHpT < 0.1)]
cutEB5 = EB[(EB.sigmaEtaEta < 0.01) & (EB.HoverE < 0.04) & (EB.iEpT < 0.07) & (EB.iTpT < 0.09)]

cutEE1 = EE[(EE.HoverE < 0.025) & (EE.iTpT < 0.05) & (EE.iEpT < 0.06) & (EE.iHpT < 0.025)]
cutEE2 = EE[(EE.sigmaEtaEta < 0.031) & (EE.iTpT < 0.05) & (EE.iEpT < 0.06) & (EE.iHpT < 0.025)]
cutEE3 = EE[(EE.sigmaEtaEta < 0.031) & (EE.HoverE < 0.025) & (EE.iEpT < 0.06) & (EE.iHpT < 0.025)]
cutEE4 = EE[(EE.sigmaEtaEta < 0.031) & (EE.HoverE < 0.025) & (EE.iTpT < 0.05) & (EE.iHpT < 0.025)]
cutEE5 = EE[(EE.sigmaEtaEta < 0.031) & (EE.HoverE < 0.025) & (EE.iTpT < 0.05) & (EE.iEpT < 0.06)]


# histograms (comparing distributions of every variable before and after performing cuts
# on all the other variables except for the one plotted (EB1/EE1, EB2/EE2...)

figcuts1 = plt.figure(figsize = (10,15))

plt.subplot(3, 2, 1)
plt.hist(EB.sigmaEtaEta, bins = 100, range = (0, 0.025))
plt.xlabel("sigmaEtaEta", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')
plt.ylim(0, 10000)

plt.subplot(3, 2, 1)
plt.hist(cutEB1.sigmaEtaEta, bins = 100, range = (0, 0.025))
plt.xlabel("sigmaEtaEta", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')
plt.ylim(0, 10000)

plt.subplot(3, 2, 2)
plt.hist(EE.sigmaEtaEta, bins = 100, range = (0, 0.055))
plt.xlabel("sigmaEtaEta", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')
plt.ylim(0, 10000)

plt.subplot(3, 2, 2)
plt.hist(cutEE1.sigmaEtaEta, bins = 100, range = (0, 0.055))
plt.xlabel("sigmaEtaEta", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')
plt.ylim(0, 10000)



plt.subplot(3, 2, 3)
plt.hist(EB.HoverE, bins = 100, range = (0, 0.15))
plt.xlabel("HoverE", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')
plt.ylim(0, 1000)

plt.subplot(3, 2, 3)
plt.hist(cutEB2.HoverE, bins = 100, range = (0, 0.15))
plt.xlabel("HoverE", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')
plt.ylim(0, 1000)

plt.subplot(3, 2, 4)
plt.hist(EE.HoverE, bins = 100, range = (0, 0.15))
plt.xlabel("HoverE", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')
plt.ylim(0, 10000)

plt.subplot(3, 2, 4)
plt.hist(cutEE2.HoverE, bins = 100, range = (0, 0.15))
plt.xlabel("HoverE", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')
plt.ylim(0, 10000)
plt.savefig('Fig3.png')
plt.show()


figcuts2 = plt.figure(figsize = (10,15))


plt.subplot(3, 2, 1)
plt.hist(EB.iTpT, bins = 100, range = (0, 0.5))
plt.xlabel("IsoTrack/pt", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')
plt.ylim(0, 10000)

plt.subplot(3, 2, 1)
plt.hist(cutEB3.iTpT, bins = 100, range = (0, 0.5))
plt.xlabel("IsoTrack/pt", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')
plt.ylim(0, 10000)

plt.subplot(3, 2, 2)
plt.hist(EE.iTpT, bins = 100, range = (0, 0.5))
plt.xlabel("IsoTrack/pt", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')
plt.ylim(0, 10000)

plt.subplot(3, 2, 2)
plt.hist(cutEE3.iTpT, bins = 100, range = (0, 0.5))
plt.xlabel("IsoTrack/pt", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')
plt.ylim(0, 10000)



plt.subplot(3, 2, 3)
plt.hist(EB.iEpT, bins = 100, range = (0, 0.5))
plt.xlabel("IsoEcal/pt", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')
plt.ylim(0, 10000)

plt.subplot(3, 2, 3)
plt.hist(cutEB4.iEpT, bins = 100, range = (0, 0.5))
plt.xlabel("IsoEcal/pt", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')
plt.ylim(0, 10000)

plt.subplot(3, 2, 4)
plt.hist(EE.iEpT, bins = 100, range = (0, 0.5))
plt.xlabel("IsoEcal/pt", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')
plt.ylim(0, 10000)

plt.subplot(3, 2, 4)
plt.hist(cutEE4.iEpT, bins = 100, range = (0, 0.5))
plt.xlabel("IsoEcal/pt", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')
plt.ylim(0, 10000)


plt.subplot(3, 2, 5)
plt.hist(EB.iHpT, bins = 100, range = (0, 0.5))
plt.xlabel("IsoHcal/pt", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')
plt.ylim(0, 10000)

plt.subplot(3, 2, 5)
plt.hist(cutEB5.iHpT, bins = 100, range = (0, 0.5))
plt.xlabel("IsoHcal/pt", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')
plt.ylim(0, 10000)

plt.subplot(3, 2, 6)
plt.hist(EE.iHpT, bins = 100, range = (0, 0.5))
plt.xlabel("IsoHcal/pt", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')
plt.ylim(0, 10000)

plt.subplot(3, 2, 6)
plt.hist(cutEE5.iHpT, bins = 100, range = (0, 0.5))
plt.xlabel("IsoHcal/pt", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')
plt.ylim(0, 10000)
plt.savefig('Fig4.png')
plt.show()


# using all the cuts at once
EBcut = EB[(EB.sigmaEtaEta < 0.01) & (EB.HoverE < 0.04) & (EB.iTpT < 0.09) & (EB.iEpT < 0.07) & (EB.iHpT < 0.1)]
EEcut = EE[(EE.sigmaEtaEta < 0.031) & (EE.HoverE < 0.025) & (EE.iTpT < 0.05) & (EE.iEpT < 0.06) & (EE.iHpT < 0.025)]
# the events outside of cut ranges are considered background events
EBbg = EB[(EB.sigmaEtaEta > 0.01) | (EB.HoverE > 0.04) | (EB.iTpT > 0.09) | (EB.iEpT > 0.07) | (EB.iHpT > 0.1)]
EEbg = EE[(EE.sigmaEtaEta > 0.031) | (EE.HoverE > 0.025) | (EE.iTpT > 0.05) | (EE.iEpT > 0.06) | (EE.iHpT > 0.025)]


# before and after distributions for the EB region
figfinal = plt.figure(figsize = (10,15))

plt.subplot(4, 2, 1)
plt.xlim(0,200)
plt.hist(EB.pt, bins = 170)
plt.xlabel("pt [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)

plt.subplot(4, 2, 1)
plt.xlim(0,200)
plt.hist(EBcut.pt, bins = 125)
plt.xlabel("pt [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)

plt.subplot(4, 2, 2)
plt.hist(EB.pt, bins = 50)
plt.xlabel("pt [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')


plt.subplot(4, 2, 2)
plt.hist(EBcut.pt, bins = 50)
plt.xlabel("pt [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')


plt.subplot(4, 2, 3)
plt.hist(EB.eta, bins = 50)
plt.xlabel("eta", fontsize = 15)
plt.ylabel("events", fontsize = 15)

plt.subplot(4, 2, 3)
plt.hist(EBcut.eta, bins = 50)
plt.xlabel("eta", fontsize = 15)
plt.ylabel("events", fontsize = 15)

plt.subplot(4, 2, 4)
plt.hist(EB.eta, bins = 50)
plt.xlabel("eta", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')

plt.subplot(4, 2, 4)
plt.hist(EBcut.eta, bins = 50)
plt.xlabel("eta", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')

plt.subplot(4, 2, 4)
plt.hist(EB.MET, bins = 50, range = [-1.5, 1.5])
plt.xlabel("eta", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')



plt.subplot(4, 2, 5)
plt.xlim(0,100)
plt.hist(EB.MET, bins = 150)
plt.xlabel("MET [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)

plt.subplot(4, 2, 5)
plt.xlim(0,100)
plt.hist(EBcut.MET, bins = 150)
plt.xlabel("MET [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)

plt.subplot(4, 2, 6)
plt.hist(EB.MET, bins = 50)
plt.xlabel("MET [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')

plt.subplot(4, 2, 6)
plt.hist(EBcut.MET, bins = 50)
plt.xlabel("MET [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')


plt.subplot(4, 2, 7)
plt.xlim(0,200)
plt.hist(EB.mt, bins = 125)
plt.xlabel("mt [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)

plt.subplot(4, 2, 7)
plt.xlim(0,200)
plt.hist(EBcut.mt, bins = 85)
plt.xlabel("mt [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)

plt.subplot(4, 2, 8)
plt.hist(EB.mt, bins = 50)
plt.xlabel("mt [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')

plt.subplot(4, 2, 8)
plt.hist(EBcut.mt, bins = 50)
plt.xlabel("mt [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')

plt.savefig('Fig5.png')
plt.show()


# before and after distributions for the EE region

figfinal2 = plt.figure(figsize = (10,15))

plt.subplot(4, 2, 1)
plt.xlim(0,200)
plt.hist(EE.pt, bins = 150)
plt.xlabel("pt [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)

plt.subplot(4, 2, 1)
plt.xlim(0,200)
plt.hist(EEcut.pt, bins = 150)
plt.xlabel("pt [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)

plt.subplot(4, 2, 2)
plt.xlim(0,400)
plt.hist(EE.pt, bins = 80)
plt.xlabel("pt [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')


plt.subplot(4, 2, 2)
plt.xlim(0,400)
plt.hist(EEcut.pt, bins = 80)
plt.xlabel("pt [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')



plt.subplot(4, 2, 3)
plt.hist(EE.eta, bins = 50)
plt.xlabel("eta", fontsize = 15)
plt.ylabel("events", fontsize = 15)

plt.subplot(4, 2, 3)
plt.hist(EEcut.eta, bins = 50)
plt.xlabel("eta", fontsize = 15)
plt.ylabel("events", fontsize = 15)

plt.subplot(4, 2, 4)
plt.hist(EE.eta, bins = 50)
plt.xlabel("eta", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')

plt.subplot(4, 2, 4)
plt.hist(EEcut.eta, bins = 50)
plt.xlabel("eta", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')



plt.subplot(4, 2, 5)
plt.xlim(0,100)
plt.hist(EE.MET, bins = 100)
plt.xlabel("MET [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)

plt.subplot(4, 2, 5)
plt.xlim(0,100)
plt.hist(EEcut.MET, bins = 100)
plt.xlabel("MET [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)

plt.subplot(4, 2, 6)
plt.hist(EE.MET, bins = 50)
plt.xlabel("MET [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')

plt.subplot(4, 2, 6)
plt.hist(EEcut.MET, bins = 50)
plt.xlabel("MET [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')


plt.subplot(4, 2, 7)
plt.xlim(0,200)
plt.hist(EE.mt, bins = 125)
plt.xlabel("mt [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)

plt.subplot(4, 2, 7)
plt.xlim(0,200)
plt.hist(EEcut.mt, bins = 125)
plt.xlabel("mt [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)

plt.subplot(4, 2, 8)
plt.xlim(0,300)
plt.hist(EE.mt, bins = 85)
plt.xlabel("mt [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')

plt.subplot(4, 2, 8)
plt.xlim(0,300)
plt.hist(EEcut.mt, bins = 85)
plt.xlabel("mt [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')

plt.savefig('Fig6.png')
plt.show()


# acceptance: events after cuts / events before cuts
print(len(EBcut))
print(len(EEcut))
print("Acceptance EB {:}".format(len(EBcut)/len(EB)))
print("Acceptance EE: {:}".format(len(EEcut)/len(EE)))

# background acceptance
print("Akceptancja całej reszty do cięć na odwrót EB {:}".format(len(EBbg)/len(EB)))
print("Akceptancja całej reszty do cięć na odwrót EE {:}".format(len(EEbg)/len(EE)))


####################

figfinalEBbg = plt.figure(figsize = (10,10))

plt.subplot(2, 2, 1)
plt.xlim(0, 100)
plt.xlabel("MET [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)
np.random.seed(0)
data = np.random.normal(50, 20, 10000)
(counts, bins) = np.histogram(EBcut.MET, bins=120)
factor = 1
plt.hist(bins[:-1], bins, weights=factor*counts)

plt.subplot(2, 2, 1)
plt.xlim(0, 100)
plt.xlabel("MET [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)
np.random.seed(0)
data = np.random.normal(50, 20, 10000)
(counts, bins) = np.histogram(EBbg.MET, bins=100)
factor = 0.22
plt.hist(bins[:-1], bins, weights=factor*counts)

plt.subplot(2, 2, 2)
plt.xlim(0, 200)
plt.ylim(0, 10000)
plt.xlabel("MET [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')
np.random.seed(0)
data = np.random.normal(50, 20, 10000)
(counts, bins) = np.histogram(EBcut.MET, bins=100)
factor = 1
plt.hist(bins[:-1], bins, weights=factor*counts)

plt.subplot(2, 2, 2)
plt.ylim(0, 10000)
plt.xlim(0, 200)
plt.xlabel("MET [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')
np.random.seed(0)
data = np.random.normal(50, 20, 10000)
(counts, bins) = np.histogram(EBbg.MET, bins=100)
factor = 0.22
plt.hist(bins[:-1], bins, weights=factor*counts)





plt.subplot(2, 2, 3)
plt.xlim(0, 150)
plt.xlabel("mt [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)
np.random.seed(0)
data = np.random.normal(50, 20, 10000)
(counts, bins) = np.histogram(EBcut.mt, bins=100)
factor = 1
plt.hist(bins[:-1], bins, weights=factor*counts)

plt.subplot(2, 2, 3)
plt.xlim(0, 150)
plt.xlabel("mt [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)
np.random.seed(0)
data = np.random.normal(50, 20, 10000)
(counts, bins) = np.histogram(EBbg.mt, bins=100)
factor = 0.09
plt.hist(bins[:-1], bins, weights=factor*counts)

plt.subplot(2, 2, 4)
plt.xlim(0, 200)
plt.ylim(0, 1000)
plt.xlabel("mt [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')
np.random.seed(0)
data = np.random.normal(50, 20, 10000)
(counts, bins) = np.histogram(EBcut.mt, bins=100)
factor = 1
plt.hist(bins[:-1], bins, weights=factor*counts)

plt.subplot(2, 2, 4)
plt.xlim(0, 200)
plt.ylim(0, 1000)
plt.xlabel("mt [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')
np.random.seed(0)
data = np.random.normal(50, 20, 10000)
(counts, bins) = np.histogram(EBbg.mt, bins=100)
factor = 0.09
plt.hist(bins[:-1], bins, weights=factor*counts)

plt.savefig('Fig7.png')
plt.show()

######################

figfinalEEbg = plt.figure(figsize = (10,10))

plt.subplot(2, 2, 1)
plt.xlim(0, 80)
plt.xlabel("MET [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)
np.random.seed(0)
data = np.random.normal(50, 20, 10000)
(counts, bins) = np.histogram(EEcut.MET, bins=120)
factor = 1
plt.hist(bins[:-1], bins, weights=factor*counts)

plt.subplot(2, 2, 1)
plt.xlim(0, 80)
plt.xlabel("MET [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)
np.random.seed(0)
data = np.random.normal(50, 20, 10000)
(counts, bins) = np.histogram(EEbg.MET, bins=100)
factor = 0.27
plt.hist(bins[:-1], bins, weights=factor*counts)

plt.subplot(2, 2, 2)
plt.xlim(0, 100)
plt.ylim(0, 10000)
plt.xlabel("MET [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')
np.random.seed(0)
data = np.random.normal(50, 20, 10000)
(counts, bins) = np.histogram(EEcut.MET, bins=100)
factor = 1
plt.hist(bins[:-1], bins, weights=factor*counts)

plt.subplot(2, 2, 2)
plt.ylim(0, 10000)
plt.xlim(0, 100)
plt.xlabel("MET [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')
np.random.seed(0)
data = np.random.normal(50, 20, 10000)
(counts, bins) = np.histogram(EEbg.MET, bins=100)
factor = 0.27
plt.hist(bins[:-1], bins, weights=factor*counts)





plt.subplot(2, 2, 3)
plt.xlim(0, 150)
plt.xlabel("mt [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)
np.random.seed(0)
data = np.random.normal(50, 20, 10000)
(counts, bins) = np.histogram(EEcut.mt, bins=190)
factor = 1
plt.hist(bins[:-1], bins, weights=factor*counts)

plt.subplot(2, 2, 3)
plt.xlim(0, 150)
plt.xlabel("mt [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)
np.random.seed(0)
data = np.random.normal(50, 20, 10000)
(counts, bins) = np.histogram(EEbg.mt, bins=100)
factor = 0.25
plt.hist(bins[:-1], bins, weights=factor*counts)

plt.subplot(2, 2, 4)
plt.xlim(0, 200)
plt.ylim(0, 1000)
plt.xlabel("mt [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')
np.random.seed(0)
data = np.random.normal(50, 20, 10000)
(counts, bins) = np.histogram(EEcut.mt, bins=160)
factor = 1
plt.hist(bins[:-1], bins, weights=factor*counts)

plt.subplot(2, 2, 4)
plt.xlim(0, 200)
plt.ylim(0, 1000)
plt.xlabel("mt [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')
np.random.seed(0)
data = np.random.normal(50, 20, 10000)
(counts, bins) = np.histogram(EEbg.mt, bins=100)
factor = 0.25
plt.hist(bins[:-1], bins, weights=factor*counts)

plt.savefig('Fig8.png')
plt.show()

#############

EBscut = EBcut[(EBcut.mt > 50) & (EBcut.MET > 25)]
EEscut = EEcut[(EEcut.mt > 50) & (EEcut.MET > 25)]

EBsbg = EBbg[(EBbg.mt > 50) & (EBbg.MET > 25)]
EEsbg = EEbg[(EEbg.mt > 50) & (EEbg.MET > 25)]


print(len(EBscut))
print(len(EEscut))

print("Tlo po cieciach EB:", len(EBsbg))
print("EE: ", len(EEsbg))

print("eb finalne: ", len(EBscut) - math.ceil(len(EBsbg)*0.22))
print("ee finalne: ", len(EEscut) - math.ceil(len(EEsbg)*0.25))

eps_av = (len(EBscut)*0.798+len(EEscut)*0.67)/(len(EBscut)+len(EEscut))
print(eps_av)

Neb = len(EBscut) - math.ceil(len(EBsbg)*0.22)
Nee = len(EEscut) - math.ceil(len(EEsbg)*0.25)


print(Neb, Nee)

Np = Neb/0.4933/0.798 + Nee/0.4933/0.67

print(Np)

sigma = Np/36000
sigmaerr = 234/36000

print("cs: ", sigma, "nb +- ", sigmaerr)

print("-----")

print(len(EB[EB.Q==1]), len(EB[EB.Q==-1]))
print(len(EE[EE.Q==1]), len(EE[EE.Q==-1]))

print("-----")

print(len(cutEB1)/len(EB))
print(len(cutEB2)/len(EB))
print(len(cutEB3)/len(EB))
print(len(cutEB4)/len(EB))
print(len(cutEB5)/len(EB))
print(len(EBcut)/len(EB))
print("---")
print(len(cutEE1)/len(EE))
print(len(cutEE2)/len(EE))
print(len(cutEE3)/len(EE))
print(len(cutEE4)/len(EE))
print(len(cutEE5)/len(EE))
print(len(EEcut)/len(EE))

figfinalEBbg2 = plt.figure(figsize = (10,10))

plt.subplot(2, 2, 1)
plt.xlim(0, 100)
plt.xlabel("MET [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)
np.random.seed(0)
data = np.random.normal(50, 20, 10000)
(counts, bins) = np.histogram(EBbg.MET, bins=100)
factor = 1
plt.hist(bins[:-1], bins, weights=factor*counts)

plt.subplot(2, 2, 2)
# plt.ylim(0, 10000)
plt.xlim(0, 200)
plt.xlabel("MET [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')
np.random.seed(0)
data = np.random.normal(50, 20, 10000)
(counts, bins) = np.histogram(EBbg.MET, bins=100)
factor = 1
plt.hist(bins[:-1], bins, weights=factor*counts)

plt.subplot(2, 2, 3)
plt.xlim(0, 150)
plt.xlabel("mt [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)
np.random.seed(0)
data = np.random.normal(50, 20, 10000)
(counts, bins) = np.histogram(EBbg.mt, bins=100)
factor = 1
plt.hist(bins[:-1], bins, weights=factor*counts)

plt.subplot(2, 2, 4)
plt.xlim(0, 200)
# plt.ylim(0, 1000)
plt.xlabel("mt [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')
np.random.seed(0)
data = np.random.normal(50, 20, 10000)
(counts, bins) = np.histogram(EBbg.mt, bins=100)
factor = 1
plt.hist(bins[:-1], bins, weights=factor*counts)

plt.savefig('Fig9.png')
plt.show()

######################

figfinalEEbg2 = plt.figure(figsize = (10,10))

plt.subplot(2, 2, 1)
plt.xlim(0, 100)
plt.xlabel("MET [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)
np.random.seed(0)
data = np.random.normal(50, 20, 10000)
(counts, bins) = np.histogram(EEbg.MET, bins=100)
factor = 0.22
plt.hist(bins[:-1], bins, weights=factor*counts)

plt.subplot(2, 2, 2)
# plt.ylim(0, 10000)
plt.xlim(0, 200)
plt.xlabel("MET [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')
np.random.seed(0)
data = np.random.normal(50, 20, 10000)
(counts, bins) = np.histogram(EEbg.MET, bins=100)
factor = 0.22
plt.hist(bins[:-1], bins, weights=factor*counts)

plt.subplot(2, 2, 3)
plt.xlim(0, 150)
plt.xlabel("mt [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)
np.random.seed(0)
data = np.random.normal(50, 20, 10000)
(counts, bins) = np.histogram(EEbg.mt, bins=100)
factor = 0.09
plt.hist(bins[:-1], bins, weights=factor*counts)

plt.subplot(2, 2, 4)
plt.xlim(0, 200)
# plt.ylim(0, 1000)
plt.xlabel("mt [GeV]", fontsize = 15)
plt.ylabel("events", fontsize = 15)
plt.yscale('log')
np.random.seed(0)
data = np.random.normal(50, 20, 10000)
(counts, bins) = np.histogram(EEbg.mt, bins=100)
factor = 0.09
plt.hist(bins[:-1], bins, weights=factor*counts)

plt.savefig('Fig91.png')
plt.show()
