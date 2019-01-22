# WSFProducer
W tagging SF producer

Install:
```
cmsrel CMSSW_9_4_2
cd CMSSW_9_4_2/src
cmsenv
git clone https://github.com/cmantill/WSFProducer.git
cd WSFProducer/
python setupPDF.py
```

Run:
```
python fitW.py --tag XXX 

-b: BATCH mode
--xMin: min mass
--xMax: max mass
--fits: run individual fits to each MC sample
--comb: run combined data and MC fit
--ind: run pass	and fail fits not simult
```

