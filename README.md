# Signalwindow


This manual is based on CMSSW_10_1_X detector geometry.

[1] Defining eta region for signal windows based on 3 out of 4 efficiency.

1)
signal_windows_/signal_windows/basicm/basic_study.h

line 388 ~ 394 

-----------------------------------------
   if (tree == 0) {
             TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("./SingleElectron_NoPU.root");                      --- <1>
                   if (!f || !f->IsOpen()) {                                                                            
                                f = new TFile("./SingleElectron_NoPU.root");                                             
                                      }
                         TDirectory * dir = (TDirectory*)f->Get("./SingleElectron_NoPU.root:/l1PiXTRKTree");            --- <2>
                               dir->GetObject("L1PiXTRKTree",tree);

------------------------------------------

<1> : location of NOPU single electron sample.     
<2> : Name of smaple's tree.


2)
signal_windows/signal_windows/basicm/basic_study.C

This code considers 1~4 pixel layers and 1~5 disks. To use additional disk, the relevant code should be added.

------------------
root basic_study.C
basic_study a
a.Loop()
------------------

Eff_nopu.root is generated.
use this root file as a input of efficiency_check.C

3)
/signal_windows/signal_windows/basicm/basic_efficiency_check.C
This code doesn't need a head files.

line 13 
-------------------------------------------
TFile* histFile = new TFile("./Eff_nopu.root");   --- <1>
-------------------------------------------
<1> : Input of 2)

------------------
root efficiency_check.C
------------------
after running the code, effcheck.png is generated.

Devide the eta region based on the highest pixel combination efficiency.


[2]
Defining Signal windows using single electron samples.


/signal_windows/signal_windows/sw_dist.h

1)
line 454
Input the location of single particle samples.

-----------------------------------------
TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("./SingleElectron_NoPU.root");
      if (!f || !f->IsOpen()) {
          f = new TFile("./SingleElectron_NoPU.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("./SingleElectron_NoPU.root:/l1PiXTRKTree");
      dir->GetObject("L1PiXTRKTree",tree);

   }
------------------------------------------


2)
at the StorePixelHit, set the pixel combinations for each eta based on [1].

-------------------------------------------------------------
void sw_dist::StorePixelHit(int region){
...
for(int a=0; a<bRecHitN; a++){                                                  
TVector3 current_hit;
current_hit.SetXYZ( bRecHitGx->at(a), bRecHitGy->at(a), bRecHitGz->at(a) );
---------------------------------------------------------------    
*these lines define only barrel pixel clusters.
to define disks of each eta regions, use the following code.

--------------------------------------------------------------
for(int a=0; a<fRecHitN; a++){                                                  
TVector3 current_hit;
current_hit.SetXYZ( fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a) );
-------------------------------------------------------------

The following is an example of defining pixel clusters for each eta region.

-----------------------------------------------------------
if( region == 1 ){                              ---<1>
    if (bRecHitLayer->at(a) == 1 ){             ---<2>
        layers[1]++;                            ---<3>
        first_layer_hits.push_back(current_hit) ---<4>
    }
------------------------------------------------------

<1> : region number. Each eta region should have four pixel clusters(layers+disks). 
<2> : if the pixel is second barrel, bRecHitLayer->at(a) == 2
        if pixel cluster is 1st Disk, fRecHitDisk.->at(a) == 1;
<3> :
    the most closest pixel cluster : layers[1]
    secondary closest pixel cluseter : layers[2]
    ...

<4> : Match pair with <3>. first, second, third, fourth.

-------------------------------------------------------------

/signal_windows/signal_windows/sw_dist.C

3)
defining eta range of each eta regions


line 122
--------------------------------------------------------------
int eta_region = 0;
if( fabs(EgEta) < 0.8 ) eta_region =1;
if( fabs(EgEta) < 1.4 && fabs(EgEta) > 0.8 ) eta_region =2;
if( fabs(EgEta) < 1.7 && fabs(EgEta) > 1.4 ) eta_region =3;
if( fabs(EgEta) < 2.1 && fabs(EgEta) > 1.7 ) eta_region =4;
if( fabs(EgEta) < 2.7 && fabs(EgEta) > 2.1 ) eta_region =5;
if( fabs(EgEta) < 3.0 && fabs(EgEta) > 2.7 ) eta_region =6;
if( fabs(EgEta) > 3. ) continue;
-------------------------------------------------------------
Set the range of each eta regions.



4)
at terminal
--------------------
root sw_dist.C
sw_dist a
a.Loop()
--------------------
from this process, sw.root is created.
sw.root is input file of following codes.




5) Signal windows defining.
-------------------------------------------------------------------
Delta Phi signal windows
------------------------------------------------------------------
signal_windows/signal_windows/fit_median/Make2Dplots.h

line 127
------------------------------------------------------------------
if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../sw.root");   --- <1>
    if (!f || !f->IsOpen()) {
        f = new TFile("../sw.root");
    }
    f->GetObject("t",tree);

}
Init(tree);
------------------------------------------------------------------
<1> : input the location of sw.root which is made from 4) 


6) Delta Phi signal windows defining
/signal_windows/signal_windows/fit_median/Make2Dplots.C

line 198
---------------------------------------------------------------------
if( eta_ == 1 && fabs(ntEgEta->at(0)) > 0.8 ) continue; 
if( eta_ == 2 && (fabs(ntEgEta->at(0)) < 0.8 || fabs(ntEgEta->at(0)) > 1.4)) continue; 
if( eta_ == 3 && (fabs(ntEgEta->at(0)) < 1.4 || fabs(ntEgEta->at(0)) > 1.7)) continue; 
if( eta_ == 4 && (fabs(ntEgEta->at(0)) < 1.7 || fabs(ntEgEta->at(0)) > 2.1)) continue; 
if( eta_ == 5 && (fabs(ntEgEta->at(0)) < 2.1 || fabs(ntEgEta->at(0)) > 2.7)) continue; 
if( eta_ == 6 && (fabs(ntEgEta->at(0)) < 2.7 || fabs(ntEgEta->at(0)) > 3.0)) continue; 
------------------------------------------------------------------------------------------------------------
the eta_ means eta region.
Set the range of each eta regions ust as 3)

7)
/signal_windows/signal_windows/fit_median/x_phi.C
----------------------
.L Make2Dplots.C+
Make2Dplots a
a.Loop(1)
    //a.Loop(2)         --- <1>
    //a.Loop(3)
    //a.Loop(4)
    //a.Loop(5)
a.Loop(6)
-----------------------------------------------------------------
<1> : a.Loop(x) -> x means the eta region to be output of signal windows.
this example print out the signal windows of regino1 and 6. 
To get other result of eta region, just delete //.



8)
at terminal
------------
./run.sh
------------
than the parameters of fitting function for delta phi signal windows are generated.
Copy the parameters and paste it at the bellow files.

ROI_...txt -> RegionOfInterst.h
EGmatching...txt ->  withEM.h 
Pixelmatching...txt   -> withoutEM.h

if you open one of this txt files
--------------------------------
if( region == 2 && i == 0 ){
    p[0] = 0.000924974;
    p[1] = -1.80158;
    p[2] = -0.77993;
    p[3] = 0.126864;
}
-------------------------------
Copy this codes.


9)
paste paremeters from 8) at following functions.

RegionOfInterst.h:
------------------------------------------
double ROI_func(int region, double eget){
    if(region == 1){
        // fit for median
        p[0] = 0.00130283;
        p[1] = -0.317612;
        p[2] = -0.942851;
        p[3] = -0.107378;
        p[4] =  1.27837;
    }
----------------------------------------

    withEM.h
---------------------------------------
double SW_func2_dphi_v2(int region, int nth, double eget){
    if( region == 1 && nth == 0 ){   --- <1>
        p[0] = -3.76446e-05;
        p[1] = 1.15386;
        p[2] = -1.08511;
        p[3] = -0.818285;
    }
----------------------------------------------------
<1> nth == x : x -> 0 : EM_P12, 1 : EM_P13, 2 : EM_P14, 3 : EM_P23, 4: EM_P24, 5: EM_P34
    
    withoutEM.h
------------------------------------------------------
double SW_func1_dphi_v2(int region, int nth, double eget){
f( region == 1 && nth == 0 ){       --- <1>
    p[0] = 0.000215938;
    p[1] = 0.0608554;
    p[2] = -0.267592;
    p[3] = 0.31785;
}
------------------------------------------------------
<1> nth == x : x -> 0 : P012, 1 : P013, 2 : P014, 3 : P023, 4: P024, 5: P034, 6: P123, 7: P124, 8: P134, 9: P234


