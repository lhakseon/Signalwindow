# Signalwindow


This manual is based on CMSSW_10_1_X detector geometry.


Table of contents

1. Defining eta region for signal windows based on 3 out of 4 efficiency.
2. Defining Signal windows using single electron samples.

이 매뉴얼은 3 out of 4 이피션시를 이용해 각 영역별로 시그널 윈도우를 만들때 고려할 픽셀 조합을 선택하는 방법과
파일업이 없는 싱글 일렉트론 샘플을 이용해 시그널 윈도우를 정의하는 방법에 대해 서술되어 있다.
각각의 스텝엔 시그널 윈도우를 얻기 위해 무엇을 해야 하는지에 대해 단계순으로 간단히 설명 되어 있다.
스텝 하단엔 수정/실행을 위한 코드의 위치가 적혀있고 그 아래에 사용자의 환경에 맞게 수정해야할 부분을  < >로 표시하였다.  



------------------------------------


## [1] Defining eta region for signal windows based on 3 out of 4 efficiency.

시그널 윈도우를 만들땐 네개의 픽셀 클러스터를 활용하게 되는데 gen particle의 에타에 따라 사용하는 픽셀 클러스터의 조합이 달라진다.
각 에타 영역에서 사용하는 픽셀 조합은 3 out of 4 이피션시를 기준으로 판별하는데 가장 높은 이피션시를 가지는 픽셀 조합을 시그널 윈도우를 만들때 활용한다.
이 파트는 시그널 윈도우를 정의할때 사용할 3 out of 4 픽셀 이피션시의 그래프를 만드는 매뉴얼로 4개의 픽셀 배럴과 2개의 픽셀 디스크에 최적화 되어있다.







### step1 : set the location of sample for efficiency check.
아래의 경로에 efficiency check에 사용할 No-pileup electron sample의 위치를 지정한다.
>signal_windows_/signal_windows/basicm/basic_study.h

line 388 ~ 394 

```
   if (tree == 0) {
       TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("./SingleElectron_NoPU.root");      --- <1>
           if (!f || !f->IsOpen()) {                                                                            
               f = new TFile("./SingleElectron_NoPU.root");                                       --- <1>               
           }
       TDirectory * dir = (TDirectory*)f->Get("./SingleElectron_NoPU.root:/l1PiXTRKTree");        --- <2>
           dir->GetObject("L1PiXTRKTree",tree);                                                   --- <2>

```

<1> : write the location of NOPU single electron sample.     
<2> : write the name of smaple's tree.


### step 2 : Run the basic_study.C 
basic_study.C를 터미널에서 실행시켜 Eff_nopu.root를 얻는다.
>signal_windows/signal_windows/basicm/basic_study.C

This code considers 1 to 4 pixel layers and 1 to 5 disks. To use additional disk, the relevant code should be added.
Run the following code as 

```
root basic_study.C
basic_study a
a.Loop()
```

**Eff_nopu.root** is generated.
use this root file as a input of efficiency_check.C

### step 3 : Set the location of Eff_nopu.root at efficiency_check.C and run the code.
Eff_nopu.root가 저장된 위치를 efficinecy_check.C에 지정한다.
>/signal_windows/signal_windows/basicm/basic_efficiency_check.C

- This code doesn't need a head files.

line 13 
```
TFile* histFile = new TFile("./Eff_nopu.root");   --- <1>
```
<1> : Set the location of input file(Eff_nopu.root) in here. 

```
root efficiency_check.C
```
After running the code **effcheck.png** is generated.

You can find the Devided eta region checking effcheck.png.
Select the most highest efficinecy of pixel combination in each eta region.

----------------------------

## [2] Defining Signal windows using single electron samples.



### step 1 : Input the location of No-PU particle sample
Input으로 사용할 No-PU particle 샘플의 경로를 sw_dist.h에 지정해준다.
>/signal_windows/signal_windows/sw_dist.h
line 454

```
TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("./SingleElectron_NoPU.root");  
      if (!f || !f->IsOpen()) {
          f = new TFile("./SingleElectron_NoPU.root");                                   
      }
      TDirectory * dir = (TDirectory*)f->Get("./SingleElectron_NoPU.root:/l1PiXTRKTree");
      dir->GetObject("L1PiXTRKTree",tree);

   }
```


### step 2 : At the StorePixelHit, set the pixel combinations.
sw_dist.h에서 각각의 에타 영역에 대한 pixel combination을 정의한다.
>/signal_windows/signal_windows/sw_dist.h
```
void sw_dist::StorePixelHit(int region){
...
for(int a=0; a<bRecHitN; a++){                                                  
TVector3 current_hit;
current_hit.SetXYZ( bRecHitGx->at(a), bRecHitGy->at(a), bRecHitGz->at(a) );
```    
- these lines define only barrel pixel clusters.
to define disks of each eta regions, use the following code.

```
for(int a=0; a<fRecHitN; a++){                                                  
TVector3 current_hit;
current_hit.SetXYZ( fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a) );
```

The following is an example of defining pixel clusters for each eta region.

```
if( region == 1 ){                              ---<1>
    if (bRecHitLayer->at(a) == 1 ){             ---<2>
        layers[1]++;                            ---<3>
        first_layer_hits.push_back(current_hit) ---<4>
    }
```

<1> : region number. Each eta region should have four pixel clusters(layers+disks). 

<2> : if the pixel is second barrel, bRecHitLayer->at(a) == 2
        if pixel cluster is 1st Disk, fRecHitDisk.->at(a) == 1

<3> :
    the most closest pixel cluster : layers[1]
    secondary closest pixel cluseter : layers[2]
    ...

<4> : Match pair with <3> like this way -> first, second, third, fourth.




### step3 : defining eta range of each eta regions
sw_dist.C에서 eta region의 범위를 정의한다.
>/signal_windows/signal_windows/sw_dist.C

line 122
```
int eta_region = 0;
if( fabs(EgEta) < 0.8 ) eta_region =1;
if( fabs(EgEta) < 1.4 && fabs(EgEta) > 0.8 ) eta_region =2;
if( fabs(EgEta) < 1.7 && fabs(EgEta) > 1.4 ) eta_region =3;
if( fabs(EgEta) < 2.1 && fabs(EgEta) > 1.7 ) eta_region =4;
if( fabs(EgEta) < 2.7 && fabs(EgEta) > 2.1 ) eta_region =5;
if( fabs(EgEta) < 3.0 && fabs(EgEta) > 2.7 ) eta_region =6;
if( fabs(EgEta) > 3. ) continue;
```
Set the range of each eta regions.



### step 4 : run sw_dist.C
터미널에서 sw_dist.C를 실행하여 sw.root 파일을 생성한다.
At terminal
```
root sw_dist.C
sw_dist a
a.Loop()
```
from this process, sw.root is created.
sw.root is input file of following codes.




### step 5 : set the location of sw.root at Make2Dplots.h
step4에서 생성된 sw.root를 Make2Dplot.h의 인풋으로 지정한다.
- this step focus on Delta Phi signal windows. 

>signal_windows/signal_windows/fit_median/Make2Dplots.h

line 127
```
if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../sw.root");   --- <1>
    if (!f || !f->IsOpen()) {
        f = new TFile("../sw.root");
    }
    f->GetObject("t",tree);

}
Init(tree);
```
<1> : input the location of sw.root


### step 6 : set the eta range of each eta region.
Make2Dplots.C에서 각각의 에타 영역의 범위를 정의한다.
>/signal_windows/signal_windows/fit_median/Make2Dplots.C

line 198
```
if( eta_ == 1 && fabs(ntEgEta->at(0)) > 0.8 ) continue; 
if( eta_ == 2 && (fabs(ntEgEta->at(0)) < 0.8 || fabs(ntEgEta->at(0)) > 1.4)) continue; 
if( eta_ == 3 && (fabs(ntEgEta->at(0)) < 1.4 || fabs(ntEgEta->at(0)) > 1.7)) continue; 
if( eta_ == 4 && (fabs(ntEgEta->at(0)) < 1.7 || fabs(ntEgEta->at(0)) > 2.1)) continue; 
if( eta_ == 5 && (fabs(ntEgEta->at(0)) < 2.1 || fabs(ntEgEta->at(0)) > 2.7)) continue; 
if( eta_ == 6 && (fabs(ntEgEta->at(0)) < 2.7 || fabs(ntEgEta->at(0)) > 3.0)) continue; 
```
**eta_** is the eta region.
Set the range of each eta regions as the same way of step 3.

### step 7 : Specify eta region to draw and get the fitting function of signal window.
x_phi.C에서 시그널 윈도우를 출력할 eta region을 정의한다. run.sh를 실행시키면 시그널 윈도우의 피팅 함수를 얻을 수 있다.
>/signal_windows/signal_windows/fit_median/x_phi.C
```
.L Make2Dplots.C+
Make2Dplots a
a.Loop(1)
    //a.Loop(2)         --- <1>
    //a.Loop(3)
    //a.Loop(4)
    //a.Loop(5)
a.Loop(6)
```
<1> : a.Loop(x) -> x means the eta region to be output of signal windows.
this example print out the signal windows of regino1 and 6. 
To get other result of eta region, just delete //.

After defined the eta region try the following code.

### step 8 : run the following code and get copy the fitting function.
아래의 코드를 실행하면 delta phi에 대한 시그널 윈도우의 피팅함수와 파라미터를 획득할 수 있다. 

```
./run.sh
```
than the parameters of fitting function for delta phi signal windows are generated.

```
ROI_...txt
EGmatching...txt
Pixelmatching...txt
```


### step 9-1 : Copy the parameters and paste it at the bellow files.
step 8에서 생성된 피팅함수를 복사해서 다음 위치의 파일에 붙여넣는다..
```
ROI_...txt -> RegionOfInterst.h
EGmatching...txt ->  withEM.h 
Pixelmatching...txt   -> withoutEM.h
```

각각의 텍스트 파일에서 아래의 텍스트를 복사해야 한다.
```
if( region == 2 && i == 0 ){
    p[0] = 0.000924974;
    p[1] = -1.80158;
    p[2] = -0.77993;
    p[3] = 0.126864;
}
```


#step 9-2.
복사한 파라미터를 붙여넣을때 고려해야 할 사항
RegionOfInterst.h:
```
double ROI_func(int region, double eget){
    if(region == 1){
        // fit for median
        p[0] = 0.00130283;
        p[1] = -0.317612;
        p[2] = -0.942851;
        p[3] = -0.107378;
        p[4] =  1.27837;
    }
```
Region of Interest는 step 9에서 i는 무시하고 오직 region만을 고려한다.

withEM.h
   

```
double SW_func2_dphi_v2(int region, int nth, double eget){
    if( region == 1 && nth == 0 ){   --- <1>
        p[0] = -3.76446e-05;
        p[1] = 1.15386;
        p[2] = -1.08511;
        p[3] = -0.818285;
    }
```

    EM-Pixel matching의 시그널 윈도우를 정의하는 코드이다.region과 함께 <1>의 nth에 픽셀 조합에 따라 다른 숫자를 기입한다.
    각각의 숫자에 따른 EM-Pixel 조합은 아래와 같다.
    <1> : nth == x : x -> 0 : EM_Pixel12, 1 : EM_Pixel13, 2 : EM_Pixel14, 3 : EM_Pixel23, 4: EM_Pixel24, 5: EM_Pixel34
    
withoutEM.h



```
double SW_func1_dphi_v2(int region, int nth, double eget){
f( region == 1 && nth == 0 ){       --- <1>
    p[0] = 0.000215938;
    p[1] = 0.0608554;
    p[2] = -0.267592;
    p[3] = 0.31785;
}
```
마찬가지로 nth에 Pixel-Pixel matching의 조합에 따라 <1>의 nth에 다른 숫자를 기입한다.
<1> : nth == x : x -> 0 : P012, 1 : P013, 2 : P014, 3 : P023, 4: P024, 5: P034, 6: P123, 7: P124, 8: P134, 9: P234
0 means beam spot.

