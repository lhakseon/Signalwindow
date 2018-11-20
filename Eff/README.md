# Manual for efficiency calculation

시그널 윈도우를 사용하여 이피션시와 백그라운드 레이트를 계산하는 방법에 대한 메뉴얼이 서술되어있다.
"/Signalwindow/README.md"의 메뉴얼을 따라하였다면 **Region of interest**, **Pixel-EM matching**, **Pixel-Pixel matching** 세가지 시그널 윈도우에 대한 피팅함수의 파라미터를 구했을 것이다. 
이피션시/레이트를 계산하기 위해선 어떻게 시그널 윈도우를 적용할수 있는지에 대해 서술되었다.


------------------------------------

## Configurations for efficiency measurement

이피션시의 측정은 200PU의 일렉트론 샘플을 사용하며 파일없이 없는 경우보다 계산에 많은 시간이 소요되기 때문에 job을 분할해서 서버에서 돌려야한다. 이것에 대한 환경설정방법이 아래의 스텝에 서술되어있다.


### Set the location of CMSSW framework at batch_script.py
서버에 job을 맞기기 위해 우선 아래의 코드에 CMSSW/root frame work의 위치를 지정해준다.
>Signalwindow/Eff/batch_scrpt.py

line 8
```
config+='source /cvmfs/cms.cern.ch/cmsset_default.sh \n'
config+='cd /u/user/moon/cmssw/CMSSW_8_0_26/src \n'
```

### Set the location of 200PU electron sample at inputlist_my.txt
아래 경로의 inputlist_my.txt 에 샘플이 분할된 모든 샘플의 위치를 아래의 예시와 같은 방식으로 기입해야 한다.
>Signalwindow/Eff/inputlist_my.txt

example)
```
dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/jhkim/SingleE_FlatPt-2to100/SingleElectron_200PU_PhaseIIFall17D/180829_174202/0000/SingleEle200PU__1.root
dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/jhkim/SingleE_FlatPt-2to100/SingleElectron_200PU_PhaseIIFall17D/180829_174202/0000/SingleEle200PU__10.root
dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/jhkim/SingleE_FlatPt-2to100/SingleElectron_200PU_PhaseIIFall17D/180829_174202/0000/SingleEle200PU__100.root
```
### Configuration setting of submit_control.py
sumit_cotrol.py의 아래 항목들을 사용자의 환경에 맞춰 수정한다.
>Signalwindow/Eff/submit_control.py
line 18 : inputlist_my.txt에 기입한 파일의 수
```
number_of_files = 267
```

line 22 : core의 갯수를 파일의 수와 동일하게 맞춰준다.
```
number_of_cores = 267
```

line 55 : 실행한 job이 저장될 workspace를 지정해준다.
```
workspace = "/u/user/syrian14/NewSW-Rate/Trackisolation/CMSSW9/PixTRK-Eff/ForKNU-eff/work/"
```
### setting pixel combination of each eta region. 
Singnalwindow 의 메뉴얼의 sw_dist.h와 동일한 방법으로 각각의 에타 영역에 사용하는 각각의 픽셀 클러스터를 정의한다.
아래의 코드에서 < > 로 처리된 부분을 이용자의 환경에 따라 수정해야 한다.
>Signalwindow/Eff/test

line 1234 : Set the pixel combination of each eta region at StorePixelHit.
line 1236 : this code is for pixel barrel cluster. 
```
if( region == 1 || region == 2 || region == 3 ){                ---<1>
    if( bRecHitLayer->at(a) == 2 ){ // Second layer             ---<2>
        if( Dphi < L2_Dphi_cut1 && Dphi > L2_Dphi_cut2){        ---<3>
            Dphi_Ele_pass = 1; el_or_po = 1;
        }
        if( Dphi > -L2_Dphi_cut1 && Dphi < -L2_Dphi_cut2){      ---<3>
            Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
        }
        if( Dphi_Ele_pass || Dphi_Pos_pass ){
            layers[2]++;                                        ---<3>
            second_layer_hits.push_back( TVector3(bRecHitGx->at(a), bRecHitGy->at(a), bRecHitGz->at(a)));
            second_layer_hits_Ele_or_Pos.push_back(el_or_po);   ---<3>
        }
    }
```
    <1> : Write the eta region which the pixel cluster is used.
    <2> : Define what pixel cluster is used in this code.
    <3> : [ L2_Dphi, layer[2], second ] means that among the 4 pixel cluster in each region, this pixel is secondary closest pixel from beam spot.


### test.C 시그널 윈도우 디테일 설정.

1. 시그널 윈도우 폭의 수정
아래의 파라미터는 피팅 함수가 가질 윈도우의 폭을 설정되어 있다. default 상태에선 두번째 배열까지의 값을 사용한다.
    line 24 ~ 29 : 
```
const double EM_PiX_dphi_width_[9] = {0.02, 0.03, 0.022, 0.024, 0.026, 0.028, 0.030, 0.032, 0.034}; //
const double EM_PiX_deta_width_[9] = {0.01, 0.015, 0.011, 0.012, 0.013, 0.014, 0.015, 0.16, 0.17}; //

const double PiX_PiX_dphi_width_[9] = {0.0017, 0.003, 0.0018, 0.0019, 0.0020, 0.0021, 0.0022, 0.0023, 0.0024};
const double PiX_PiX_deta_width_[9] = {0.0017, 0.003, 0.0018, 0.0019, 0.0020, 0.0021, 0.0022, 0.0023, 0.0024};
```

2. eta region의 범위 설정
아래의 코드에서 각 eta region의 범위를 정의해야 한다.
line 245

```
if( fabs(EgEta) <= 0.8 ) eta_region =1;
if( fabs(EgEta) <= 1.4 && fabs(EgEta) > 0.8 ) eta_region =2;
if( fabs(EgEta) <= 1.7 && fabs(EgEta) > 1.4 ) eta_region =3;
if( fabs(EgEta) <= 2.1 && fabs(EgEta) > 1.7 ) eta_region =4;
if( fabs(EgEta) <= 2.7 && fabs(EgEta) > 2.1 ) eta_region =5;
if( fabs(EgEta) <= 3.0 && fabs(EgEta) > 2.7 ) eta_region =6;

if( fabs(EgEta) > 3. ) continue;
```

3. 시그널 윈도우의 넓이 지정
아래의 코드를 수정해 각 에타 영역별 시그널 윈도우의 폭을 적용할 수 있다.
line 293
```
       if(eta_region == 1) SetSingalBoundary(eta_region, EM_PiX_dphi_width_[nth_eg_pix_deta], EM_PiX_deta_width_[nth_eg_pix_deta], PiX_PiX_dphi_width_[nth_eg_pix_deta], PiX_PiX_deta_width_[nth_eg_pix_deta]);
       else if(eta_region == 2) SetSingalBoundary(eta_region, EM_PiX_dphi_width_[nth_eg_pix_deta+1], EM_PiX_deta_width_[nth_eg_pix_deta+1], PiX_PiX_dphi_width_[nth_eg_pix_deta+1], PiX_PiX_deta_width_[nth_eg_pix_deta+1]);
       else SetSingalBoundary(eta_region, EM_PiX_dphi_width_[nth_eg_pix_deta], EM_PiX_deta_width_[nth_eg_pix_deta], PiX_PiX_dphi_width_[nth_eg_pix_deta], PiX_PiX_deta_width_[nth_eg_pix_deta]);
```

### Run the Efficiency measure code.
파이선 명령어로 아래의 코드를 실행하여 job을 실행할 수 있다. 
코드가 실행완료 되면 workspace/.../output/Tree/ 경로 상에 Tree_SE_PU200_100.root 와 같은 파일이 생성된다.
모든 파일을 hadd 를 이용해 묶고 efficiency drawing code의 input으로 사용하면 efficeincy를 구할수 있다.

```
python submit_control.py
```
위의 코드가 전부 돌아가는 데에 최대 한시간 정도가 소요된다.
코드가 전부 완료되면 아래의 경로에서 hadd 명령어를 이용해 결과 파일을 묶는다.
```
hadd Result_all.root workspace/.../output/Tree/*.root
```

작성중 







