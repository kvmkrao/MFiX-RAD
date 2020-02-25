!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: RADIATION SPECTRAL                                     !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Authos: V Kotteda, M Stoellinger                   Date: 1-july-19  !
!   Comments                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
module rad_spectral
use rad_param
use rad_config
use rad_fields
use rad_gas_species, only : gasInfoType, rad_gas_species_init, getGasInfo
use rad_solid_species, only : solidParInfoType, rad_solid_species_init, getSolidInfo, smax
use rad_solid_species, only : getDesSolidInfo
implicit none
private

! ---- Interfaces ------
public :: rad_spectral_init
public :: rad_spectral_calc
public :: rad_spectral_srad
public :: rad_spectral_final

interface rad_spectral_init
    module procedure init1
end interface

interface rad_spectral_final
    module procedure finalize1
end interface

interface rad_spectral_calc
    module procedure calculateSpectralFields
end interface

interface rad_spectral_srad
    module procedure calcRadSources
end interface
! ---- Data members ----
character(len=strLen) :: modelName_
integer :: modelId_

! ---- Parameters ------
! model ids
integer, parameter :: unknown=0, gray=1
integer, parameter :: graycon=2,graywsgg=3, nongraywsgg=4 

! smith 1982    http://doi.org/10.1115/1.3245174
! P_H20/p_CO2 = 1 
real(dp), dimension(3) :: ki1
real(dp), dimension(3,4) :: bei1
! P_H20/p_CO2 = 2
real(dp), dimension(3) :: ki2
real(dp), dimension(3,4) :: bei2

! BORDBER2014    https://doi.org/10.1016/j.combustflame.2014.03.013  
real(dp), dimension(4,5,5) :: cijk
real(dp), dimension(4,5) :: dik
! JOHANSSON2011  https://doi.org/10.1016/j.combustflame.2011.02.001
real(dp), dimension(4,3) :: c1ik,c2ik,c3ik
real(dp), dimension(4)   :: k1j, k2j
! Dorigon2013    https://doi.org/10.1016/j.ijheatmasstransfer.2013.05.010 
real(dp), dimension(4) :: kpi1, kpi2
real(dp), dimension(4) :: b1i0, b1i1, b1i2,b1i3,b1i4 
real(dp), dimension(4) :: b2i0, b2i1, b2i2,b2i3,b2i4 
! TAYLOR74       https://doi.org/10.1016/0017-9310(74)90067-2 
real(dp), dimension(4)  :: kf2i, bf21i, bf22i 
real(dp), dimension(4)  :: kf1i, bf11i, bf12i

!   Shan 2017
real(dp), dimension(4,3) :: d55ik
real(dp), dimension(4,3) :: d51ik
real(dp), dimension(4,3) :: d52ik

real(dp), dimension(4,3) :: d105ik
real(dp), dimension(4,3) :: d101ik
real(dp), dimension(4,3) :: d102ik

real(dp), dimension(4,3) :: d155ik
real(dp), dimension(4,3) :: d151ik
real(dp), dimension(4,3) :: d152ik

real(dp), dimension(4,8,3) :: c55ijk
real(dp), dimension(4,8,3) :: c51ijk
real(dp), dimension(4,8,3) :: c52ijk

real(dp), dimension(4,8,3) :: c105ijk
real(dp), dimension(4,8,3) :: c101ijk
real(dp), dimension(4,8,3) :: c102ijk

real(dp), dimension(4,8,3) :: c155ijk
real(dp), dimension(4,8,3) :: c151ijk
real(dp), dimension(4,8,3) :: c152ijk

! smith1982 
data ki1  / 0.4303, 7.055, 178.1 /
data bei1 /5.150,  0.7749, 1.907, &
          -2.303,  3.399, -1.824, &
           0.9779,-2.297,  0.5608, &
          -1.494,  3.770, -0.5122 /
data ki2 / 0.4201, 6.516, 131.9 /
data bei2 / 6.508, -0.2504, 2.718, &
            -5.551, 6.112, -3.118, &
             3.029,-3.882, 1.221, &
            -5.353, 6.528,-1.612 /

!  TAYLOR74  
data kf1i  / 0.0, 0.91, 9.4, 130.0/
data kf2i  / 0.0, 0.69, 7.4, 80.0 /
data bf21i / 0.364, 0.266, 0.252, 0.118 /
data bf22i / 4.74, 7.19, -7.41, -4.52 /

data bf11i/ 0.4092, 0.284, 0.211,  0.0958/ 
data bf12i / 7.53,    2.58, -6.54, - 3.57/

! JOHANSSON2011 
data c1ik / 0.358 , 0.0731,-0.0466, &
            0.392, -0.212,  0.0191, &
            0.142, -0.0831, 0.0148, &
            0.0798,-0.0370, 0.0023/
data c2ik / -0.165,-0.0554, 0.0930, &
            -0.291, 0.644, -0.209, &
            0.348,  -0.294, 0.0662, &
            0.0866, -0.106,  0.0305 /
data c3ik / 0.0598,  0.0028,-0.0256, &
            0.0784, -0.197,  0.0662, &
            -0.122,  0.118, -0.0295, &
            -0.0127, 0.0169,-0.0051 /

data k1j / 0.055, 0.88, 10.0, 135 /
data k2j / 0.012,-0.021,-1.6, -35 /

! Dorigon2013
data kpi2 / 0.1921, 1.719, 11.37, 111.016/
data kpi1 / 0.01873, 1.723, 12.48,  144.9/

data b1i0 / 0.07197 , 0.11070, 0.20910, 0.07092/
data b1i1 / 87.24, 33.97, -6.423, 6.586/
data b1i2 / -96.90, -24.67, -3.200, -12.78/ 
data b1i3 / 46.51, 4.647, 1.718, 5.577/ 
data b1i4 / -79.17, -1.039, -2.105, -7.709/ 

data b2i0 / 0.05617, 0.14260, 0.13620, 0.12220/
data b2i1 / 78.44  , 17.95  , 25.74  , -2.327/
data b2i2 / -85.63 , -1.077 , -37.11 , -7.492/
data b2i3 / 42.46  , -6.971 , 15.75  , 4.275 /
data b2i4 / -74.4  , 17.74  , -22.67 , -6.608/

! Borddar et al. (2014), A line by line based weighted sum of gray gases model for inhomogeneous
! co2-h2o mixture in oxy-fired combustion 
data dik  / 0.0340429, 0.3509457, 4.5707400, 109.81690, &
            0.0652305, 0.7465138, 2.1680670,-50.923590, &
           -0.0463685,-0.5293090,-1.4989010, 23.432360, &
           0.0138684, 0.1594423,  0.4917165, -5.1638920, &
          -0.0014450,-0.0166326, -0.0542999, 0.4393889/

data cijk /0.7412956, 0.1552073, 0.2550242,-0.0345199, &
          -0.9412652, 0.6755648,-0.6065428, 0.4112046, &
           0.8531866,-1.1253940, 0.8123855,-0.5055995, &
          -0.3342806, 0.6040543,-0.4532290, 0.2317509, &
           0.0431436,-0.1105453, 0.0869309,-0.0375491, &
          -0.5244441,-0.4862117, 0.3805403, 0.2656726, &
           0.2799577, 1.4092710, 0.3494024,-0.5728350, &
           0.0823075,-0.5913199,-1.1020091, 0.4579559, &
           0.1474987,-0.0553385, 0.6784475,-0.1656759, &
          -0.0688622, 0.0464663,-0.1306996, 0.0229520, &
           0.5822860, 0.3668088,-0.4249709,-0.1225365, &
          -0.7672319,-1.3834490, 0.1853509, 0.2924490, &
           0.5289430, 0.9085441, 0.4046178,-0.2616436, &
          -0.4160689,-0.1733014,-0.3432603, 0.1052608, &
           0.1109773,-0.0016129, 0.0741446,-0.0160047, &
          -0.2096994,-0.1055508, 0.1429446, 0.0300151, &
           0.3204027, 0.4575210,-0.1013694,-0.0798076, &
          -0.2468463,-0.3334201,-0.0811822, 0.0764841, &
           0.1697627, 0.0791608, 0.0883088,-0.0321935, &
          -0.0420861,-0.0035398,-0.0202929, 0.0050463, &
           0.0242031, 0.0105857,-0.0157408,-0.0028205, &
          -0.0391017,-0.0501976, 0.0130244, 0.0079966, &
           0.0310940, 0.0384236, 0.0062981,-0.0079084, &
          -0.0204066,-0.0098934,-0.0084152, 0.0033870, &
           0.0049188, 0.0006121, 0.0020110,-0.0005364 /

data d55ik/ 0.006029 ,  -0.199053 ,  -1.418661 ,   4.459971 , &
           -0.011827 ,   0.328104 ,   1.337215 ,  -1.175645 , &
            0.084915 ,   0.382368 ,   1.617068 ,  16.216748 /

data c55ijk/ 1.872342 ,  -2.064489 ,   2.349744 ,  -0.843703 , &
           -15.340838 ,  14.403190 , -18.160530 ,   7.530383 , &
            51.170876 , -38.566516 ,  57.766407 , -28.512409 , &
           -89.520361 ,  48.979413 , -99.317202 ,  59.715968 , &
            88.568872 , -28.818367 , 102.518367 , -75.000915 , &
           -50.317662 ,   5.067608 , -64.856238 ,  55.881863 , &
            16.119345 ,   1.523995 ,  21.494650 , -21.314173 , &
            -2.067575 ,  -0.608189 ,  -1.935853 ,   2.038943 , &
            -1.116641 ,   0.328847 ,  -0.810192 ,   0.507078 , &
             8.782650 ,  -0.115866 ,   5.482606 ,  -4.624084 , &
           -27.872438 ,  -9.488601 , -14.898147 ,  18.076640 , &
            46.107317 ,  35.946352 ,  22.562267 , -39.569415 , &
           -43.438223 , -55.682381 , -24.349207 ,  52.572509 , &
            24.507541 ,  42.198873 ,  19.850975 , -41.799412 , &
            -8.533322 , -15.663322 ,  -7.978735 ,  16.990026 , &
             1.167545 ,   2.428273 ,   0.263832 ,  -1.568924 , &
            -0.265666 ,   0.900522 ,  -0.383530 ,   0.001026 , &
             2.749892 ,  -8.419774 ,   3.438275 ,   0.023706 , &
           -11.703235 ,  32.348946 , -12.384558 ,  -0.317123 , &
            26.431399 , -65.580982 ,  22.630073 ,   1.394674 , &
           -34.091088 ,  74.808992 , -21.434981 ,  -2.977644 , &
            25.042637 , -47.316301 ,   9.069436 ,   3.350396 , &
            -9.527931 ,  15.164369 ,  -0.767674 ,  -2.028522 , &
             1.571371 ,  -1.719171 ,   0.038905 ,   0.746250 /

data d51ik /-0.005802 ,  -0.198322 ,  -0.440286 ,  -1.750137 , &
             0.009889 ,   0.345448 ,   0.684101 ,   3.849764 , &
             0.077015 ,   0.373514 ,   1.699031 ,  15.256571 /

data c51ijk/-0.074991 ,  -1.066106 ,   1.258173 ,  -0.115450 , &
             0.771284 ,   9.589939 , -11.338655 ,   1.052139 , &
            -3.146666 , -35.791135 ,  42.309353 ,  -3.980068 , &
             6.498031 ,  71.610615 , -84.405279 ,   8.102800 , &
            -7.093954 , -82.653536 ,  96.855079 ,  -9.639142 , &
             3.811333 ,  54.733867 , -63.658308 ,   6.778503 , &
            -0.759248 , -19.069505 ,  21.963283 ,  -2.587078 , &
            -0.000101 ,   2.648971 ,  -2.946827 ,   0.315838 , &
             0.573212 ,   1.375637 ,  -2.625496 ,   0.801207 , &
            -5.282765 , -12.481586 ,  23.730554 ,  -7.249380 , &
            19.919157 ,  47.064320 , -88.810794 ,  27.189005 , &
           -39.359187 , -95.305948 , 177.687296 , -54.645464 , &
            43.267625 , 111.464788 ,-204.470563 ,  63.510348 , &
           -25.869751 , -74.654150 , 134.767958 , -42.807243 , &
             7.602800 ,  26.034831 , -46.593826 ,  15.480762 , &
            -0.855869 ,  -3.564916 ,   6.259163 ,  -2.077193 , &
            -0.623759 ,   0.127532 ,   0.797016 ,  -0.328101 , &
             5.754569 ,  -1.033602 ,  -7.391167 ,   2.955915 , &
           -22.019647 ,   3.378640 ,  28.436028 , -11.006391 , &
            45.160053 ,  -5.612633 , -58.660422 ,  21.835991 , &
           -53.528305 ,   4.694200 ,  70.041519 , -24.787007 , &
            36.699034 ,  -1.306354 , -48.688538 ,  16.130152 , &
           -13.376344 ,  -0.536333 ,  18.422713 ,  -5.955664 , &
             2.066210 ,   0.463133 ,  -2.706017 ,   1.431161 /

data d52ik /-0.002965 ,  -0.045327 ,  -0.073402 ,  -0.393744 , &
             0.011573 ,   0.144236 ,   0.262384 ,   1.662631 , &
             0.072494 ,   0.421731 ,   1.753864 ,  16.087311 /

data c52ijk/ -0.030185 ,  -0.125040 ,   0.280741 ,  -0.144880 , &
              0.271846 ,   1.142654 ,  -2.524148 ,   1.301160 , &
             -0.996758 ,  -4.313443 ,   9.369023 ,  -4.835003 , &
              1.893863 ,   8.681252 , -18.521680 ,   9.598821 , &
             -1.959008 , -10.004766 ,  20.961671 , -10.960490 , &
              1.060943 ,   6.539091 , -13.518233 ,   7.193341 , &
             -0.266817 ,  -2.197653 ,   4.560319 ,  -2.515905 , &
              0.023091 ,   0.284247 ,  -0.600422 ,   0.346696 , &
              0.243135 ,   0.312006 ,  -0.949575 ,   0.516194 , &
             -2.137348 ,  -2.959072 ,   8.582667 ,  -4.637560 , &
              7.653834 ,  11.628540 , -32.049222 ,  17.237544 , &
            -14.265800 , -24.433711 ,  63.793541 , -34.228685 , &
             14.649410 ,  29.439613 , -72.766689 ,  39.105248 , &
             -8.090036 , -20.062026 ,  47.367690 , -25.727239 , &
              2.190301 ,   6.946991 , -16.145885 ,   9.068609 , &
             -0.229987 ,  -0.915464 ,   2.144804 ,  -1.253010 , &
             -0.338488 ,   0.250096 ,   0.098527 ,  -0.013658 , &
              3.108590 ,  -2.108831 ,  -1.057788 ,   0.095075 , &
            -11.904232 ,   7.336728 ,   4.614787 ,  -0.199995 , &
             24.670834 , -13.555506 , -10.650267 ,  -0.076810 , &
            -30.045036 ,  14.070604 ,  14.231053 ,   0.939441 , &
             21.669709 ,  -7.703702 , -11.428345 ,  -1.364690 , &
             -8.456276 ,   1.679655 ,   5.377736 ,   0.385316 , &
              1.417136 ,   0.178406 ,  -0.938062 ,   0.576120 /

data d105ik / 0.060686 ,  -1.041861 ,   0.260359 ,  -3.612278 , &
             -0.031126 ,   0.833667 ,  -0.158378 ,   6.622774 , &
              0.084301 ,   0.386560 ,   2.060522 ,  16.716146 /

data c105ijk/ -0.909657 ,   4.296794 ,  -3.746867 ,   0.815741 , &
               8.637419 , -35.979317 ,  30.230812 ,  -6.586534 , &
             -33.429591 , 122.704169 , -98.380670 ,  21.241418 , &
              67.000407 ,-218.686791 , 164.864805 , -34.243695 , &
             -72.539891 , 217.267368 ,-149.859338 ,  27.071711 , &
              39.962788 ,-119.734546 ,  71.253715 ,  -6.632197 , &
              -9.309476 ,  35.029262 , -16.918071 ,  -2.485306 , &
               0.696841 ,  -4.560086 ,   2.576325 ,   0.071107 , &
               1.760717 ,  -3.968878 ,   2.945659 ,  -0.682742 , &
             -16.074113 ,  33.714187 , -24.089843 ,   5.651252 , &
              60.090827 ,-117.015837 ,  79.929260 , -18.918699 , &
            -117.625088 , 213.283419 ,-137.882765 ,  32.423882 , &
             127.560653 ,-218.505599 , 131.160159 , -28.897851 , &
             -74.585899 , 125.838802 , -67.164311 ,  10.749282 , &
              21.268034 , -39.062852 ,  17.844960 ,   0.254896 , &
              -2.445269 ,   5.293862 ,  -2.758516 ,   0.208695 , &
              -0.821507 ,   1.068790 ,  -0.386203 ,   0.081014 , &
               7.357867 ,  -8.965258 ,   2.792418 ,  -0.669905 , &
             -27.150088 ,  30.563729 ,  -7.561958 ,   2.277838 , &
              53.159329 , -54.186722 ,   8.727294 ,  -4.132201 , &
             -59.264879 ,  52.860964 ,  -1.757563 ,   4.307290 , &
              37.460477 , -27.460053 ,  -5.181011 ,  -2.427853 , &
             -12.324019 ,   6.741839 ,   4.291999 ,   0.267030 , &
               1.710346 ,  -0.411902 ,  -0.693280 ,   0.546410 /

data d101ik / -0.009646,  -0.007814 ,  -0.072726 ,  -2.130523 , &
              0.022721 ,  -0.038210 ,   0.284140 ,   4.700758 , &
              0.074961 ,   0.563987 ,   1.922534 ,  17.306716 /

data c101ijk/-0.760853 ,   1.066721 ,  -0.657585 ,   0.093095 , &
              6.602811 ,  -8.925041 ,   5.343182 ,  -0.763699 , &
            -23.666550 ,  30.895458 , -17.972107 ,   2.622282 , &
             45.150600 , -57.239488 ,  32.480355 ,  -4.909423 , &
            -49.201235 ,  61.206427 , -34.094475 ,   5.389770 , &
             30.503419 , -37.798273 ,  20.693001 ,  -3.330156 , &
            -10.003940 ,  12.573294 ,  -6.684299 ,   0.948994 , &
              1.360284 ,  -1.707883 ,   0.883381 ,  -0.131783 , &
              1.715170 ,  -2.000487 ,   1.089512 ,  -0.177509 , &
            -14.788359 ,  16.653052 ,  -8.800853 ,   1.494798 , &
             52.454482 , -57.153659 ,  29.408760 ,  -5.306318 , &
            -98.621945 , 104.621899 , -52.861865 ,  10.349833 , &
            105.574534 ,-110.338720 ,  55.380547 , -11.906655 , &
            -64.305196 ,  67.371247 , -33.644218 ,   7.709517 , &
             20.849815 , -22.362893 ,  10.846993 ,  -2.272727 , &
             -2.839756 ,   3.053501 ,  -1.427663 ,   0.322070 , &
             -0.835934 ,   0.892113 ,  -0.230450 ,   0.009059 , &
              7.223641 ,  -7.198260 ,   1.369831 ,  -0.047386 , &
            -25.772675 ,  23.584818 ,  -2.403849 ,   0.126432 , &
             49.120209 , -40.217788 ,  -0.687044 ,  -0.428745 , &
            -54.106484 ,  37.792759 ,   7.191027 ,   1.232177 , &
             34.684967 , -18.710344 ,  -9.300879 ,  -1.733481 , &
            -11.941293 ,   4.005852 ,   5.232539 ,   0.672266 , &
              1.741729 ,  -0.004771 ,  -0.935471 ,   0.540445 /


data d102ik/  0.002665 ,   0.064218 ,   0.073152 ,   0.041387 , &
             -0.004594 ,  -0.201348 ,  -0.094342 ,   0.199615 , &
              0.089965 ,   0.655093 ,   2.155139 ,  19.635947 /

data c102ijk/ -0.078113 ,   0.156871 ,  -0.113234 ,   0.018625 , &
               0.728704 ,  -1.354971 ,   0.931663 ,  -0.150687 , &
              -2.783381 ,   4.819754 ,  -3.137677 ,   0.498896 , &
               5.652734 ,  -9.165955 ,   5.604287 ,  -0.877760 , &
              -6.645659 ,  10.115214 ,  -5.722679 ,   0.883319 , &
               4.600539 ,  -6.517691 ,   3.284553 ,  -0.481465 , &
              -1.762495 ,   2.267722 ,  -0.933177 ,   0.099119 , &
               0.286706 ,  -0.310722 ,   0.090825 ,  -0.012669 , &
               0.304736 ,  -0.417701 ,   0.251263 ,  -0.043449 , &
              -2.696167 ,   3.437977 ,  -1.924215 ,   0.335944 , &
               9.731093 , -11.544717 ,   5.917116 ,  -1.063790 , &
             -18.583432 ,  20.596517 ,  -9.476378 ,   1.821269 , &
              20.510439 , -21.431666 ,   8.576329 ,  -1.849758 , &
             -13.483895 ,  13.373863 ,  -4.265942 ,   1.019481 , &
               5.058790 ,  -4.718817 ,   0.913895 ,  -0.138257 , &
              -0.834315 ,   0.667710 ,  -0.025564 ,   0.020289 , &
              -0.108240 ,   0.219176 ,   0.063450 ,  -0.050531 , &
               1.005557 ,  -1.553255 ,  -1.095288 ,   0.498456 , &
              -3.932456 ,   4.051580 ,   6.253365 ,  -1.992710 , &
               8.579563 ,  -4.265938 , -17.196464 ,   4.068156 , &
             -11.597965 ,  -0.023082 ,  25.623447 ,  -4.318269 , &
               9.766545 ,   4.006457 , -21.270707 ,   2.107865 , &
              -4.391714 ,  -3.332652 ,   9.414515 ,  -0.612328 , &
               0.809866 ,   0.983858 ,  -1.545014 ,   0.723112 /

data d155ik /-0.143614 ,  -1.802803 ,  -3.283256 , -23.817563 , &
              0.134556 ,   1.287954 ,   3.009716 ,  22.338812 , &
              0.057303 ,   0.338638 ,   1.530256 ,  14.918505 /

data c155ijk/-1.754607 ,  -1.224826 ,   2.370708 ,   0.021255 , &
             16.519883 ,  10.985880 , -21.316703 ,  -0.394066 , &
            -63.030064 , -42.191925 ,  79.815220 ,   2.261013 , &
            124.491516 ,  88.765832 ,-159.318733 ,  -6.140006 , &
           -134.533741 ,-107.534331 , 179.857018 ,   8.489466 , &
             77.291876 ,  71.369057 ,-112.011508 ,  -5.095180 , &
            -21.438325 , -21.630346 ,  33.820255 ,   0.452434 , &
              2.354487 ,   1.896572 ,  -3.225453 ,  -0.157336 , &
              2.464773 ,   0.307103 ,  -1.673627 ,  -0.138626 , &
            -22.058418 ,  -3.500796 ,  15.403436 ,   1.390185 , &
             80.443980 ,  16.811722 , -58.704026 ,  -5.847156 , &
           -153.044066 , -42.311522 , 118.334016 ,  13.298540 , &
            161.368298 ,  57.932896 ,-133.681130 , -17.096309 , &
            -92.681882 , -40.724862 ,  82.747054 ,  11.220146 , &
             26.836505 ,  11.830914 , -24.831266 ,  -2.619183 , &
             -3.220318 ,  -0.785983 ,   2.363577 ,   0.429849 , &
             -0.584246 ,  -0.360568 ,   0.945748 ,  -0.108381 , &
              5.083391 ,   3.755755 ,  -8.868441 ,   0.920781 , &
            -18.078771 , -16.185004 ,  34.423441 ,  -3.118916 , &
             33.812807 ,  37.096918 , -71.110167 ,   5.224266 , &
            -35.726107 , -48.411005 ,  83.634770 ,  -4.199651 , &
             21.370262 ,  35.624625 , -55.587800 ,   1.130176 , &
             -6.687048 , -13.452362 ,  19.272417 ,  -0.131511 , &
              0.909929 ,   2.108001 ,  -2.477977 ,   0.620463 /

data d151ik /-0.013837 ,   0.035234 ,  -0.404531 ,  -2.250395 , &
              0.029591 ,  -0.101124 ,   0.896260 ,   4.171519 , &
              0.077341 ,   0.573668 ,   1.867302 ,  18.610360 /

data c151ijk /-0.712182 ,   0.864633 ,  -0.550187 ,   0.175554 , &
               6.332043 ,  -7.441526 ,   4.676306 ,  -1.550435 , &
             -22.801996 ,  25.913119 , -16.157049 ,   5.654211 , &
              42.779277 , -47.010162 ,  29.311847 , -11.023906 , &
             -44.815957 ,  47.692148 , -30.036153 ,  12.362220 , &
              26.110538 , -26.951350 ,  17.187416 ,  -7.828978 , &
              -7.873737 ,   7.858573 ,  -4.931180 ,   2.441523 , &
               0.967371 ,  -0.886526 ,   0.515495 ,  -0.303053 , &
               1.165242 ,  -1.047146 ,   0.594778 ,  -0.277382 , &
             -10.148375 ,   8.741385 ,  -4.884761 ,   2.453530 , &
              35.687119 , -29.226159 ,  16.189656 ,  -8.998575 , &
             -65.128921 ,  50.331911 , -28.077437 ,  17.759164 , &
              66.188118 , -48.020499 ,  27.585616 , -20.310899 , &
             -37.576141 ,  25.532949 , -15.105655 ,  13.127280 , &
              11.275263 ,  -7.149336 ,   3.925405 ,  -4.059032 , &
              -1.433045 ,   0.767334 ,  -0.312377 ,   0.513723 , &
              -0.195087 ,  -0.205808 ,   0.541769 ,  -0.077578 , &
               1.675330 ,   2.241516 ,  -5.222595 ,   0.678200 , &
              -5.757358 , -10.192325 ,  20.969668 ,  -2.391506 , &
              10.283294 ,  24.719199 , -45.062086 ,   4.214929 , &
             -10.565463 , -34.240927 ,  55.474690 ,  -3.560545 , &
               6.612726 ,  27.075821 , -38.961176 ,   0.860058 , &
              -2.297575 , -11.334466 ,  14.581940 ,   0.091141 , &
               0.363071 ,   2.027117 ,  -2.075237 ,   0.614955 /

data d152ik/  -0.003428 ,   0.001297 ,  -0.071853 ,  -0.094278 , &
               0.012911 ,  -0.005177 ,   0.292608 ,   0.232440 , &
               0.083612 ,   0.511658 ,   2.138276 ,  20.393321 /

data c152ijk/ -0.242884 ,   0.183512 ,  -0.035515 ,   0.013459 , &
               2.016718 ,  -1.511894 ,   0.282032 ,  -0.118957 , &
              -6.855786 ,   5.112449 ,  -0.936982 ,   0.443676 , &
              12.307149 ,  -9.184796 ,   1.722146 ,  -0.909287 , &
             -12.549651 ,   9.472921 ,  -1.920979 ,   1.100246 , &
               7.264954 ,  -5.613106 ,   1.265140 ,  -0.754659 , &
              -2.212853 ,   1.754457 ,  -0.409911 ,   0.238694 , &
               0.270865 ,  -0.212977 ,   0.048743 ,  -0.032812 , &
               0.683112 ,  -0.326348 ,  -0.088549 ,  -0.017898 , &
              -5.502306 ,   2.431700 ,   0.906492 ,   0.161849 , &
              18.003358 ,  -7.239387 ,  -3.630656 ,  -0.651964 , &
             -30.833607 ,  11.159928 ,   7.258833 ,   1.528255 , &
              29.789186 ,  -9.829099 ,  -7.659542 ,  -2.179601 , &
             -16.391264 ,   5.240072 ,   4.264466 ,   1.708188 , &
               4.853715 ,  -1.636828 ,  -1.284500 ,  -0.531257 , &
              -0.596433 ,   0.200934 ,   0.157561 ,   0.081072 , &
              -0.182255 ,  -0.245485 ,   0.710423 ,  -0.174968 , &
               1.344585 ,   2.621568 ,  -6.619574 ,   1.538403 , &
              -4.019807 , -11.378428 ,  25.569912 ,  -5.527582 , &
               6.460108 ,  26.065817 , -52.808655 ,  10.331220 , &
              -6.432837 , -34.213100 ,  62.604675 , -10.429869 , &
               4.273433 ,  26.030454 , -42.409021 ,   5.204830 , &
              -1.536910 , -10.742858 ,  15.270577 ,  -1.233805 , &
               0.222966 ,   1.919968 ,  -2.078423 ,   0.777365 /
contains

subroutine init1(config)
    use rad_spectral_gray, only : rad_spectral_gray_init
    implicit none
    type(configuration), intent(in) :: config
    modelName_ = config%SpectralModelName
    ! parse options
    select case (modelName_)
        case ("GRAY")
            modelId_ = gray
            call rad_spectral_gray_init(config)
        case("GRAY-CONST")
            modelId_ = graycon
            call rad_spectral_gray_init(config)
        case("GRAY-WSGG")
            modelId_ = graywsgg
            call rad_spectral_gray_init(config)
        CASE("NONGRAY-WSGG")
            modelId_ = nongraywsgg
            call rad_spectral_gray_init(config)
        case default
            modelId_ = unknown
    end select
    call rad_gas_species_init(config)
    call rad_solid_species_init(config)
end subroutine

subroutine finalize1()
    use rad_spectral_gray, only : rad_spectral_gray_final
    implicit none
    select case (modelId_)
        case(gray)
            call rad_spectral_gray_final()
        case(graycon) 
            call rad_spectral_gray_final()
        case(graywsgg)
            call rad_spectral_gray_final()
        CASE(nongraywsgg)
            call rad_spectral_gray_final()
        case default
    end select
end subroutine

subroutine calculateSpectralFields()
    use functions, only : fluid_at,wall_at, funijk
    use discretelement
    use mfix_pic, only: des_stat_wt, mppic, epg_p
    use functions, only: IS_NORMAL
    USE compar 
    use des_thermo, only : DES_T_s
    use fldvar, only: ep_g 
    use geometry, only : vol, DX, DY, DZ 
    use run, only: UNITS
    use indices, only: i_of, j_of, k_of
    use functions, only: funijk
    USE compar
    implicit none
    include 'usrnlst.inc'
    type(gasInfoType) :: gasinfo
    type(solidParInfoType), dimension(mphas) :: solidinfos
    real(dp),dimension(nq) :: kg          ! gas absorption coefficient cm^-1
    real(dp),dimension(mphas,nq) :: ks    ! solid absorption coefficient cm^-1
    real(dp),dimension(mphas,nq) :: scats ! solid scattering coefficient cm^-1
    real(dp),dimension(0:mphas,nq) :: emiss !spectral emission

    integer :: ijk, m, i, j, k, NP, lM, iq
    real(dp) :: p_conversion,cellvol,areabvol,sbcon
    real(dp), dimension(nq) :: ai,ki
    real(dp) :: tmp, rmrk,spl,ttr


    !---M.S-rad: initialize the gasinfo 
    gasinfo%T = 0.
    gasinfo%P = 0.
    gasinfo%C(:) = 0.
    gasinfo%volFrac = 0.
    !---M.S-rad: initialize the solidsinfos
    solidinfos%volFrac = 0.
    solidinfos%radius = 0.
    solidinfos%T = 0.
    solidinfos%totSpecies = 0
    do m=1, mphas
        solidinfos(m)%speciesId(:) = 0
        solidinfos(m)%speciesMassFrac(:) = 0.
    end do
    !---M.S-rad: Determine the pressure conversion from the current MFIX unit to bar
    if (UNITS=='CGS') then
        p_conversion = 1.01325e-6  ! 1atm = 1.01325e+6 barye 
    else if (UNITS=='SI') then
        p_conversion = 1.01325e-5  ! 1atm = 1.01325e+5 Pa 
    else 
        print*, "Radiative property calculation: unknown units: ", UNITS
        print*, "Known units are CGS and SI"
    endif

    sbcon = sigma  ! SI system
    if(UNITS=='CGS') sbcon = sigma * 1.d3

    ! gas phase   (gray medium, WSGG model)  
     if(modelId_.eq.3) then
 
      do ijk = 1, dimension_3
          if (fluid_at(ijk)) then
            call getGasInfo(ijk, gasinfo)
            !---M.S-rad: convert pressure to bar
            gasinfo%P = gasinfo%P * p_conversion
            rmrk       = gasinfo%C(H2O)/gasinfo%C(CO2)
            if(gasinfo%C(H2O).eq.0) rmrk = 0.0d0 
            spl        = rad_spl
            if(UNITS=='CGS') spl = spl*1.e-2 ! cm to m  
            tmp        = 0.0d0 
            ai         = 0.0d0
            ki         = 0.0d0 
            ttr        = gasinfo%T
            if(rmrk.le.0.0d0.or.gasinfo%C(CO2).eq.0.0d0) then 
               tmp = 0.0d0 
            else if(rad_wsgg.eq.'SMITH82') then
                if (rmrk.eq.2) then
                  do i=1,3
                    ki(i) = ki2(i)
                    ai(i) = bei2(i,1)*1.0D-1 + bei2(i,2)*1.0D-4*ttr+ &
                    bei2(i,3)*1.0D-7*ttr**2.0 +bei2(i,4)*1.0D-11*ttr**3.0
                    tmp= tmp + ai(i)*(1.0d0-exp(-ki(i)*(gasinfo%P * 0.986923)*(gasinfo%C(H2O)+gasinfo%C(CO2))*spl))
                  end do
                else 
                  do i=1,3
                    ki(i) = ki1(i)
                    ai(i) = bei1(i,1)*1.0D-1 + bei1(i,2)*1.0D-4*ttr+ &
                    bei1(i,3)*1.0D-7*ttr**2.0 +bei1(i,4)*1.0D-11*ttr**3.0
                    tmp= tmp + ai(i)*(1.0d0-exp(-ki(i)*(gasinfo%P * 0.986923)*(gasinfo%C(H2O)+gasinfo%C(CO2))*spl))
                  end do
                end if
            else if(rad_wsgg.eq.'BORDBAR14') then
                do i=1,4
                   ai(i)  = 0.0d0 
                   ki(i) = 0.0d0
                   do k=1,5   ! absorption co-efficient 
                      ki(i) = ki(i) + dik(i,k) * (rmrk)**(k-1)
                      do j=1,5 ! weighting factors 
                        ai(i)= ai(i)+cijk(i,j,k)*(rmrk**(k-1))*(ttr/1200.0)**(j-1)
                      end do
                   end do
                tmp= tmp + ai(i)*(1.0d0-exp(-ki(i)*(gasinfo%P * 0.986923)*(gasinfo%C(H2O)+gasinfo%C(CO2))*spl)) 
                end do
            else if(rad_wsgg.eq.'JOHANSSON11') then
                do i=1,4
                  ai(i) = 0.0d0
                  ki(i) = k1j(i) +k2j(i) *rmrk
                  do k=1,3
                    ai(i)=ai(i)+(c1ik(i,k)+c2ik(i,k)*rmrk+c3ik(i,k)*rmrk*rmrk)*(ttr/1200.0)**(k-1)
                  end do
                  tmp= tmp+ai(i)*(1.0d0-exp(-ki(i)*(gasinfo%P * 0.986923)*(gasinfo%C(H2O)+gasinfo%C(CO2))*spl))
               end do
            else if (rad_wsgg.eq.'TAYLOR74') then  
              if(rmrk.eq.2) then 
                 do i=1,4
                    ai(i) = bf21i(i) + bf22i(i) * ttr* 1.0d-5
                    tmp   = tmp + ai(i)*(1.0-exp(-kf2i(i)*(gasinfo%P * 0.986923)*(gasinfo%C(H2O)+gasinfo%C(CO2))*spl))
                    ki(i) = kf2i(i)
                 end do
              else 
                 do i=1,4
                    ai(i) = bf11i(i) + bf12i(i) * ttr* 1.0D-5
                    tmp   = tmp + ai(i)*(1.0d0-exp(-kf1i(i)*(gasinfo%P * 0.986923)*(gasinfo%C(H2O)+gasinfo%C(CO2))*spl))
                    ki(i) = kf1i(i) 
                 end do
              end if

            else if(rad_wsgg.eq.'DORIGON13') then 
              if (rmrk.eq.2) then  ! 400 K to 2500 K , 0.001 atmm to 10 atm m 
                do i=1,4 
                  ai(i) = b2i0(i) + b2i1(i)*1.0e-5*ttr  + b2i2(i)*1.0e-8*ttr*ttr + b2i3(i)*1.0e-11*ttr**3.0 + b2i4(i) * 1.0E-15*ttr**4.0
                  tmp = tmp + ai(i)*(1.0-exp(-kpi2(i)*(gasinfo%P * 0.986923)*(gasinfo%C(H2O)+gasinfo%C(CO2))*spl))
                  ki(i) = kpi2(i)
               end do 
              else   ! rmrk.eq.1 
               do i=1,4
                  ai(i) = b1i0(i) + b1i1(i)*1.0e-5*ttr  + b1i2(i)*1.0e-8*ttr*ttr + b1i3(i)*1.0e-11*ttr**3.0 + b1i4(i) * 1.0E-15*ttr**4.0
                  tmp = tmp + ai(i)*(1.0-exp(-kpi1(i)*(gasinfo%P * 0.986923)*(gasinfo%C(H2O)+gasinfo%C(CO2))*spl))
                  ki(i) = kpi1(i)
               end do 
              end if
            end if  ! rad_wsgg loop end
            kg(1) =  -log(1.0d0 - tmp) / spl
            k_g(ijk,:) = kg(1)   
            E(ijk,0,:) = sbcon*(ttr**4.0)/Pi  
            if(UNITS=='CGS') k_g(ijk,:) = k_g(ijk,:) *1.0D-2 ! meters to centi meters
      end if     ! fluid_ijk loop  
     end do      ! do loop 
     else if (modelId_.eq.4) then
      do ijk = 1, dimension_3
          if (fluid_at(ijk)) then
          call getGasInfo(ijk, gasinfo)
          !---M.S-rad: convert pressure to bar
          gasinfo%P = gasinfo%P * p_conversion
          rmrk       = gasinfo%C(H2O)/gasinfo%C(CO2)
          if(gasinfo%C(H2O).eq.0) rmrk = 0.0d0 
          spl        = rad_spl
          if(UNITS=='CGS') spl = spl*1.e-2 ! cm to m 
          ttr        = gasinfo%T
            if(gasinfo%C(H2O).eq.0.0.or.gasinfo%C(CO2).eq.0) goto 999
            if(rad_wsgg.eq.'SMITH82') then
                if (rmrk.eq.2) then
                  do i=1,nq  
                    ki(i) = ki2(i)
                    ai(i) = bei2(i,1)*1.0D-1 + bei2(i,2)*1.0D-4*ttr+ &
                    bei2(i,3)*1.0D-7*ttr**2.0 +bei2(i,4)*1.0D-11*ttr**3.0
                  end do
                else 
                  do i=1,nq  
                    ki(i) = ki1(i)
                    ai(i) = bei1(i,1)*1.0D-1 + bei1(i,2)*1.0D-4*ttr+ &
                    bei1(i,3)*1.0D-7*ttr**2.0 +bei1(i,4)*1.0D-11*ttr**3.0
                  end do
                end if
            else if(rad_wsgg.eq.'BORDBAR14') then
                do i=1,nq 
                   ki(i) = 0.0d0 
                   ai(i) = 0.0d0
                   do k=1,5   ! absorption co-efficient 
                      ki(i) = ki(i) + dik(i,k) * (rmrk)**(k-1)
                      do j=1,5 ! weighting factors 
                        ai(i)= ai(i)+cijk(i,j,k)*(rmrk**(k-1))*(ttr/1200.0)**(j-1)
                      end do
                   end do
                end do
            else if(rad_wsgg.eq.'SHAN17') then
               if(gasp.le.5.0) then
                  if(rmrk.le.0.5) then
                   do i=1,nq
                     ai(i) = 0.0d0
                     ki(i) = d55ik(i,1)*rmrk*rmrk + d55ik(i,2)*rmrk + d55ik(i,3)
                     do j=1,8
                       ai(i)=ai(i)+(c55ijk(i,j,1)*rmrk*rmrk+c55ijk(i,j,2)*rmrk+c55ijk(i,j,3))*(ttr/1200.0)**(8-j)
                     end do
                   end do 
                  else if(rmrk.gt.0.5.and.rmrk.le.1.0) then
                   do i=1,nq
                    ai(i) = 0.0d0
                    ki(i) = d51ik(i,1)*rmrk*rmrk + d51ik(i,2)*rmrk + d51ik(i,3)
                    do j=1,8
                      ai(i)=ai(i)+(c51ijk(i,j,1)*rmrk*rmrk+c51ijk(i,j,2)*rmrk+c51ijk(i,j,3))*(ttr/1200.0)**(8-j)
                    end do
                   end do 
                  else
                   do i=1,nq
                    ai(i) = 0.0d0
                    ki(i) = d52ik(i,1)*rmrk*rmrk + d52ik(i,2)*rmrk + d52ik(i,3)
                    do j=1,8
                      ai(i)=ai(i)+(c52ijk(i,j,1)*rmrk*rmrk+c52ijk(i,j,2)*rmrk+c52ijk(i,j,3))*(ttr/1200.0)**(8-j)
                    end do
                   end do 
                  end if
                else if((gasp.gt.5.0).and.(gasp.le.10.0)) then
                  if(rmrk.le.0.5) then
                  do i=1,nq 
                     ai(i) = 0.0d0
                     ki(i) = d105ik(i,1)*rmrk*rmrk + d105ik(i,2)*rmrk + d105ik(i,3)
                     do j=1,8
                       ai(i)=ai(i)+(c105ijk(i,j,1)*rmrk*rmrk+c105ijk(i,j,2)*rmrk+c105ijk(i,j,3))*(ttr/1200.0)**(8-j)
                     end do
                  end do 
                  else if(rmrk.gt.0.5.and.rmrk.le.1.0) then
                   do i=1,nq 
                     ai(i) = 0.0d0
                     ki(i) = d101ik(i,1)*rmrk*rmrk + d101ik(i,2)*rmrk + d101ik(i,3)
                     do j=1,8
                       ai(i)=ai(i)+(c101ijk(i,j,1)*rmrk*rmrk+c101ijk(i,j,2)*rmrk+c101ijk(i,j,3))*(ttr/1200.0)**(8-j)
                     end do
                   end do 
                  else
                   do i=1,nq 
                     ai(i) = 0.0d0 
                     ki(i) = d102ik(i,1)*rmrk*rmrk + d102ik(i,2)*rmrk + d102ik(i,3)
                     do j=1,8
                       ai(i)=ai(i)+(c102ijk(i,j,1)*rmrk*rmrk+c102ijk(i,j,2)*rmrk+c102ijk(i,j,3))*(ttr/1200.0)**(8-j)
                     end do
                   end do 
                  end if
                else
                  if(rmrk.le.0.5) then
                   do i=1,nq 
                     ai(i) = 0.0d0 
                     ki(i) =d155ik(i,1)*rmrk*rmrk + d155ik(i,2)*rmrk + d155ik(i,3)
                     do j=1,8
                       ai(i)=ai(i)+(c155ijk(i,j,1)*rmrk*rmrk+c155ijk(i,j,2)*rmrk+c155ijk(i,j,3))*(ttr/1200.0)**(8-j)
                     end do
                   end do 
                  else if(rmrk.gt.0.5.and.rmrk.le.1.0) then
                   do i=1,nq 
                     ki(i) = d151ik(i,1)*rmrk*rmrk + d151ik(i,2)*rmrk + d151ik(i,3)
                     do j=1,8
                       ai(i)=ai(i)+(c151ijk(i,j,1)*rmrk*rmrk+c151ijk(i,j,2)*rmrk+c151ijk(i,j,3))*(ttr/1200.0)**(8-j)
                     end do
                   end do 
                  else
                   do i=1,nq 
                     ai(i) = 0.0d0 
                     ki(i) = d152ik(i,1)*rmrk*rmrk + d152ik(i,2)*rmrk + d152ik(i,3)
                     do j=1,8
                       ai(i)=ai(i)+(c152ijk(i,j,1)*rmrk*rmrk+c152ijk(i,j,2)*rmrk+c152ijk(i,j,3))*(ttr/1200.0)**(8-j)
                     end do
                   end do 
                  end if
                end if
            else if(rad_wsgg.eq.'JOHANSSON11') then
                do i=1,nq 
                  ai(i) = 0.0d0 
                  ki(i) = k1j(i) +k2j(i) *rmrk
                  do k=1,3
                    ai(i)=ai(i)+(c1ik(i,k)+c2ik(i,k)*rmrk+c3ik(i,k)*rmrk*rmrk)*(ttr/1200.0)**(k-1)
                  end do
              end do
            else if (rad_wsgg.eq.'TAYLOR74') then  
              if(rmrk.eq.2) then 
                 do i=1,nq
                    ai(i) = bf21i(i) + bf22i(i) * ttr* 1.0d-5
                    ki(i) = kf2i(i)
                 end do
              else
                 do i=1,nq
                    ai(i) = bf11i(i) + bf12i(i) * ttr* 1.0d-5
                    ki(i) = kf1i(i)
                 end do
              end if  
            else if(rad_wsgg.eq.'DORIGON13') then 
              if (rmrk.eq.2) then  ! 400 K to 2500 K , 0.001 atmm to 10 atm m 
               do i=1,nq  
                 ai(i) = b2i0(i) + b2i1(i)*1.0e-5*ttr  + b2i2(i)*1.0e-8*ttr*ttr + b2i3(i)*1.0e-11*ttr**3.0 + b2i4(i) * 1.0E-15*ttr**4.0
                 ki(i) = kpi2(i)
               end do 
              else                 ! rmrk.eq.1 
               do i=1,nq 
                  ai(i) = b1i0(i) + b1i1(i)*1.0e-5*ttr  + b1i2(i)*1.0e-8*ttr*ttr + b1i3(i)*1.0e-11*ttr**3.0 + b1i4(i) * 1.0E-15*ttr**4.0
                  ki(i) = kpi1(i)
               end do 
              end if
            end if  ! rad_wsgg loop end
999	   continue 
           do i=1,nq 
              k_g(ijk,i) = ki(i) 
              a_g(ijk,i) = ai(i)
              E(ijk,0,i) = ai(i)*sbcon*(ttr**4.0)/Pi  
              if(UNITS=='CGS') k_g(ijk,i) = k_g(ijk,i) *1.0d-2 ! 1/meters to 1/cm 
          end do ! nq loop end
      end if     ! fluid_ijk loop  end  
     end do      ! do loop ijk end 

     else 
       do ijk = 1, dimension_3
          if (fluid_at(ijk)) then
              call getGasInfo(ijk, gasinfo)
              !---M.S-rad: convert pressure to bar
              gasinfo%P = gasinfo%P * p_conversion
              ! find absorption coeffient and emission term 
              call spectralgasCalc(gasinfo, kg, emiss) 
              ! put cell calculations into fields
              k_g(ijk,:) = kg
              E(ijk,0,:) = emiss(0,:)
          endif
       enddo
    end if 


    if(discrete_element.or.mppic) then
      if(modelId_.eq.2) then 
      ! gray constant DEM 
        do np=1,max_pip
          if(is_normal(np)) then
           lM = pijk(np,5)
             ijk = pijk(np,4)
               if(fluid_at(ijk)) then
                   cellvol        = vol(ijk)    
                   areabvol       = Pi*(des_radius(np)**2.0)/cellvol
                   k_s(ijk,lm,:)  = k_s(ijk,lm,:) + const_emisp*areabvol
                   E(ijk,lM,:)    = (sbcon/Pi)*(DES_T_s(NP))**4.d0
                   scat(ijk,lM,:) = scat(ijk,lM,:)+ (1.0d0-const_sfp)* (1.0d0-const_emisp)*areabvol
               end if! fluid_at loop ends here
          endif      ! normal   loop ends here
        end do       ! particle loop ends here
      else 
     ! gray  and non-gray DEM/PIC models 
         do ijk = ijkstart3, ijkend3
          if(fluid_at(ijk)) then
              if(pinc(ijk) > 0) then
                   do m=1, mphas
                      call getDesSolidInfo(ijk,m, solidinfos(m))
                   end do
                   call spectralsldCalc(solidinfos, ks, scats, emiss)
                   do m=1,mphas
                      k_s(ijk,m,:)  = ks(m,:)
                      E(ijk,m,:)    = emiss(m,:)
                      if(modelId_.eq.4) E(ijk,m,:) = E(ijk,m,:)*1.0/float(nq) ! Non-gray model 
                      scat(ijk,m,:) = scats(m,:)
                   end do
              end if
          end if
         end do
      end if 
   ! solid phase in TFM (gray and non-gray models) 
    else 
      do ijk = 1, dimension_3
        if (fluid_at(ijk)) then
            call getGasInfo(ijk, gasinfo)
            !---M.S-rad: convert pressure to bar
            gasinfo%P = gasinfo%P * p_conversion
              do m=1, smax
                 call getSolidInfo(ijk,m, solidinfos(m))
              end do
              call spectralsldCalc(solidinfos, ks, scats, emiss)
            ! put cell calculations into fields
             do m=1, smax
                   E(ijk,m,:)    = emiss(m,:) 
                   if(modelId_.eq.4) E(ijk,m,:) = E(ijk,m,:)*1.0/float(nq) ! Non-gray model 
                   k_s(ijk,m,:)  = ks(m,:)
                   scat(ijk,m,:) = scats(m,:)
             end do
        endif
      enddo
    end if 

    ! unit conversions
    if (UNITS=='CGS') then
        scat = scat/1.d2
    endif

    call extrapolation()

end subroutine

subroutine calcRadSources()
    use discretelement, only: discrete_element
    use mfix_pic,       only: mppic
    use fldvar,         only: T_g
    use indices,        only: i_of, j_of, k_of
    use functions,      only: fluid_at, wall_at
    use functions,      only: funijk
    use geometry,       only: dx,dy,dz
    USE compar

    implicit none
    type(gasInfoType) :: gasinfo
    real(dp),dimension(1) :: T
    integer iq, m,ijk,i,j,k
    double precision :: sbcon
    real(dp), dimension(dimension_3)  :: tdis
    include 'usrnlst.inc'
    radiationOn = RAD_ON 
    Srad = 0.d0

    sbcon = sigma  ! SI system
    if(UNITS=='CGS') sbcon = sigma * 1.d3

   if(modelId_.eq.3) then
     do ijk = 1, dimension_3
       do iq = 1, nq
          Srad(ijk,0) = Srad(ijk,0)+k_g(ijk,iq)*(G(ijk,iq)-4.d0*pi*E(ijk,0,iq))
          do m=1, mphas
             Srad(ijk,m) = Srad(ijk,m)+k_s(ijk,m,iq)*(G(ijk,iq)-4.d0*pi*E(ijk,m,iq))
          end do
       end do
      end do 
   else if(modelId_.eq.4) then 
     do ijk = 1, dimension_3
       do iq = 1, nq
          Srad(ijk,0) = Srad(ijk,0)+ k_g(ijk,iq)*(G(ijk,iq)-a_g(ijk,iq)*4.d0*sbcon*T_g(ijk)**4.0)  
          do m=1, mphas
             Srad(ijk,m) = Srad(ijk,m)+k_s(ijk,m,iq)*(G(ijk,iq)-4.d0*pi*E(ijk,m,iq))
          end do
       end do
      end do
   else 
       do iq = 1, nq
        Srad(:,0) = Srad(:,0)+k_g(:,iq)*(G(:,iq)-4.d0*pi*E(:,0,iq)) 
    end do
   end if 

   if(discrete_element.or.mppic) then
       do iq = 1, nq
           do m=1, mphas
               Srad(:,m) = Srad(:,m)+k_s(:,m,iq)*(G(:,iq)-4.d0*pi*E(:,m,iq))
           end do
       end do 
   else
       do iq = 1,nq   
           do m=1, smax
               Srad(:,m) = Srad(:,m)+k_s(:,m,iq)*(G(:,iq)-4.d0*pi*E(:,m,iq))
           end do
       end do
   end if 

      ! by far Srad units are W/(m^3) in SI or erg/(cm^3.s) in CGS
      ! C_pg uses J/kg.K in SI or cal/g.K in CGS
      ! 1 erg = 2.39005736d-8 cal
      if (UNITS=='CGS') Srad = Srad * 2.39005736d-8
    
end subroutine


subroutine spectralsldCalc(solidinfos, ks, scats, emiss)
    use rad_spectral_gray, only :rad_spectral_graysld_calc,rad_spectral_graysldconst_calc
    implicit none
    type(solidParInfoType),dimension(mphas), intent(in) :: solidinfos
    real(dp),dimension(mphas,nq), intent(out) :: ks
        ! solid absorption coefficient cm^-1
    real(dp),dimension(mphas,nq), intent(out) :: scats
        ! solid scattering coefficient cm^-1
    real(dp),dimension(0:mphas,nq), intent(out) :: emiss
        !spectral emission 

    select case(modelId_)
        case (gray)
            call rad_spectral_graysld_calc(solidinfos, ks, scats, emiss)
        case(graycon)
            call rad_spectral_graysldconst_calc(solidinfos, ks, scats, emiss)
        case(graywsgg)
            call rad_spectral_graysld_calc(solidinfos, ks, scats, emiss)
        case(nongraywsgg)
            call rad_spectral_graysld_calc(solidinfos, ks, scats, emiss)
        case default
    end select
end subroutine

subroutine spectralgasCalc(gasinfo, kg, emiss)
    use rad_spectral_gray, only :rad_spectral_graygas_calc,rad_spectral_graygasconst_calc
    implicit none
    type(gasInfoType), intent(in) :: gasinfo
    real(dp),dimension(nq), intent(out) :: kg
        ! gas absorption coefficient cm^-1
        ! gas emission coefficient cm^-1
    real(dp),dimension(0:mphas,nq), intent(out) :: emiss
        !spectral emission 

    select case(modelId_)
        case (gray)
            call rad_spectral_graygas_calc(gasinfo, kg, emiss)
        case(graycon)
            call rad_spectral_graygasconst_calc(gasinfo, kg, emiss)
        case default
    end select
end subroutine


subroutine extrapolation()
    use rad_util, only : extrapolateGhostCells
    implicit none
    integer :: iq, is
    do iq=1, nq
        call extrapolateGhostCells(k_g(:,iq))
        do is = 1, size(k_s, 3)
            call extrapolateGhostCells(k_s(:,is,iq))
            call extrapolateGhostCells(scat(:,is,iq))
        end do
    end do

end subroutine

end module rad_spectral
