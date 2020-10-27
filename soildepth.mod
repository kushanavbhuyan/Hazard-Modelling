############################################
# Model: SOIL DEPTH steady state           #
# Date:                                    #
# Version: 1.0                             #
# Author: Bastian van den Bout             #
############################################


binding

DEMO = Dem.map;
DEM = Dem.map;

timer
1 2 1;

initial

#dzx1 = max(0.0,DEM - cover(shift(DEM,0,-1),DEM));
#dzx2 = max(0.0,DEM - cover(shift(DEM,0,1),DEM));
#dzy1 = max(0.0,DEM - cover(shift(DEM,-1,0),DEM));
#dzy2 = max(0.0,DEM - cover(shift(DEM,1,0),DEM));

#dzt = max(0.001,dzx1 + dzx2 + dzy1 + dzy2);

#wx1 = dzx1/dzt;
#wx2 = dzx2/dzt;
#wy1 = dzy1/dzt;
#wy2 = dzy2/dzt;

sdx1 = scalar(0.0);
sdx2 = scalar(0.0);
sdy1 = scalar(0.0);
sdy2 = scalar(0.0);

dzx1 = scalar(0.0);
dzx2 = scalar(0.0);
dzy1 = scalar(0.0);
dzy2 = scalar(0.0);

dzt = scalar(0.0);

wx1 = scalar(0.0);
wx2 = scalar(0.0);
wy1 = scalar(0.0);
wy2 = scalar(0.0);

w_sum_in = scalar(0.0);
w_sum_out = scalar(0.0);

added = scalar(0.0);
store = scalar(0.0);
storeold = scalar(0.0);

R = scalar(2.0*2.0);

sd_prev = scalar(1.0);
sd = scalar(0.0);
error = scalar(0.0);
prod = scalar(0.0);
mov = scalar(0.0);
acc = scalar(0.0);

cellsize = 10.0;

QS = 0.5;
k = 0.001;
sdh = 1;


AH = cellsize * cellsize;
LH = cellsize;
LG = cellsize;


loops = 0;
loops_sd = 0;

dynamic 


repeat{   
DEM = DEMO + 0.00001 * store;

dzx1 = max(0.0,DEM - cover(shift(DEM,0,-1),DEM));
dzx2 = max(0.0,DEM - cover(shift(DEM,0,1),DEM));
dzy1 = max(0.0,DEM - cover(shift(DEM,-1,0),DEM));
dzy2 = max(0.0,DEM - cover(shift(DEM,1,0),DEM));

dzt = max(0.001,dzx1 + dzx2 + dzy1 + dzy2);

wx1 = dzx1/dzt;
wx2 = dzx2/dzt;
wy1 = dzy1/dzt;
wy2 = dzy2/dzt;


store = store + R; 
storeold = store;

store = store*0.75;
store = store + 0.25*cover(shift(wx2,0,-1),0.0) * cover(shift(storeold,0,-1),0.0);
store = store + 0.25*cover(shift(wx1,0,1),0.0) * cover(shift(storeold,0,1),0.0);
store = store + 0.25*cover(shift(wy2,-1,0),0.0) * cover(shift(storeold,-1,0),0.0);
store = store + 0.25*cover(shift(wy1,1,0),0.0) * cover(shift(storeold,1,0),0.0);


loops = loops + 1;


} until loops gt 200;

report sd_store.map = store;





acc = LH*AH*(sqrt(store/mapmaximum(store)));

repeat{

#calculate error in soil depth

sdx1 = cover(shift(sd,0,-1),sd);
sdx2 = cover(shift(sd,0,1),sd);
sdy1 = cover(shift(sd,-1,0),sd);
sdy2 = cover(shift(sd,1,0),sd);

dzx1 = max(0.0,cover(shift(DEM,0,-1),DEM) - DEM);
dzx2 = max(0.0,cover(shift(DEM,0,1),DEM) - DEM);
dzy1 = max(0.0,cover(shift(DEM,-1,0),DEM) - DEM);
dzy2 = max(0.0,cover(shift(DEM,1,0),DEM) - DEM);

dzt = max(0.001,dzx1 + dzx2 + dzy1 + dzy2);

wx1_in = dzx1*sdx1*sdx1*sdx1;
wx2_in = dzx2*sdx2*sdx2*sdx2;
wy1_in = dzy1*sdy1*sdy1*sdy1;
wy2_in = dzy2*sdy2*sdy2*sdy2;

w_sum_in = wx1_in +wx2_in +wy1_in +wy2_in;

dzx1 = max(0.0,DEM - cover(shift(DEM,0,-1),DEM));
dzx2 = max(0.0,DEM - cover(shift(DEM,0,1),DEM));
dzy1 = max(0.0,DEM - cover(shift(DEM,-1,0),DEM));
dzy2 = max(0.0,DEM - cover(shift(DEM,1,0),DEM));

dzt = max(0.001,dzx1 + dzx2 + dzy1 + dzy2);

wx1_out = dzx1*sdx1*sdx1*sdx1;
wx2_out = dzx2*sdx2*sdx2*sdx2;
wy1_out = dzy1*sdy1*sdy1*sdy1;
wy2_out = dzy2*sdy2*sdy2*sdy2;

w_sum_out = wx1_out +wx2_out +wy1_out +wy2_out;

mov = LH*k*(w_sum_in-w_sum_out)/LG;

#adapt soil depth estimate

sd = sdh *max(0.01,ln(max(1.0,AH/(( mov+ QS* (acc))))));
report sd_mov.map = mov;
report sd_acc.map = acc;

loops_sd = loops_sd + 1;
} until loops_sd gt 100;

sd = mapmaximum(sd)-sd;
report sd.map = sd;