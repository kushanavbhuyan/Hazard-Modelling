#! --matrixtable --lddout


############################################
# Model:                                   #
# Date:                                    #
# Version: 1.0                             #
# Author:                                  #
############################################


binding







initial


#dem
report grad.map = sin(atan(slope(dem.map)));

report dem_filter.map = windowaverage(dem.map,200);
report factor_noise.map =2.0* max(0.0,-0.5 +(1.0 - grad.map));
report dem_fixed.map = dem.map * (1.0-factor_noise.map) + dem_filter.map * factor_noise.map;

report ldd_all.map = lddcreate(srtm.map,1e10,1e10,1e10,1e10);
report accuflux_all.map = accuflux(ldd.map,1.0);

report ldd.map = lddcreate(dem.map,1e10,1e10,1e10,1e10);
report accuflux.map = accuflux(ldd.map,1.0);
report pits.map = pit(ldd.map);
report catchments.map = catchment(ldd.map,pits.map);


report outpoint.map = pit(ldd.map);
report outpoint_act.map = pit(ldd.map);

report id.map = dem.map * 0.0 + 1.0;
#channel

channelmask1 = scalar(if(accuflux.map gt 13000,1.0,0.0));
channelldd1 = lddcreate(-accuflux.map * if(channelmask1 gt 0.5,1.0),1e31,1e31,1e31,1e31);

report channelaccuflux.map = accuflux(channelldd1,1.0);
report channelmask.map = scalar(if(accuflux.map gt 13000 ,1.0,0.0));
report channelldd.map = lddcreate(-accuflux.map * if(channelmask.map gt 0.5,1.0),1e31,1e31,1e31,1e31);

report chandepth.map = channelmask.map*0.015*9.68 * accuflux.map ** 0.32;
report channelwidth.map = channelmask.map*11.3 * accuflux.map ** 0.083;
report chanside.map = dem.map * 0.0;
report channelgrad.map = grad.map;
report chanman.map = dem.map * 0.0 + 0.075;
#land use
report cover.map = max(0.01,1.0-exp(-2.0 * ndvi.map/(1-max(ndvi.map,0.02))));
report lai.map = max(0.01,ln(1.0- cover.map)/(-0.4));


############################
#SURFACE MAPS 
############################


report n.map = lookupscalar(lu.tbl,1,luclass.map); # %
report rr.map = lookupscalar(lu.tbl,2,luclass.map);  # %
report ch.map = lookupscalar(lu.tbl,3,luclass.map);  # %


#soil stuff

report soildepth1.map = sd.map * 1000.0;



 # sl2 is 5cm depth
  S = sand.map /10.0; # %
  C = clay.map /10.0;  # %
  OM = organic.map;  # g/kg
  Gravel = gravel.map/10.0; #%
 
 # prep data
  S = S/100.0;
  C = C/100.0;
  OM = OM /100*1.72; #conversion org carbon to org matter 
  Gravel = Gravel/100;
  Densityfactor = 0.9;
 

  #output 
  report sand1.map = S;
  report clay1.map = C;
  report orgmat1.map = OM;
  report grav1.map = Gravel;


  # multiple regression eq
  M1500 =-0.024*S+0.487*C+0.006*OM+0.005*S*OM-0.013*C*OM+0.068*S*C+0.031; #W18) 
  M1500adj =M1500+0.14*M1500-0.02; #X18) 
  M33  =-0.251*S+0.195*C+0.011*OM+0.006*S*OM-0.027*C*OM+0.452*S*C+0.299; #Y18) 
  M33adj = M33+(1.283*M33*M33-0.374*M33-0.015); #Z18) 
  PM33    = 0.278*S+0.034*C+0.022*OM-0.018*S*OM-0.027*C*OM-0.584*S*C+0.078; #AA18)
  PM33adj = PM33+(0.636*PM33-0.107); #AB18)
  SatPM33 = M33adj + PM33adj; #AC18)
  SatSadj = -0.097*S+0.043; #AD18)
  SadjSat = SatPM33  + SatSadj; #AE18)
  Dens_om = (1-SadjSat)*2.65; #AF18)
  Dens_comp = Dens_om * Densityfactor; #AG18)
  PORE_comp =(1-Dens_om/2.65)-(1-Dens_comp/2.65); #AI18)
  M33comp = M33adj - 0.2*PORE_comp; #AJ18)
    	
  #output 
  report thetas1.map = cover(1-(Dens_comp/2.65),0.5); #AH18)
  PoreMcomp = thetas1.map-M33comp; #AK18)
  LAMBDA = (ln(M33comp)-ln(M1500adj))/(ln(1500)-ln(33)); #AL18)
  GravelRedKsat =(1-Gravel)/(1-Gravel*(1-1.5*(Dens_comp/2.65))); #AM18)
  report Ksat1.map =cover(1930*(PoreMcomp)**(3-LAMBDA)*GravelRedKsat,10.0); #AN18)
  report BD1.map = Gravel*2.65+(1-Gravel)*Dens_comp;    #U18
  report WP1.map = M1500adj;
  report FC1.map = M33adj;
  report PAW1.map = (M33adj - M1500adj)*(1-Gravel);

  bB = (ln(1500)-ln(33))/((ln(FC1.map)-ln(WP1.map)));
  aA = exp(ln(33) + bB*ln(FC1.map));

  report psi1.map = cover(max(10, aA*(FC1.map + 0.7 * (thetas1.map - FC1.map))**-bB),0.0);
  report thetai1.map = (FC1.map + 0.7 * (thetas1.map - FC1.map));

  report zero.map = dem.map * 0.0;
  report soildensity.map = zero.map + 2200;
  report soilifa.map = zero.map + 0.4;
  report soilrocksize.map = zero.map + 0.05;


  report coh.map = zero.map + 20.0;
  report aggrstab.map = zero.map + 12;
  report d50.map = zero.map + 60;
  report d90.map = zero.map + 90;



