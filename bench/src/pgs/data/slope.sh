#!/bin/sh
#run java to get new synthetic
#/users/elias.arias/Home/git/ea/bench/src/pgs/data/bin2mdvd par=constant2D.par
/users/elias.arias/Home/git/ea/bench/src/pgs/data/bin2mdvd par=complex2D.par
/users/elias.arias/Home/git/ea/bench/src/pgs/data/slopexy par=slopexy.par
/users/elias.arias/Home/git/ea/bench/src/pgs/data/mdvd2bin par=tempx.par
/users/elias.arias/Home/git/ea/bench/src/pgs/data/mdvd2bin par=tempy.par
sdw -a par=sdw_test2D.par
#/users/elias.arias/Home/git/ea/bench/src/pgs/data/mdvd2bin par=constant2D_out.par
/users/elias.arias/Home/git/ea/bench/src/pgs/data/mdvd2bin par=complex2D_out.par
sh jy ../slopeDemo2D.py
