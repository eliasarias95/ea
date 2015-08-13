#!/bin/sh

#run java to get new synthetic
#/users/elias.arias/Home/git/ea/bench/src/pgs/data/bin2mdvd par=constant2D.par
#/users/elias.arias/Home/git/ea/bench/src/pgs/data/bin2mdvd par=complex2D.par
#/users/elias.arias/Home/git/ea/bench/src/pgs/data/bin2mdvd par=complex3D.par

#run slopexy
#/users/elias.arias/Home/git/ea/bench/src/pgs/data/slopexy par=slopexy.par

#convert slopexy output to binary
#/users/elias.arias/Home/git/ea/bench/src/pgs/data/mdvd2bin par=constant2D_sxy.par
#/users/elias.arias/Home/git/ea/bench/src/pgs/data/mdvd2bin par=complex2D_sxy.par
#/users/elias.arias/Home/git/ea/bench/src/pgs/data/mdvd2bin par=complex3Dx_sxy.par
#/users/elias.arias/Home/git/ea/bench/src/pgs/data/mdvd2bin par=complex3Dy_sxy.par
#/users/elias.arias/Home/git/ea/bench/src/pgs/data/mdvd2bin par=Santos_2D_sxy.par
#/users/elias.arias/Home/git/ea/bench/src/pgs/data/mdvd2bin par=Santos_2D.par
#/users/elias.arias/Home/git/ea/bench/src/pgs/data/mdvd2bin par=Triton_2D_sxy.par
#/users/elias.arias/Home/git/ea/bench/src/pgs/data/mdvd2bin par=Triton_2D.par

#run sdw
sdw -a par=sdw.par

#convert sdw output to binary
#/users/elias.arias/Home/git/ea/bench/src/pgs/data/mdvd2bin par=constant2D_sdw.par
#/users/elias.arias/Home/git/ea/bench/src/pgs/data/mdvd2bin par=complex2D_sdw.par
#/users/elias.arias/Home/git/ea/bench/src/pgs/data/mdvd2bin par=complex3Dx_sdw.par
#/users/elias.arias/Home/git/ea/bench/src/pgs/data/mdvd2bin par=complex3Dy_sdw.par
/users/elias.arias/Home/git/ea/bench/src/pgs/data/mdvd2bin par=Santos_2D_sdw.par
#/users/elias.arias/Home/git/ea/bench/src/pgs/data/mdvd2bin par=Triton_2D_sdw.par
sh jy ../slopeDemo2D.py
