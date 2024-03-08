
build:

%.png : %.npy ; qeasy 0 1 $^ $@

build: x_sq3.png x_sq10.png

#x_sq3.png  : x.png ; ./krt square3 $^ | qeasy 0 1 - $@
#x_sq10.png : x.png ; ./krt square10 $^ | qeasy 0 1 - $@
#x_sq40.png : x.png ; ./krt square40 $^ | qeasy 0 1 - $@
x_sq%.png : x.png ; ./krt square$* $^ | qeasy 0 1 - $@


#%SCRIPT ./krt square3 f/x.png  |qeasy 0 1 - f/x_sq3.png
#%SCRIPT ./krt square10 f/x.png |qeasy 0 1 - f/x_sq10.png
#%SCRIPT ./krt square40 f/x.png |qeasy 0 1 - f/x_sq40.png
#%SCRIPT ./krt square5  f/adelson.png |qauto -p 0 - f/adelson_sq5.png
#%SCRIPT ./krt gauss5   f/adelson.png |qauto -p 0 - f/adelson_g5.png
#%SCRIPT ./krt landc25  f/adelson.png |qauto -p 0 - f/adelson_l25.png
#%SCRIPT ./krt square5 -h gap2 f/adelson.png|qauto -p 0 - f/adelson_sq5_gap2.png
#%SCRIPT ./krt gauss5  -h gap2 f/adelson.png|qauto -p 0 - f/adelson_g5_gap2.png
#%SCRIPT ./krt landc25 -h gap2 f/adelson.png|qauto -p 0 - f/adelson_l25_gap2.png
#%SCRIPT ./krt square5 -h gap10 f/adelson.png|qauto -p 0 - f/adelson_sq5_gap10.png
#%SCRIPT ./krt gauss5  -h gap10 f/adelson.png|qauto -p 0 - f/adelson_g5_gap10.png
#%SCRIPT ./krt landc25 -h gap10 f/adelson.png|qauto -p 0 - f/adelson_l25_gap10.png
#%SCRIPT ./krt square5 -h gap20 f/adelson.png|qauto -p 0 - f/adelson_sq5_gap20.png
#%SCRIPT ./krt gauss5  -h gap20 f/adelson.png|qauto -p 0 - f/adelson_g5_gap20.png
#%SCRIPT ./krt landc25 -h gap20 f/adelson.png|qauto -p 0 - f/adelson_l25_gap20.png
#%SCRIPT ./krt square5 -h gap200 f/adelson.png|qauto -p 0 - f/adelson_sq5_gap200.png
#%SCRIPT ./krt gauss5  -h gap200 f/adelson.png|qauto -p 0 - f/adelson_g5_gap200.png
#%SCRIPT ./krt landc25 -h gap200 f/adelson.png|qauto -p 0 - f/adelson_l25_gap200.png
#%SCRIPT ./krt square3  f/adelson.png |qauto -p 0 - f/adelson_sq3.png
#%SCRIPT ./krt square13  f/adelson.png |qauto -p 0 - f/adelson_sq13.png
#%SCRIPT ./krt square33  f/adelson.png |qauto -p 0 - f/adelson_sq33.png
#%SCRIPT ./krt gauss3   f/adelson.png |qauto -p 0 - f/adelson_g3.png
#%SCRIPT ./krt gauss13   f/adelson.png |qauto -p 0 - f/adelson_g13.png
#%SCRIPT ./krt gauss33   f/adelson.png |qauto -p 0 - f/adelson_g33.png
#%SCRIPT ./krt landc3   f/adelson.png |qauto -p 0 - f/adelson_l3.png
#%SCRIPT ./krt landc13   f/adelson.png |qauto -p 0 - f/adelson_l13.png
#%SCRIPT ./krt landc33   f/adelson.png |qauto -p 0 - f/adelson_l33.png
#%SCRIPT ./krt square3  -h gap20 f/adelson.png |qauto -p 0 -  f/adelson_sq3_gap20.png
#%SCRIPT ./krt square13 -h gap20  f/adelson.png |qauto -p 0 - f/adelson_sq13_gap20.png
#%SCRIPT ./krt square33 -h gap20  f/adelson.png |qauto -p 0 - f/adelson_sq33_gap20.png
#%SCRIPT ./krt gauss3   -h gap20 f/adelson.png |qauto -p 0 -  f/adelson_g3_gap20.png
#%SCRIPT ./krt gauss13  -h gap20  f/adelson.png |qauto -p 0 - f/adelson_g13_gap20.png
#%SCRIPT ./krt gauss33  -h gap20  f/adelson.png |qauto -p 0 - f/adelson_g33_gap20.png
#%SCRIPT ./krt landc3   -h gap20 f/adelson.png |qauto -p 0 -  f/adelson_l3_gap20.png
#%SCRIPT ./krt landc13  -h gap20  f/adelson.png |qauto -p 0 - f/adelson_l13_gap20.png
#%SCRIPT ./krt landc33  -h gap20  f/adelson.png |qauto -p 0 - f/adelson_l33_gap20.png
#%SCRIPT ./krt square13 -h gap2 f/adelson.png |qauto -p 0 - f/adelson_sq13_gap2.png
#%SCRIPT ./krt square13 -h logistic2 f/adelson.png |qauto -p 0 - f/adelson_sq13_logis2.png
#%SCRIPT ./krt square13 -h arctan2 f/adelson.png |qauto -p 0 - f/adelson_sq13_atan2.png
#%SCRIPT ./krt square13 -h gap10 f/adelson.png |qauto -p 0 - f/adelson_sq13_gap10.png
#%SCRIPT ./krt square13 -h logistic10 f/adelson.png |qauto -p 0 - f/adelson_sq13_logis10.png
#%SCRIPT ./krt square13 -h arctan10 f/adelson.png |qauto -p 0 - f/adelson_sq13_atan10.png
#%SCRIPT ./krt square13 -h gap20 f/adelson.png |qauto -p 0 - f/adelson_sq13_gap20.png
#%SCRIPT ./krt square13 -h logistic20 f/adelson.png |qauto -p 0 - f/adelson_sq13_logis20.png
#%SCRIPT ./krt square13 -h arctan20 f/adelson.png |qauto -p 0 - f/adelson_sq13_atan20.png
#%SCRIPT ./krt square13 -h gap200 f/adelson.png |qauto -p 0 - f/adelson_sq13_gap200.png
#%SCRIPT ./krt square13 -h logistic200 f/adelson.png |qauto -p 0 - f/adelson_sq13_logis200.png
#%SCRIPT ./krt square13 -h arctan200 f/adelson.png |qauto -p 0 - f/adelson_sq13_atan200.png
#%SCRIPT plambda zero:500x120 ":i 50 / floor" -o f/mach.npy
#%%SCRIPT plambda zero:500x120 ":i 50 / floor randg 0.01 * +" -o f/mach.npy
#%%SCRIPT plambda zero:500x120 ":i 50 / floor randg 0.1 * +" -o f/mach.npy
#%SCRIPT plambda f/mach.npy '7 - 14 * 127 +' -o f/mach.png
#%SCRIPT export X="set term pngcairo size 500,200;unset key;unset xtics"
#%SCRIPT cat f/mach.npy|(echo "$X";cline 0)|gnuplot >f/p_mach.png
#%SCRIPT ./krt square21 f/mach.npy|qeasy 0 1 - f/mach_sq21.png
#%SCRIPT ./krt square21 f/mach.npy|(echo "$X";cline 0)|gnuplot >f/p_mach_sq21.png
#%SCRIPT ./krt gauss7 f/mach.npy|qeasy 0 1 - f/mach_g7.png
#%SCRIPT ./krt gauss7 f/mach.npy|(echo "$X";cline 0)|gnuplot >f/p_mach_g7.png
#%SCRIPT ./krt cauchy3 f/mach.npy|qeasy 0 1 - f/mach_c3.png
#%SCRIPT ./krt cauchy3 f/mach.npy|(echo "$X";cline 0)|gnuplot >f/p_mach_c3.png
#%SCRIPT ./krt landc40 f/mach.npy|qeasy 0 1 - f/mach_la.png
#%SCRIPT ./krt landc40 f/mach.npy|(echo "$X";cline 0)|gnuplot >f/p_mach_la.png
#%SCRIPT plambda f/mach.npy ":i 0.004 * +" -o f/umach.npy
#%SCRIPT plambda f/umach.npy '7 - 14 * 127 +' -o f/umach.png
#%SCRIPT cat f/umach.npy|(echo "$X";cline 0)|gnuplot >f/p_umach.png
#%SCRIPT ./krt square21 f/umach.npy|qeasy 0 1 - f/umach_sq21.png
#%SCRIPT ./krt square21 f/umach.npy|(echo "$X";cline 0)|gnuplot >f/p_umach_sq21.png
#%SCRIPT ./krt gauss7 f/umach.npy|qeasy 0 1 - f/umach_g7.png
#%SCRIPT ./krt gauss7 f/umach.npy|(echo "$X";cline 0)|gnuplot >f/p_umach_g7.png
#%SCRIPT ./krt cauchy3 f/umach.npy|qeasy 0 1 - f/umach_c3.png
#%SCRIPT ./krt cauchy3 f/umach.npy|(echo "$X";cline 0)|gnuplot >f/p_umach_c3.png
#%SCRIPT ./krt landc40 f/umach.npy|qeasy 0 1 - f/umach_la.png
#%SCRIPT ./krt landc40 f/umach.npy|(echo "$X";cline 0)|gnuplot >f/p_umach_la.png
#%SCRIPT plambda f/mach.npy ":i 0.004 * -" -o f/dmach.npy
#%SCRIPT plambda f/dmach.npy '7 - 14 * 127 +' -o f/dmach.png
#%SCRIPT cat f/dmach.npy|(echo "$X";cline 0)|gnuplot >f/p_dmach.png
#%SCRIPT ./krt square21 f/dmach.npy|qeasy 0 1 - f/dmach_sq21.png
#%SCRIPT ./krt square21 f/dmach.npy|(echo "$X";cline 0)|gnuplot >f/p_dmach_sq21.png
#%SCRIPT ./krt gauss7 f/dmach.npy|qeasy 0 1 - f/dmach_g7.png
#%SCRIPT ./krt gauss7 f/dmach.npy|(echo "$X";cline 0)|gnuplot >f/p_dmach_g7.png
#%SCRIPT ./krt cauchy3 f/dmach.npy|qeasy 0 1 - f/dmach_c3.png
#%SCRIPT ./krt cauchy3 f/dmach.npy|(echo "$X";cline 0)|gnuplot >f/p_dmach_c3.png
#%SCRIPT ./krt landc40 f/dmach.npy|qeasy 0 1 - f/dmach_la.png
#%SCRIPT ./krt landc40 f/dmach.npy|(echo "$X";cline 0)|gnuplot >f/p_dmach_la.png
#%SCRIPT plambda zero:256x256 "randg" -o f/randg.npy
#%SCRIPT ./krt square7 f/randg.npy |plambda "48 * 1 + round"| ghisto -p | gnuplot > f/randg_k7_h.png
#%SCRIPT ./krt gauss1.5 f/randg.npy |plambda '1000 * round 1000 /'| ghisto -p | gnuplot > f/randg_g15_h.png
