#!/bin/sh
# Plotting trajectories for all parameter sets, locked to the stimulus
# TC (VPM and POm) and RE populations

pra=tlp
prg=plp

fla=d1
flb=d2

grfl=ptrint.gr

tlp.ex $fla
plp.ex $fla
stm.ex $fla

cat > ptrint.xx.fsed <<EOF
s/scan=u/scan=e/
s/INPUT ff 2.0 5.0 8.0 11.0/INPUT ff parmin=1.0 parmax=11.5 npar=21/
EOF

sed -f ptrint.xx.fsed $pra.n.$fla > $pra.n.$flb

tlp.ex $flb
plp.ex $flb

awk '{if ($2 == "1") print $1, 1000.0 * $7}' $prg.res.$flb > ptrint.xx.ult.in
awk '{if ($2 == "1") print $1, 1000.0 * $6}' $prg.res.$flb > ptrint.xx.ult.t5
awk '{if ($2 == "2") print $1, 1000.0 * $7}' $prg.res.$flb > ptrint.xx.upt.in
awk '{if ($2 == "2") print $1, 1000.0 * $6}' $prg.res.$flb > ptrint.xx.upt.t5

awk '{if ($1 != "#") print $1, $2}' delay_exp > delay.exp.vpm.xx
awk '{if ($1 != "#") print $1, $3}' delay_exp > delay.exp.pom.xx

awk '{print $1, $2 * 40.0}' intg.res.a3 > intg.res.a3.xx
awk '{print $1, $2 * 40.0}' intg.res.a4 > intg.res.a4.xx

intnor.ex ptrint ult 40.0
intnor.ex ptrint upt 40.0

xmgrace \
        -graph 0 plp.xx.$fla.1.ult plp.xx.$fla.2.ult \
                 plp.xx.$fla.3.ult plp.xx.$fla.4.ult \
        -graph 1 plp.xx.$fla.1.upt plp.xx.$fla.2.upt \
                 plp.xx.$fla.3.upt plp.xx.$fla.4.upt \
        -graph 2 stm.flt.$fla \
        -graph 3 stm.fpt.$fla \
        -graph 4 intn.res.ult intg.res.a3.xx \
                 ptrint.xx.ult.t5 delay.exp.vpm.xx \
        -graph 5 intn.res.upt intg.res.a4.xx \
                 ptrint.xx.upt.t5 delay.exp.pom.xx \
        -graph 6 plp.xx.$fla.1.ura plp.xx.$fla.2.ura \
                 plp.xx.$fla.3.ura plp.xx.$fla.4.ura \
        -graph 7 plp.xx.$fla.1.urb plp.xx.$fla.2.urb \
                 plp.xx.$fla.3.urb plp.xx.$fla.4.urb \
        -hdevice EPS -p $grfl -printfile ptrint.$fla.eps

/bin/rm $pra.out.$fla $pra.tmp.$fla $pra.col.$fla
/bin/rm $pra.out.$flb $pra.tmp.$flb $pra.col.$flb
/bin/rm plp.xx.$fla.*.*
/bin/rm stm.flt.$fla stm.fpt.$fla stm.ffr.$fla
/bin/rm ptrint.xx.fsed
/bin/rm plp.xx.$flb.*.*
/bin/rm ptrint.xx.ult.in ptrint.xx.ult.t5 ptrint.xx.upt.in ptrint.xx.upt.t5
/bin/rm delay.exp.vpm.xx delay.exp.pom.xx
/bin/rm intg.res.a3.xx  intg.res.a4.xx
/bin/rm intn.out.ult intn.res.ult intn.out.upt intn.res.upt
