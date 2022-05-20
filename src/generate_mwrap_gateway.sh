LLNAME=flattrirouts
MWNAME=${LLNAME}.mw
GATEWAY=${LLNAME}.c
../../mwrap/mwrap -c99complex -mex -${LLNAME} -c ${GATEWAY} ${MWNAME} 
../../mwrap/mwrap -c99complex -list -mex ${LLNAME} -mb ${MWNAME}  
