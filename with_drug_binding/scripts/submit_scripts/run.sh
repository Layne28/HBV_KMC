#!/bin/bash

##seed kappa_length kappa_phi kappa_theta phiCD phiAB gb0 dimer_conc dimer_binding_rate mu_conformation dgirr drugconc drugaffinity drugbindingrate gAA gDB gBC gCD gDC gother 
#./source/assemble 104788256 4200.000 40.000 800.000 0.24 0.48 -9.8 -11.5 0.02 -4.5 0.100 0.0 0.0 0.0 0.3 0.1 -0.1 0 0.0 -.95
seed=104788
kappa_length=4200.000
kappa_phi=40.000
kappa_theta=800.00
phiCD=0.0 
phiAB=0.48 
gb0=-9.8
dimer_conc=-11.5
dimer_binding_rate=0.0002
mu_conformation=-4.5
dgirr=0.100
drugconc=0.0
drugaffinity=0.0
drugbindingrate=0.0
gAA=0.3
gDB=0.1
gBC=-0.1
gCD=10
gDC=0.0
gother=-0.95

./source/assemble ${seed} ${kappa_length} ${kappa_theta} ${kappa_phi} ${phiCD} ${phiAB} ${gb0} ${dimer_conc} ${dimer_binding_rate} ${mu_conformation} ${dgirr} ${drugconc} ${drugaffinity} ${drugbindingrate} ${gAA} ${gDB} ${gBC} ${gCD} ${gDC} ${gother}
