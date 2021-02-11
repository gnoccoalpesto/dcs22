# PROJECT 3

distributed task assignment for robotic networks

# TASK 1

1)N agents

solve constraint coupled linear program problem (1):

min{z1...zN} sum{1,N}(ci'zi)
subj: 	sum{1,N} Hizi=b
	zi € Pi, i€{1...N}

ci € R^ni, Hi € R^Sxni, b € R^S
compact polyhedron Pi={zi € R^ni | Dizi=<di, Gizi=gi} i€{1...N}

using Distributed Dual Subgradient:

contained in "Distributed optimization for smart cyber-physical networks"

2) montecarlo simulation and plots to show convergence varying N and problem size

# TASK 2

N robots, N tasks scattered in a finite environment

each robot has to serve a task, each task must be served by at most one robot

min(tot robots' travel distance)

bipartite assignment graph Ga={Va,Ua;Ea}
tasks Va={1...N}, agents Ua={1...N}
edge (i,k) € Ea exist iff agent i assignable to tank k
cost cik of ende (i,k) [note cik=0 iff edge not exist]

linear optimization problem:

min{x}sum{(i,k) € Ea} (cik xik)
subj: 	0=<x=<1
	sum{k|(i,k) € Ea}(xik=1}, i€{1...N}
	sum{i|(i,k) € Ea}(xik=1}, k€{1...N}

x: stack(xic)

recaste in form 1) by: 
zi=stack(xic), ci=stack(cik), k|(i,k) € Ea
Pi= {xi |0=<xi=<1, sum{k|(i,k)€ Ea}(xik=1), i€{1...N}
Hi, b vector defined suitably

tasks assigned by evaluating ci vectors
note: tasks spawned contemporarly
