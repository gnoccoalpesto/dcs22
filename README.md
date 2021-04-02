# PROJECT 3

distributed task assignment for robotic networks

# TASK 1

1)N agents

solve constraint coupled linear program problem (1):

min{z1...zN} sum{1,N}(ci'zi)
subj: 	sum{1,N} Hizi=b
	zi € Pi, i€{1...N}

ci € R^ni, Hi € R^Sxni, b € R^S
compact polyhedron Pi={zi € R^ni | Dizi=<di, Gizi=gi } i€{1...N}

using Distributed Dual Subgradient:

contained in "Distributed optimization for smart cyber-physical networks"

2) montecarlo simulation and plots to show convergence varying N and problem size
------------------------------------------------------------------------
# TASK 2
{-> possible further improvements}

N robots, N tasks scattered in a finite environment

each robot has to serve a task, each task must be served by at most one robot

min(tot robots' TRAVEL* distance) *->compenetration?

bipartite assignment graph Ga={Va,Ua;Ea}
tasks Va={1...N}, agents Ua={1...N}
edge (i,k) € Ea exist iff agent i assignable to tank k

-> many ramndomly select robots that can't execute task

cost cik of ende (i,k) [note cik=0 iff edge not exist]
1) euclidean distance
->2) manhattan? non compenetration (each slab occupied by a robot)
3) 

linear optimization problem:

min{x}sum{(i,k) € Ea} (cik xik)
subj: 	0=<x=<1 LB=0, UB=1
	sum{k|(i,k) € Ea}(xik=1}, i€{1...N} NeighIn(k)
	sum{i|(i,k) € Ea}(xik=1}, k€{1...N} NeighOut(i)

x=[...;...,xik,...;...]
X=[xt=1,...]

recast in form 1) by: 
zi=[...xik,...]   (==xi)

ci=[...,cik...]
ci=ci_average+noise(ci) random stochastic disturbance


k|(i,k) € Ea

(at most k=1...N)

Pi= {xi |0=<xi=<1, sum{k|(i,k)€ Ea}(xik=1), i€{1...N}
Pi=

Hi, b vector defined suitably
1) Hi=I
2) Hi=...?

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
xik=1 for one i, one k, which is the task k assigned eventually to agent i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
(tasks assigned by evaluating ci vectors)


ALGORITHM:
1)spawn robots
2) spawn task
3) computation

posiotions can be random, but solution must be unique
 (avoid rng(time) to have a nice solution)
