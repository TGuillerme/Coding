##########################
#Absolute distance in a set of trees
##########################
#Measure the absolute distance between two trees based on their shape
#v1.0
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
=======

>>>>>>> FETCH_HEAD
=======
#Dummy update: Adam want's to know the truth
>>>>>>> parent of 2cf3386... Adam is no longer the truth
=======
#Dummy update: Adam want's to know the truth
>>>>>>> parent of a2cdcf9... Merge pull request #3 from kanead/patch-1
#----
#guillert(at)tcd.ie 25/07/2013
##########################

distrees<-function(t0,trees,method="topo"){ #add method="brlen" and options

warning('In developement')

	#Library
	require(ape)
	require(phangorn)
	#require(phyloclust)

	#Input
	if(class(t0) !="phylo") {
		stop("t0 is not a phylo object")
	}
	if(class(trees) !="multiPhylo") {
		stop("t0 is not a multiPhylo object")
	}
	if(class(method) !="character") {
		stop("method should be 'topo' or 'brlen'")
	}
	if(method=="topo") {ok<-"ok"} else {
		if (method=="brlen") {ok<-"ok"} else {
		stop("method should be 'topo' or 'brlen'")
			}
		}
#		##Binary
#		if(is.binary.tree(t0)) {ok<-"ok"} else {
#			warning("t0 is not a binary tree")}
#		test<-NULL
#		for (i in 1:length(trees)){test[i]<-is.binary.tree(trees[[i]])}
#		if(any(test, FALSE))
#			{warning("trees contain at least one non-binary tree")} else {ok<-"ok"}
#		
#		##Tip labels
#		test<-NULL
#		for (i in 1:length(trees)){test[i]<-(sort(t0$tip.label)==sort(trees[[i]]$tip.label))}
#		if(any(test, FALSE))
#			{stop("trees don't have the same tip labels")} else {ok<-"ok"}

	#Loading the functions
	distree.obs<-function(x,y){
		Symmetry.obs<-abs((sd(c(abs((balance(x)[,1]-1)-balance(x)[,2]), abs((balance(x)[1,1]-1)-balance(x)[1,2])))/mean(c(abs((balance(x)[,1]-1)-balance(x)[,2]), abs((balance(x)[1,1]-1)-balance(x)[1,2]))))-(sd(c(abs((balance(y)[,1]-1)-balance(y)[,2]), abs((balance(y)[1,1]-1)-balance(y)[1,2])))/mean(c(abs((balance(y)[,1]-1)-balance(y)[,2]), abs((balance(y)[1,1]-1)-balance(y)[1,2])))))
		Topology.obs<-dist.topo(x,y)
		return(c(Symmetry.obs, Topology.obs))}

	distree.max<-function(n){
		Symmetry.max<-(sum(((seq(from=(n-0), to=(n-(n-1)))-((n*(n/2)+(n/2)-1)/n))^2))/n)/((n*(n/2)+(n/2)-1)/n)
		Topology.max<-(2*(n-2)-2)
		return(c(Symmetry.max, Topology.max))}
	
	distree.abs<-function(x,y){
		n<-length(x$tip.label)
		obs<-distree.obs(x,y)
		max<-distree.max(n)
		abs<--1*(obs-max)/max
		return(abs)}	

	#Measuring the distance
	Symmetry<-NULL
	Topology<-NULL
	for (i in 1:length(trees)){
		dist<-distree.abs(t0,trees[[i]])
		Symmetry[i]<-dist[1]
		Topology[i]<-dist[2]}
		
	#Output
	output<-data.frame("Symmetry"=Symmetry, "Topology"=Topology)
	return(output)

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
=======
	
=======
	#Dummy
	Truth<-"adam"
>>>>>>> parent of 2cf3386... Adam is no longer the truth
=======
	#Dummy
	Truth<-"adam"
>>>>>>> parent of a2cdcf9... Merge pull request #3 from kanead/patch-1

>>>>>>> FETCH_HEAD
}




#Testing

Sym.max<-NULL
Top.max<-NULL
Sym.obs<-NULL
Top.obs<-NULL
Sym.abs<-NULL
Top.abs<-NULL

for (n in 5:100){Sym.max[n]<-distree.max(n)[[1]]}
for (n in 5:100){Top.max[n]<-distree.max(n)[[2]]}
for (n in 5:100){Sym.obs[n]<-distree.obs(rtree(n),rtree(n))[[1]]}
for (n in 5:100){Top.obs[n]<-distree.obs(rtree(n),rtree(n))[[2]]}
for (n in 5:100){Sym.abs[n]<-distree.abs(rtree(n),rtree(n))[[1]]}
for (n in 5:100){Top.abs[n]<-distree.abs(rtree(n),rtree(n))[[2]]}



plot(Sym.max, main="Symmetric distance variaton", xlab="n",ylab="Symetric distance")
points(Sym.obs, col="blue")
points(Sym.abs, col="red")
legend(0,15, c("Maximal","Random","Absolute"), col=c("black","blue","red"),pch=21)

plot(Top.max, main="Topological distance variaton", xlab="n",ylab="Topological distance")
points(Top.obs, col="blue")
points(Top.abs, col="red")
legend(0,150, c("Maximal","Random","Absolute"), col=c("black","blue","red"),pch=21)



trees.rand<-NULL
for (i in 1:99){trees.rand[[i]]<-rmtree(125,51)}

for (i in 1:length(trees.rand)){
	for (j in 1:length(trees.rand[[1]])){
		trees.rand[[i]][[j]]$edge.length<-NULL
	}
}

for (i in 1:length(trees.rand)){
	Sym.rand<-distrees(trees.rand[[i]][[1]], trees.rand[[i]])[,1]
	Topo.rand<-distrees(trees.rand[[i]][[1]], trees.rand[[i]])[,2]}
