# parallelBlast
BLAST is a widely used sequence alignment software. With the increasing of avaliable sequence data, we are dealing with larger data with BlAST. But BLAST is always ran under single thread model, some times it takes ages to finish the alignment.
ParallelBlast could split the input files into smaller files, run BLAST for them simultaneously, and finally merge the results together, whith could take advantage of the current multiple-core-CPU.
You are encouraged to import the source code into eclipse to modify and complie.
