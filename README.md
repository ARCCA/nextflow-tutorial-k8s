# nextflow-tutorial

This are initial notes on how to use Nextflow on a Kubernetes cluster running on Sparrow.
The workflow assumes that your pipeline is accessible through GitHub and any required
programs are available in a container.

This example uses a repo hosted in the ARCCA GitHub site and a basic Nextflow container
(https://hub.docker.com/r/nextflow/nextflow). If you need additonal tools, you need
to prepare a custom Docker container (the ARCCA team can help you).


## Running the tutorial pipeline

Create a directory to work on, e.g:
```
$ mkdir /scratch/$USER/nf-work
$ cd /scratch/$USER/nf-work
```

Load required files:
```
$ module purge
$ module use /apps/local/modules/projects
$ module load nextflow/21.10.6
$ module load scwXXXX/kubernetes
$ module list
```

Clone ARCCA Nextflow turorial for k8s:
```
$ git clone https://github.com/ARCCA/nextflow-tutorial-k8s.git
```

At this point you should be able to run the pipeline test:
```
$ nextflow kuberun -latest ARCCA/nextflow-tutorial-k8s -profile k8s -n nextflow
Pod started: lonely-perlman
N E X T F L O W  ~  version 21.10.6
Pulling ARCCA/nextflow-tutorial-k8s ...
 Already-up-to-date
Launching `ARCCA/nextflow-tutorial-k8s` [lonely-perlman] - revision: c32a5e1ac4 [master]
[e6/f1caf5] Submitted process > splitLetters
[77/82b6b1] Submitted process > convertToUpper (1)
[c3/c14fe2] Submitted process > convertToUpper (2)
HELLO
WORLD!
```

The above command launches a pod on the Kubernetes cluster running on Sparrow and executes
the tutorial pipeline.


# Accessing the data
To access the data, launch an interactive pod:
```
$ nextflow kuberun login -profile k8s
Pod started: festering-gutenberg
bash-4.2#
bash-4.2# pwd
/scw1162-data/hawk.username
bash-4.2# ls
nextflow.config  work
```

At this point you can explore the data that was produced by the nextflow pipeline.
Note: the pod will be terminated when you disconnect.
```
bash-4.2# ls work/e6/f1caf5513c2d317c64346dd1760f0d/
chunk_aa  chunk_ab
```

To copy data back to Hawk you need to have a running pod (see above how to start an
interactive pod). You need to use the name of the running pod (festering-gutenberg in
the previous example):
```
$ mkdir results
$ kubectl exec -i festering-gutenberg -n nextflow -- cat /scw1162-data/c.c1045890/work/e6/f1caf5513c2d317c64346dd1760f0d/chunk_aa > results/chunk_aa
$ kubectl exec -i festering-gutenberg -n nextflow -- cat /scw1162-data/c.c1045890/work/e6/f1caf5513c2d317c64346dd1760f0d/chunk_ab > results/chunk_ab
$ ls results
chunk_aa  chunk_ab
```

## Making changes (untested).
Assuming you have a GitHub account, create a fork of this repo. Go to 
https://github.com/ARCCA/nextflow-tutorial-k8s and click on 'fork', follow the 
instructions to create the fork in your own repo.

- clone your fork repo on Hawk:
  ```
  $ git clone your-repo
  $ cd your-repo
  ```

- make changes, e.g.:
  ```
  vim main.nf
  git add main.nf
  git commit -m 'explain the changes made'
  git push
  ```

- run kubernetes cluster pointing to your own repo:
  ```
  nextflow kuberun -r your-branch -latest your-repo -profile test,k8s
  ```

