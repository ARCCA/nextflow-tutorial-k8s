# nextflow-tutorial

This are initial notes on how to use Nextflow on a Kubernetes (k8s) cluster running
on Sparrow.
The workflow assumes that your pipeline is accessible through GitHub and any required
programs are available in a container.

This example uses a repo hosted in the ARCCA GitHub site and a basic Nextflow container
(https://hub.docker.com/r/nextflow/nextflow). If you need additonal tools, you need
to prepare a custom Docker container (the ARCCA team can help you).

## Installing dependencies
Nextflow+k8s relies on creating a container (with any software required by the pipeline)
by fetching an image from some place accessible by Nextflow (e.g. [DockerHub][docker-hub]

In this example, we create a Docker container in the user's laptop (it is assumed that
Docker is already available) using the Dockerfile provided in this repo. 

To build an image  using the command line:
```
docker build --tag munozcriollojj/nf-pipeline-test:latest -f Dockerfile .
```

In the above command the tag used assumes that the image created will be pushed to the
repo `nf-pipeline-test` located in the user's (munozcriollojj) Docker repo.


## Setting up the environment
Create a directory to work on, e.g:
```
username@cl2(hawk)$ mkdir /scratch/$USER/nf-work
username@cl2(hawk)$ cd /scratch/$USER/nf-work
```

Load required files:
```
module purge
module use /apps/local/modules/projects
module load nextflow/21.10.6
module load scwXXXX/kubernetes
module list
```

### Transferring data to and from the Kubernetes cluster
The k8s cluster has a storage volume attached. To access the data in this volume, we need
to create an interactive pod:
```
username@cl2(hawk)$ nextflow kuberun login -profile k8s
Pod started: nostalgic-nightingale
bash-4.2#
bash-4.2# pwd
/scw1162-data/hawk.username
bash-4.2# ls
nextflow.config  work
```

At this point you can explore the data that was produced by the nextflow pipeline.
**Note: the pod will be terminated when you disconnect** (i.e. when you type `exit`
or <kbd>Ctrl</kbd> + <kbd>d</kbd>).
```
$ bash-4.2# ls -l output/
total 4
lrwxrwxrwx 1 root root 89 May 17 11:51 all.markdup.genecount.txt -> /scw1162-data/c.c1045890/work/b4/eceec0dfe946c50859b51a384c76ed/all.markdup.genecount.txt
...
```

To copy data from Hawk to the k8s storage volume, we can use the pod created in the previous step (i.e. `nostalgic-nightingale`). Make sure to double check your username in the commmand below:
```
username@cl2(hawk)$ cat testData/nSEP01-00303-S01-WB-c-0hr_S15_L001_R1_001.fastq.gz | kubectl exec -i nostalgic-nightingale -- tee /scw1162-data/username/data/nSEP01-00303-S01-WB-c-0hr_S15_L001_R1_001.fastq.gz > /dev/null
```

Now, this limits the transfer to one file at a time. This can be automatised to copy all
the files in a directory at the same time (see `transfer-data.sh` for a bit of
inspiration on how to do this).

To copy data back to Hawk you need to have a running pod (see above how to start an
interactive pod). You need to use the name of the running pod (`nostalgic-nightingale` in
the previous example):
```
username@cl2(hawk)$ mkdir output
username@cl2(hawk)$ kubectl exec -i nostalgic-nightingale -- cat /scw1162-data/username/work/b4/eceec0dfe946c50859b51a384c76ed/all.markdup.genecount.txt > output/all.markdup.genecount.txt
username@cl2(hawk)$ ls output/
all.markdup.genecount.txt
```

### Running the tutorial pipeline
Clone ARCCA Nextflow turorial for k8s and change to the scw1162 example branch:
```
username@cl2(hawk)$ git clone https://github.com/ARCCA/nextflow-tutorial-k8s.git
username@cl2(hawk)$ git checkout scw1162-example
```

At this point you should be able to run the pipeline test:
```
username@cl2(hawk)$ nextflow kuberun -r scw1162-example -latest ARCCA/nextflow-tutorial-k8s -profile k8s -with-trace trace.txt
```

The above command launches a pod on the Kubernetes cluster running on Sparrow and executes
the example pipeline.



## Remove completed pods
To remove pods completed with an `Error` status:
```
username@cl2(hawk)$ kubectl delete pod --field-selector=status.phase==Failed
```

To remove pods completed with a `Completed` status:
```
username@cl2(hawk)$ kubectl delete pod --field-selector=status.phase==Succeeded
```


### Making changes (untested).
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




[dockerhub]: https://hub.docker.com/
