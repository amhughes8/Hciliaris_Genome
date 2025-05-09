# Running `blobtools` in the browser
This is a Markdown file showing how to run `blobtools` in a browser using one node

1. Login to Discovery
```
ssh hughes.annab@login.discovery.neu.edu
```

2. Start an interactive SLURM batch job via `srun` and activate the `conda` environment
```
# start interactive job
srun --partition=short --nodes=1 --cpus-per-task=1 --pty /bin/bash

# activate the conda environment
module load anaconda3/2022.05
source activate /work/gatins/hci_genome/env

# now change directory to where blobtools lives
cd /work/gatins/hci_genome/processing/blobtools2
```

3. Start a `screen` session for the API
```
# this starts a screen session named blob_api
screen -S blob_api # you will enter a blank terminal

# run the blobtoolkit-api
BTK_API_PORT=8880 BTK_PORT=8881 BTK_FILE_PATH=/work/gatins/hci_genome/processing/blobtools2/BlobDirs ./blobtoolkit-api

# now to exit this screen execute control + A + D
# you will return to your normal terminal
```

4. Start a `screen` session for the viewer
```
# this starts a screen session named blob_viewer
screen -S blob_viewer

# run the blobtoolkit-viewer
BTK_API_PORT=8880 BTK_PORT=8881 ./blobtoolkit-viewer

# exit the screen session with control + A + D
```

5. Use `ssh ProxyJump` to connect to the viewer port
```
# open a new terminal tab in iTerm2 (or whatever you use) and ssh tunnel by ProxyJump
ssh -J hughes.annab@login.discovery.neu.edu -L8880:localhost:8880 -L8881:localhost:8881 hughes.annab@<nodename>

# open you browser and type this in the search bar
localhost:8881

# now you can use blobtools!

# -J specifies the bastion node (from where the tunnel will jump) since Discovery is the main source
# -L8880:localhost:8880 points to the node specified at the end (localhost) and the port the API is running on
# -L8881:localhost:8881 points to the port the viewer is running on
```

6. How to rejoin and kill screen sessions when done
```
# when you finish with blobtools go in the terminal where you executed the ProxyJump and exit
exit
# and again
exit

# now in your working node terminal kill the viewer, then the API like so
screen -r blob_viewer # rejoin the viewer screen
control + c # kills the viewer
exit # this kills that screen session; you'll return to the main session

# join the API screen session
screen -r blob_api
control + c
exit
```
**And that's it! That's how you do it! Enjoy `blobbing`!**