


# make a singulatiry definition file for pascalX
# definition file:
Bootstrap: docker
From: ubuntu:20.04

%labels
    Maintainer Alanna Cote
    Application PascalX
    Description "PascalX (gene & pathway enrichment) packaged in a reproducible SIF (Ubuntu 20.04, Py3.8)"

%help
    PascalX container image.
    After building the SIF, typical usage on HPC:
      apptainer exec pascalx_ubuntu20.04.sif pascalx -h

%environment
    # Runs at container runtime
    export LD_LIBRARY_PATH=/opt/PascalX/build/lib:$LD_LIBRARY_PATH
    export PATH=/usr/local/bin:/usr/bin:/bin

%post
    set -ex

    # Base OS deps
    apt-get update
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
        python3 python3-dev python3-setuptools python3-pip \
        g++ make libboost-all-dev wget unzip git ca-certificates \
        pkg-config libfreetype6-dev libpng-dev

    update-ca-certificates

    # Tell the runtime linker where PascalX libs will live
    echo "/opt/PascalX/build/lib" > /etc/ld.so.conf.d/pascalx.conf

    # Upgrade pip tooling first (prevents easy_install path)
    python3 -m pip install --upgrade pip setuptools wheel

    # Fetch PascalX (clean, shallow clone)
    cd /opt
    rm -rf PascalX
    git clone --depth 1 https://github.com/BergmannLab/PascalX.git
    cd PascalX

    # Build core libs/binaries
    make all
    ldconfig

    # Ensure build-time tests can find libruben.so
    export LD_LIBRARY_PATH=/opt/PascalX/build/lib:$LD_LIBRARY_PATH
    make test

    # ---- Python bindings (pin deps to Py3.8-compatible wheels) ----
    # Pre-install wheels that are known-good on Python 3.8
    python3 -m pip install --no-cache-dir "numpy<1.25" "matplotlib<3.8"

    # Now install the PascalX Python package with pip (NOT setup.py)
    cd python
    python3 -m pip install --no-cache-dir .

    # Optional: Jupyter for ad-hoc exploration (kept unpinned; can remove)
    python3 -m pip install --no-cache-dir jupyter

%runscript
    # Default entrypoint: run whatever the user passes, in a login shell
    exec /usr/bin/env bash -lc "$*"
    








2) Build the .sif remotely (no Linux needed)
Use Sylabs Remote Builder in the browser:
Go to the Remote Builder, create/login to your Sylabs account. 
cloud.sylabs.io
Upload Singularity.pascalx and start a build (it compiles on an x86_64 host—exactly what you need for Minerva). When it finishes, download the .sif. (CLI users can also singularity build --remote image.sif Singularity.pascalx from any Linux box with SingularityCE installed.) 
Sylabs
+1
Note: Apptainer itself does not support build --remote; that’s a SingularityCE/Sylabs feature. That’s why the web builder (or a Linux VM with SingularityCE) is the right choice here.
apptainer.org
3) Copy the image to Minerva and test
# from your Mac
scp pascalx_ubuntu20.04.sif <UNI>@minerva.hpc.mssm.edu:/sc/arion/projects/<proj>/containers/

# on Minerva
module load apptainer/1.3.6
apptainer exec /sc/arion/projects/<proj>/containers/pascalx_ubuntu20.04.sif pascalx -h
You should see the PascalX CLI help. Minerva supports Apptainer and documents using .sif images in batch jobs. It also auto-binds your project paths, so your data under /sc/arion/ is visible in the container. 
Icahn School of Medicine at Mount Sinai
4) Persist PascalX resources (faster subsequent runs)
Create a persistent ref/annotation dir once (e.g., /sc/arion/projects/<proj>/ref/pascalx/) so PascalX’s imported panels/cache live outside scratch and get reused. (Their docs describe one-time reference imports/caching.) 
Bergmann Lab
5) Wire into Nextflow on Minerva
Add these bits to nextflow.config:
singularity.enabled = true
singularity.autoMounts = true
process.beforeScript = 'module load apptainer/1.3.6'

process.executor = 'lsf'
process.clusterOptions = '-q premium -R "rusage[mem=64G]" -n 8'  // edit to your queue

params.pascalx_sif = '/sc/arion/projects/<proj>/containers/pascalx_ubuntu20.04.sif'
params.ref_dir     = '/sc/arion/projects/<proj>/ref/pascalx'
Example process:
process PASCALX_GENES {
  container params.pascalx_sif
  cpus 8; memory '64 GB'; time '24h'

  input:
    path gwas_file
  output:
    path "genescore.txt"

  script:
  """
  REF="${params.ref_dir}/EUR.simulated"
  ANN="${params.ref_dir}/ensemble.txt"
  pascalx -p ${task.cpus} ${ANN} ${REF} genescore.txt genescoring \
          -sh True -cr 0 -cp 4 ${gwas_file}
  """
}
This aligns with Minerva’s Apptainer usage and typical Nextflow-on-LSF patterns; you can also borrow the public nf-core Minerva profile if useful. 
Icahn School of Medicine at Mount Sinai
+1
Optional alternatives
Small Linux VM route: spin up an Ubuntu x86_64 VM, install SingularityCE or Apptainer, run singularity|apptainer build ..., then scp the .sif to Minerva. (Apptainer prefers local unprivileged builds; SingularityCE supports --remote.) 
apptainer.org
Quick sanity checklist
 Singularity.pascalx saved on your Mac (above).
 Built .sif via Sylabs Remote Builder and downloaded it. 
cloud.sylabs.io
 Copied .sif to /sc/arion/projects/<proj>/containers/ and verified pascalx -h. 
Icahn School of Medicine at Mount Sinai
 Nextflow points to the .sif; Apptainer module autoloads per job; refs cached under /sc/arion/projects/<proj>/ref/pascalx/.
Follow-up questions for you
Which LSF queue / default resources should I bake into process.clusterOptions (queue, mem, walltime)?
Are you planning CPU-only initially, or should we add --nv (GPU) container options now for future PascalX runs? 
Bergmann Lab
Do you want a tiny one-time “setup” Nextflow process to build/import your preferred 1KG reference into params.ref_dir so all jobs reuse it automatically?