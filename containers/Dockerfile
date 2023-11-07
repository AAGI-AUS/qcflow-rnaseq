FROM mambaorg/micromamba:0.25.1

LABEL image.author.name "Kristina K. Gagalova"
LABEL image.author.email "kristina.gagalova@curtin.edu.au"

COPY --chown=$MAMBA_USER:$MAMBA_USER env.yml /tmp/env.yml

RUN micromamba create -n qcflow-rnaseq

RUN micromamba install -y -n qcflow-rnaseq -f /tmp/env.yml && \
    micromamba clean --all --yes

ENV PATH /opt/conda/envs/qcflow-rnaseq/bin:$PATH
