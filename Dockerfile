# Dockerfile build from pixi env
# Based on https://github.com/pavelzw/pixi-docker-example
FROM ghcr.io/prefix-dev/pixi:latest AS build

# copy source code, pixi.toml and pixi.lock to the container
WORKDIR /app
COPY install_env.sh .
COPY pixi.* .
ADD resources/ /
# install dependencies to `/app/.pixi/envs/prod`
# use `--locked` to ensure the lockfile is up to date with pixi.toml
RUN pixi install --frozen --locked --run-post-link-scripts
RUN export TAR="tar --no-same-permissions --no-same-owner" && \
	export TMPDIR=/tmp && \
	pixi run Rscript -e 'remotes::install_github(repo = "markowetzlab/QDNAseq.hg38",quiet=TRUE,upgrade=FALSE,force=FALSE)'
RUN export TAR="tar --no-same-permissions --no-same-owner" && \
	export TMPDIR=/tmp && \
	pixi run Rscript -e 'remotes::install_github(repo = "markowetzlab/QDNAseqmod",quiet=TRUE,upgrade=FALSE,force=FALSE)'
RUN export TAR="tar --no-same-permissions --no-same-owner" && \
	export TMPDIR=/tmp && \
	pixi run Rscript -e 'remotes::install_github(repo = "markowetzlab/r-swgs-absolutecn",quiet=TRUE,upgrade=FALSE,force=FALSE)'

FROM debian:12-slim AS production
RUN TZ=Etc/UTC && \
ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && \
echo $TZ > /etc/timezone
RUN apt-get update && apt-get install -y locales \
    && locale-gen en_US.UTF-8 \
    && update-locale LANG=en_US.UTF-8

ENV LANG=en_US.UTF-8
ENV LANGUAGE=en_US:en
ENV LC_ALL=en_US.UTF-8

WORKDIR /app
COPY . .
# only copy the production environment into prod container
# please note that the "prefix" (path) needs to stay the same as in the build container
COPY --from=build /app/.pixi/envs/default /app/.pixi/envs/default

ENV PATH=/app/.pixi/envs/default/bin:$PATH
ENV CONDA_PREFIX=/app/.pixi/envs/default
