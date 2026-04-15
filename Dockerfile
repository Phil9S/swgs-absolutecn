# Dockerfile build from pixi env
# Based on https://github.com/pavelzw/pixi-docker-example
FROM ghcr.io/prefix-dev/pixi:0.67.0 AS build

# copy source code, pixi.toml and pixi.lock to the container
WORKDIR /app
COPY install_env.sh .
COPY pixi.* .
ADD resources/ /
# install dependencies to `/app/.pixi/envs/prod`
# use `--locked` to ensure the lockfile is up to date with pixi.toml
RUN pixi install --frozen --locked --run-post-link-scripts
#RUN pixi run ./install_env.sh
#RUN pixi run export TAR="tar --no-same-permissions --no-same-owner" && export TMPDIR=/tmp
RUN export TAR="tar --no-same-permissions --no-same-owner" && export TMPDIR=/tmp && pixi run Rscript -e 'remotes::install_github(repo = "markowetzlab/QDNAseq.hg38",quiet=TRUE,upgrade=FALSE,force=FALSE)'
RUN export TAR="tar --no-same-permissions --no-same-owner" && export TMPDIR=/tmp && pixi run Rscript -e 'remotes::install_github(repo = "markowetzlab/QDNAseqmod",quiet=TRUE,upgrade=FALSE,force=FALSE)'

#RUN pixi shell-hook -s bash > /shell-hook
#RUN echo "#!/bin/bash" > /app/entrypoint.sh
#RUN cat /shell-hook >> /app/entrypoint.sh
# extend the shell-hook script to run the command passed to the container
#RUN echo 'exec "$@"' >> /app/entrypoint.sh

FROM debian:12-slim AS production
RUN TZ=Etc/UTC && \
ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && \
echo $TZ > /etc/timezone
WORKDIR /app
COPY . .
# only copy the production environment into prod container
# please note that the "prefix" (path) needs to stay the same as in the build container
COPY --from=build /app/.pixi/envs/default /app/.pixi/envs/default
#COPY --from=build --chmod=0755 /app/entrypoint.sh /app/entrypoint.sh

ENV PATH=/app/.pixi/envs/default/bin:$PATH
ENV CONDA_PREFIX=/app/.pixi/envs/default

#ENTRYPOINT [ "/app/entrypoint.sh" ]
