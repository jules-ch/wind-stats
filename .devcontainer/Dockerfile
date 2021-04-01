# See here for image contents: https://github.com/microsoft/vscode-dev-containers/tree/v0.154.2/containers/python-3/.devcontainer/base.Dockerfile

# [Choice] Python version: 3, 3.9, 3.8, 3.7, 3.6
ARG VARIANT="3"
FROM mcr.microsoft.com/vscode/devcontainers/python:0-${VARIANT}

# python
ENV POETRY_VERSION=1.1.5 \
    POETRY_HOME="/usr/local/poetry" \
    POETRY_VIRTUALENVS_IN_PROJECT=true

ENV PATH="$POETRY_HOME/bin:$PATH"

RUN curl -sSL https://raw.githubusercontent.com/sdispater/poetry/master/get-poetry.py | python \
    && chmod 755 $POETRY_HOME/bin/poetry