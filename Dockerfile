FROM mambaorg/micromamba:latest
ENV PIP_NO_CACHE_DIR=1
WORKDIR /tmp/env
COPY WY/environment.yml WY/requirements.txt ./
RUN micromamba create -y -f environment.yml \
    && micromamba install -y -n npy -c conda-forge vina \
    && micromamba clean -afy
ENV PATH=/opt/conda/envs/npy/bin:$PATH
COPY WY /app/WY
WORKDIR /app/WY
USER root
RUN useradd -m -u 1000 appuser \
    && chown -R appuser:appuser /app/WY
USER appuser

EXPOSE 8501
HEALTHCHECK --interval=30s --timeout=5s --start-period=90s --retries=3 \
    CMD python -c "import urllib.request; urllib.request.urlopen('http://127.0.0.1:8501/_stcore/health', timeout=4)" || exit 1

CMD ["streamlit", "run", "App.py", "--server.port=8501", "--server.address=0.0.0.0", "--server.headless=true"]
