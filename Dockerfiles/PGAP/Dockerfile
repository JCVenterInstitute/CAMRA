FROM ncbi/pgap:2024-07-18.build7555 AS app

ARG PGAP_VER="2024-07-18.build7555"

LABEL base.image="ncbi/pgap:2024-07-18.build7555"
LABEL dockerfile.version="1"
LABEL software="NCBI Prokaryotic Genome Annotation Pipeline"
LABEL software.version="2024-07-18.build7555"
LABEL description="The NCBI Prokaryotic Genome Annotation Pipeline is designed to annotate bacterial and archaeal genomes (chromosomes and plasmids)."

RUN pip install --no-input cwltool
RUN pip install --no-input pyaml
RUN pip install --no-input --force-reinstall "urllib3==1.26.6"

RUN curl -L -o ${PGAP_VER}.tgz https://github.com/ncbi/pgap/archive/refs/tags/${PGAP_VER}.tar.gz && \
    tar -zxvf ${PGAP_VER}.tgz && \
    rm -v ${PGAP_VER}.tgz && \
    if [ ! -d "/pgap/pgap-${PGAP_VER}" ]; then mv -v pgap-${PGAP_VER} /pgap; fi

WORKDIR /data