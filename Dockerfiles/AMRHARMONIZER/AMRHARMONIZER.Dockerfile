FROM python:3.11.9-alpine3.20

RUN apk add --no-cache bash

RUN pip install pandas==2.2.2

WORKDIR /data

COPY amrharmonization.py /usr/src/app/

ENTRYPOINT ["python", "/usr/src/app/amrharmonization.py"]