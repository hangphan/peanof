#################################################################
# Source Image
FROM python:3.7.5-slim


# Set noninterative mode
ENV DEBIAN_FRONTEND noninteractive

################## BEGIN INSTALLATION ######################
LABEL version="1"
LABEL software="peanof"
LABEL software.version="1.0"
LABEL description="Tools flagging outliers in paediatric heigh/weight measurements based on their SDS values (standard z-score)"
LABEL website="https://github.com/hangphan/peanof"
LABEL documentation="https://github.com/hangphan/peanof"
LABEL license="https://github.com/hangphan/peanof"
LABEL tags="Clinical, EHR, EPR"

# Maintainer
MAINTAINER Hang Phan <hangphan@gmail.com>

RUN python -m pip install \
        pandas \
        ggplot \
	statsmodels \
	sklearn \
	matplotlib \
	numpy 

RUN mkdir /data /config
# Add user 'user' with password 'user'

#RUN groupadd fuse && \
#    useradd --create-home --shell /bin/bash --user-group --uid 1000 --groups sudo,fuse user && \
#    echo `echo "xuser\xuser\n" | passwd user` && \
#    chown xuser:xuser /data && \
#    chown xuser:xuser /config

# Change user
RUN apt-get -y update && apt-get -y install git && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /usr/share/man/?? /usr/share/man/??_*

RUN sed -i 's/pd.tslib/pd/g'  /usr/local/lib/python3.7/site-packages/ggplot/utils.py 
RUN sed -i 's/pd.tslib/pd/g'  /usr/local/lib/python3.7/site-packages/ggplot/stats/smoothers.py
RUN sed  -i 's/from pandas.lib import Timestamp/from pandas import Timestamp/g'  /usr/local/lib/python3.7/site-packages/ggplot/stats/smoothers.py


RUN useradd -ms /bin/bash xuser
RUN chown xuser:xuser /data
USER xuser
ENV PEANOF_FOLDER=/home/xuser/peanof
# get peanof scripts
RUN echo 'Downloading peanof 3'
RUN echo "cache-bust" &&\
   git clone https://github.com/hangphan/peanof ${PEANOF_FOLDER}

# ENV path for scripts
ENV PATH ${PEANOF_FOLDER}/:$PATH
RUN chmod +x ${PEANOF_FOLDER}/peanof.py
 

WORKDIR /data/
VOLUME /data/

CMD ["peanof.py", "-h"]