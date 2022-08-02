# Choose any base image you prefer, like bioconda, ubuntu, alpine, etc..
# or using a multistage deployment
FROM centos:7.6.1810

# Use a single LABEL (multiple ones are discouraged because each one creates a new layer)
LABEL	description="Sample file validator" vendor="ISCIII" version="0.1" maintainer="Sara Monzon <smonzon@isciii.es>"

# First, install all the dependencies needed by the results checking procedure
RUN echo "Install basic development tools" && \
    yum -y update && yum -y install wget curl epel-release && \
    yum -y groupinstall "Development Tools" && \
    echo "Install python2.7 setuptools and pip" && \
    yum -y install python-setuptools python-pip python-wheel python3-pip python3-wheel python3-setuptools && \
    pip install --upgrade 'pip<21' wheel

RUN echo "Installing SCI-F" && \
    pip install 'setuptools<45' 'decorator<5' 'traitlets<5' scif==0.0.77 'ipython<6'

COPY ./scif_app_recipes /opt/scif_app_recipes
COPY ./containers/checkFormat/checkTreeFormat.py /opt/

# Now, we get the validation program in some way (copying, git clone, etc...)
# which should also write the results in the canonical format expected by the
# community metrics docker instances.
# The return value tells whether the validation and the translation to
# the canonical format whether worked or failed

RUN	echo "Installing checkTreeFormat script" && \
         scif install /opt/scif_app_recipes/checkTreeFormat_v1.0_centos7.scif

ENV PATH /scif/apps/checktreeformat/bin:$PATH

#RUN find /scif/apps -maxdepth 2 -name "bin" | while read in; do echo "export PATH=\$PATH:$in" >> /etc/profile.d/myenv.sh;done
#RUN find /scif/apps -maxdepth 2 -name "bin" | xargs -I % echo "ln -s %/* /usr/local/bin/" | bash

#ENTRYPOINT	["/opt/docker-entrypoint.sh"]
#CMD [checkTreeFormat.py]

#RUN find /scif/apps -maxdepth 2 -name "bin" | while read in; do echo "export PATH=\$PATH:$in" >> /etc/bashrc;done
#RUN if [ -z "${LD_LIBRARY_PATH-}" ]; then echo "export LD_LIBRARY_PATH=/usr/local/lib" >> /etc/bashrc;fi
#RUN find /scif/apps -maxdepth 2 -name "lib" | while read in; do echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$in" >> /etc/bashrc;done
