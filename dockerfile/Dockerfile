FROM python:3.10-slim-buster

RUN apt-get update \
    # Install aws-lambda-cpp build dependencies
    && apt-get install -y \
      g++ \
      make \
      cmake \
      unzip \
    # cleanup package lists, they are not used anymore in this image
    && rm -rf /var/lib/apt/lists/* \
    && apt-cache search linux-headers-generic

ARG FUNCTION_DIR="/function"
ARG BINS_DIR="/function/bin"

# Copy function code
RUN mkdir -p ${FUNCTION_DIR}
RUN mkdir -p ${BINS_DIR}

# Update pip
RUN pip install --upgrade --ignore-installed pip wheel six setuptools \
    && pip install --upgrade --no-cache-dir --ignore-installed \
        awslambdaric \
        boto3 \
        redis \
        httplib2 \
        requests \
        numpy \
        scipy \
        pandas \
        pika \
        kafka-python \
        cloudpickle \
        ps-mem \
        tblib

# Set working directory to function root directory
WORKDIR ${FUNCTION_DIR}

# Add Lithops
COPY lithops_lambda.zip ${FUNCTION_DIR}
RUN unzip lithops_lambda.zip \
    && rm lithops_lambda.zip \
    && mkdir handler \
    && touch handler/__init__.py \
    && mv entry_point.py handler/

# install vsearch
RUN wget https://github.com/torognes/vsearch/archive/v2.21.1.tar.gz -O vsearch.tar.gz \
    && tar -xvzf vsearch.tar.gz \
    && cd vsearch-2.21.1 \
    && ./autogen.sh \
    && ./configure CFLAG="-O3" CXXFLAGS="-O3" \
    && make \
    && cp bin/vsearch ${BINS_DIR} \
    && cd ..

# install mmseqs
RUN git clone https://github.com/soedinglab/MMseqs2.git \
    && cd MMseqs2 \
    && mkdir build \
    && cd build \
    && cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. .. \
    && make \
    && make install \
    && mv bin/mmseqs ${BINS_DIR} \
    && cd ../..

ENTRYPOINT [ "/usr/local/bin/python", "-m", "awslambdaric" ]
CMD [ "handler.entry_point.lambda_handler" ]

