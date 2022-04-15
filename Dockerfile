FROM ubuntu

RUN apt-get update -y && apt-get install -y  \
      ca-certificates \
      build-essential libboost-all-dev libboost-system-dev libboost-thread-dev libboost-program-options-dev cmake llvm libfreetype6-dev catch2 \
      curl && \
      curl -LO https://storage.googleapis.com/kubernetes-release/release/$(curl -s https://storage.googleapis.com/kubernetes-release/release/stable.txt)/bin/linux/amd64/kubectl  && \
      chmod +x ./kubectl  && \
      mv ./kubectl /usr/local/bin && \
      rm -rf /var/cache/apk/*

ADD target/cross/x86_64-unknown-linux-musl/release/cheminee /usr/bin/cheminee