FROM alpine

RUN apk add --update \
      ca-certificates \
      bash boost boost-dev cmake llvm flex-dev automake build-base libressl-dev openssl freetype \
      curl && \
      curl -LO https://storage.googleapis.com/kubernetes-release/release/$(curl -s https://storage.googleapis.com/kubernetes-release/release/stable.txt)/bin/linux/amd64/kubectl  && \
      chmod +x ./kubectl  && \
      mv ./kubectl /usr/local/bin && \
      rm -rf /var/cache/apk/*

ADD target/x86_64-unknown-linux-musl/release/cheminee /usr/bin/cheminee