FROM tjhei/dealii:v9.0.1-full-v9.0.1-r5-gcc5

LABEL maintainer <rene.gassmoeller@mailbox.org>

ADD make_logo.sh build_aspect.sh ./

RUN bash ./build_aspect.sh

CMD bash make_logo.sh

