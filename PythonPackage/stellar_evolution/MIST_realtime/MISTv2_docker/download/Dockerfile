FROM nugrid/nudome:18.0
RUN cd /home/user \
    && wget -q https://zenodo.org/records/15213406/files/mesa11701_MIST2.tgz\
    && echo "MIST MESA downloaded" \
    && tar -xzf mesa11701_MIST2.tgz \
    && echo "MIST MESA unpacked" \
    && mv mesa11701_MIST2 mesa \
    && rm mesa11701_MIST2.tgz \
