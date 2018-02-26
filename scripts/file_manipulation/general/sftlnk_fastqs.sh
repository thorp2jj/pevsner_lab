#!/user/bin/env bash

find /mnt/data/PWS -type f -iname "*.fastq.gz" -exec ln -s {} . \;
