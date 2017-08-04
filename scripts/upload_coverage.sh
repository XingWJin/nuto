#!/bin/bash
set -ev
if [[ "$COVERAGE" == "TRUE" ]]; then
     # get coverage information
    docker exec docker_container lcov --capture --directory /build --output-file coverage.info

    # filter out system stuff that we don't control
    EXTERNALFILES=$(pwd)/external/*
    docker exec docker_container lcov -r coverage.info '/usr/include/*' -o coverage.info
    docker exec docker_container lcov -r coverage.info "$EXTERNALFILES" -o coverage.info

    # upload to codecov
    docker exec docker_container bash <(curl -s https://codecov.io/bash) -f coverage.info
fi
