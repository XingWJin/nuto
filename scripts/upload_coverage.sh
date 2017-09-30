#!/bin/bash
set -ev
if [[ "$COVERAGE" == "TRUE" ]]; then
     # get coverage information
    docker exec dock lcov --capture --directory build --output-file coverage.info

    # filter out system stuff that we don't control
    EXTERNALFILES=$(pwd)/external/*
    docker exec dock lcov -r coverage.info '/usr/include/*' -o coverage.info
    docker exec dock lcov -r coverage.info "$EXTERNALFILES" -o coverage.info

    # upload to codecov
    docker exec dock bash <(curl -s https://codecov.io/bash) -f coverage.info
fi
