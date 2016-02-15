#!/bin/bash

export PERL5LIB=$PWD/bioperl-live-bioperl-release-1-2-3:$PWD/ensembl-test/modules:$PWD/ensembl/modules:$PWD/ensembl-variation/modules:$PWD/ensembl-funcgen/modules:$PWD/Bio-HTS/lib

export PATH=$PATH:$PWD/tabix

export HTSLIB_DIR=$PWD/htslib

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/htslib

echo "Running test suite"
echo "Using $PERL5LIB"
if [ "$COVERALLS" = 'true' ]; then
  PERL5OPT='-MDevel::Cover=+ignore,bioperl,+ignore,ensembl' perl $PWD/ensembl-test/scripts/runtests.pl -verbose $PWD/scripts/variant_effect_predictor/t $SKIP_TESTS
else
  perl $PWD/ensembl-test/scripts/runtests.pl $PWD/scripts/variant_effect_predictor/t $SKIP_TESTS
fi

rt=$?
if [ $rt -eq 0 ]; then
  if [ "$COVERALLS" = 'true' ]; then
    echo "Running Devel::Cover coveralls report"
    cover --nosummary -report coveralls
  fi
  exit $?
else
  exit $rt
fi
