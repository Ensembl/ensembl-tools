#!/bin/bash

export PERL5LIB=$PWD/bioperl-live-bioperl-release-1-2-3:$PWD/ensembl-test/modules:$PWD/ensembl/modules:$PWD/ensembl-variation/modules:$PWD/ensembl-funcgen/modules

export PATH=$PATH:$PWD/tabix

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
