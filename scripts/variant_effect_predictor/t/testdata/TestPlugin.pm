package TestPlugin;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub get_header_info {
  return {
    TestPlugin => 'Test plugin',
  };
}

sub run {
  my ($self, $tva) = @_;
  return { 'TestPlugin' => $tva->transcript->stable_id() };
}

1;
