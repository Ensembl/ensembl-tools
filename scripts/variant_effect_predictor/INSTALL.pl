#!/usr/bin/env perl

use Getopt::Long;
use File::Path;
use File::Copy;
use File::Basename;
use Archive::Extract;
use Net::FTP;
use Cwd;

$| = 1;
our $VERSION = 75;
our $have_LWP;
our $use_curl = 0;
have_LWP();

# CONFIGURE
###########

my ($DEST_DIR, $ENS_CVS_ROOT, $API_VERSION, $BIOPERL_URL, $CACHE_URL, $FASTA_URL, $FTP_USER, $help, $UPDATE, $SPECIES, $AUTO, $QUIET, $PREFER_BIN, $CONVERT, $TEST);

GetOptions(
  'DESTDIR|d=s'  => \$DEST_DIR,
  'VERSION|v=i'  => \$API_VERSION,
  'BIOPERL|b=s'  => \$BIOPERL_URL,
  'CACHEURL|u=s' => \$CACHE_URL,
  'CACHEDIR|c=s' => \$CACHE_DIR,
  'FASTAURL|f=s' => \$FASTA_URL,
  'HELP|h'       => \$help,
  'UPDATE|n'     => \$UPDATE,
  'SPECIES|s=s'  => \$SPECIES,
  'AUTO|a=s'     => \$AUTO,
  'QUIET|q'      => \$QUIET,
  'PREFER_BIN|p' => \$PREFER_BIN,
  'CONVERT|t'    => \$CONVERT,
  'TEST'         => \$TEST,
) or die("ERROR: Failed to parse arguments");

if(defined($help)) {
  usage();
  exit(0);
}

my $default_dir_used;

# check if $DEST_DIR is default
if(defined($DEST_DIR)) {
  print "Using non-default installation directory $DEST_DIR - you will probably need to add $DEST_DIR to your PERL5LIB\n";
  $default_dir_used = 0;
}
else {
  $DEST_DIR ||= '.';
  $default_dir_used = 1;
}

my $lib_dir = $DEST_DIR;

$DEST_DIR       .= '/Bio';
$ENS_GIT_ROOT ||= 'https://github.com/Ensembl/';
$BIOPERL_URL  ||= 'http://bioperl.org/DIST/BioPerl-1.6.1.tar.gz';
$API_VERSION  ||= $VERSION;
$CACHE_URL    ||= "ftp://ftp.ensembl.org/pub/release-$API_VERSION/variation/VEP";
$CACHE_DIR    ||= $ENV{HOME} ? $ENV{HOME}.'/.vep' : 'cache';
$FTP_USER     ||= 'anonymous';
$FASTA_URL    ||= "ftp://ftp.ensembl.org/pub/release-$API_VERSION/fasta/";
$PREFER_BIN     = 0 unless defined($PREFER_BIN);

my $dirname = dirname(__FILE__) || '.';

# using PREFER_BIN can save memory when extracting archives
$Archive::Extract::PREFER_BIN = $PREFER_BIN == 0 ? 0 : 1;

$QUIET = 0 unless $UPDATE || $AUTO;

# set up the URLs
my $ensembl_url_tail = '/archive/release/';
my $archive_type = '.zip';

# auto?
if($AUTO) {
  
  # check
  die("ERROR: Failed to parse AUTO string - must contain any of a (API), c (cache), f (FASTA)\n") unless $AUTO =~ /^[acf]+$/i;
  
  # require species
  if($AUTO =~ /[cf]/i) {
    die("ERROR: No species specified") unless $SPECIES;
    $SPECIES = [split /\,/, $SPECIES];
  }
}

our $prev_progress = 0;

print "\nHello! This installer is configured to install v$API_VERSION of the Ensembl API for use by the VEP.\nIt will not affect any existing installations of the Ensembl API that you may have.\n\nIt will also download and install cache files from Ensembl's FTP server.\n\n" unless $QUIET;

# UPDATE
########

if($UPDATE) {
  eval q{ use JSON; };
  die("ERROR: Updating requires JSON Perl module\n$@") if $@;
  
  print "Checking for newer version of the VEP\n";
  
  eval q{
    use HTTP::Tiny;
  };
  die("ERROR: Updating requires HTTP::Tiny Perl module\n$@") if $@;
  my $http = HTTP::Tiny->new();
 
  my $server = 'http://beta.rest.ensembl.org';
  my $ext = '/info/software?';
  my $response = $http->get($server.$ext, {
    headers => { 'Content-type' => 'application/json' }
  });
  
  die "ERROR: Failed to fetch software version number!\n" unless $response->{success};
  
  if(length $response->{content}) {
    my $hash = decode_json($response->{content});  
    die("ERROR: Failed to get software version from JSON response\n") unless defined($hash->{release});
    
    if($hash->{release} > $VERSION) {
      
      print "Ensembl reports there is a newer version of the VEP ($hash->{release}) available - do you want to download? ";
      
      my $ok = <>;
      
      if($ok !~ /^y/i) {
        print "OK, bye!\n";
        exit(0);
      }
      
      my $url = $ENS_GIT_ROOT.'ensembl-tools'.$ensembl_url_tail.$hash->{release}.$archive_type;
      
      my $tmpdir = '.'.$$.'_tmp';
      mkdir($tmpdir);
      
      print "Downloading version $hash->{release}\n";
      download_to_file($url, $tmpdir.'/variant_effect_predictor'.$archive_type);
      
      print "Unpacking\n";
      unpack_arch($tmpdir.'/variant_effect_predictor'.$archive_type, $tmpdir);
      unlink($tmpdir.'/variant_effect_predictor'.$archive_type);
      
      opendir NEWDIR, $tmpdir.'/ensembl-tools-release-'.$hash->{release}.'/scripts/variant_effect_predictor';
      my @new_files = grep {!/^\./} readdir NEWDIR;
      closedir NEWDIR;
      
      foreach my $new_file(@new_files) {
        if(-e $new_file) {
          print "Backing up $new_file to $new_file\.bak\_$VERSION\n";
          move($new_file, "$new_file\.bak\_$VERSION");
          move($tmpdir.'/ensembl-tools-release-'.$hash->{release}.'/scripts/variant_effect_predictor/'.$new_file, $new_file);
        }
        else {
          print "Copying file $new_file\n";
          move($tmpdir.'/ensembl-tools-release-'.$hash->{release}.'/scripts/variant_effect_predictor/'.$new_file, $new_file);
        }
      }
      
      rmtree($tmpdir);
      
      print "\nLooks good! Rerun INSTALL.pl to update your API and/or get the latest cache files\n";
      exit(0);
    }
    else {
      print "Looks like you have the latest version - no need to update!\n\n";
      print "There may still be post-release patches to the API - run INSTALL.pl without --UPDATE/-n to re-install your API\n";
      exit(0);
    }
  }
}

# CHECK EXISTING
################

goto CACHE if $AUTO && $AUTO !~ /a/i;
print "Checking for installed versions of the Ensembl API..." unless $QUIET;

# test if the user has the API installed
my $has_api = {
  'ensembl' => 0,
  'ensembl-variation' => 0,
  'ensembl-functgenomics' => 0,
};

eval q{
  use Bio::EnsEMBL::Registry;
};

my $installed_version;

unless($@) {
  $has_api->{ensembl} = 1;
  
  $installed_version = Bio::EnsEMBL::Registry->software_version;
}

eval q{
  use Bio::EnsEMBL::Variation::Utils::VEP;
};

$has_api->{'ensembl-variation'} = 1 unless $@;

eval q{
  use Bio::EnsEMBL::Funcgen::RegulatoryFeature;
};

$has_api->{'ensembl-functgenomics'} = 1 unless $@;


print "done\n";

my $total = 0;
$total += $_ for values %$has_api;

my $message;

if($total == 3) {
  
  if(defined($installed_version)) {
    if($installed_version == $API_VERSION) {
      $message = "It looks like you already have v$API_VERSION of the API installed.\nYou shouldn't need to install the API";
    }
    
    elsif($installed_version > $API_VERSION) {
      $message = "It looks like this installer is for an older distribution of the API than you already have";
    }
    
    else {
      $message = "It looks like you have an older version (v$installed_version) of the API installed.\nThis installer will install a limited set of the API v$API_VERSION for use by the VEP only";
    }
  }
  
  else {
    $message = "It looks like you have an unidentified version of the API installed.\nThis installer will install a limited set of the API v$API_VERSION for use by the VEP only"
  }
}

elsif($total > 0) {
  $message = "It looks like you already have the following API modules installed:\n\n".(join "\n", grep {$has_api->{$_}} keys %$has_api)."\n\nThe VEP requires the ensembl, ensembl-variation and optionally ensembl-functgenomics modules";
}

if(defined($message)) {
  print $message unless $QUIET;
  
  my $ok;
  
  if($AUTO) {
    $ok = $AUTO =~ /a/i ? 'y' : 'n';
  }
  else {
    print "\n\nSkip to the next step (n) to install cache files\n\nDo you want to continue installing the API (y/n)? ";
    
    $ok = <>;
  }
  
  if($ok !~ /^y/i) {
    print " - skipping API installation\n" unless $QUIET;
    goto CACHE;
  }
}



# SETUP
#######

print "\nSetting up directories\n" unless $QUIET;

# check if install dir exists
if(-e $DEST_DIR) {
  my $ok;
  
  if($AUTO) {
    $ok = 'y';
  }
  else {
    print "Destination directory $DEST_DIR already exists.\nDo you want to overwrite it (if updating VEP this is probably OK) (y/n)? ";
    
    $ok = <>;
  }
  
  if($ok !~ /^y/i) {
    print "Exiting\n";
    exit(0);
  }
  
  else {
    unless($default_dir_used || $AUTO) {
      print "WARNING: You are using a non-default install directory.\nPressing \"y\" again will remove $DEST_DIR and its contents!!!\nAre you really, really sure (y/n)? ";
      $ok = <>;
      
      if($ok !~ /^y/i) {
        print "Exiting\n";
        exit(0);
      }
    }
    
    # try to delete the existing dir
    rmtree($DEST_DIR) or die "ERROR: Could not delete directory $DEST_DIR\n";
  }
}

mkdir($DEST_DIR) or die "ERROR: Could not make directory $DEST_DIR\n";
mkdir($DEST_DIR.'/tmp') or die "ERROR: Could not make directory $DEST_DIR/tmp\n";


# API
#####

print "\nDownloading required files\n" unless $QUIET;

foreach my $module(qw(ensembl ensembl-variation ensembl-funcgen)) {
  my $url = $ENS_GIT_ROOT.$module.$ensembl_url_tail.$API_VERSION.$archive_type;
  
  print " - fetching $module\n" unless $QUIET;
  
  my $target_file = $DEST_DIR.'/tmp/'.$module.$archive_type;
  
  if(!-e $target_file) {
    download_to_file($url, $target_file);
  }
  
  print " - unpacking $target_file\n" unless $QUIET;
  unpack_arch("$DEST_DIR/tmp/$module$archive_type", "$DEST_DIR/tmp/");
  
  print " - moving files\n" unless $QUIET;
  
  if($module eq 'ensembl') {
    move("$DEST_DIR/tmp/$module\-release\-$API_VERSION/modules/Bio/EnsEMBL", "$DEST_DIR/EnsEMBL") or die "ERROR: Could not move directory\n".$!;
  }
  elsif($module eq 'ensembl-variation') {
    move("$DEST_DIR/tmp/$module\-release-$API_VERSION/modules/Bio/EnsEMBL/Variation", "$DEST_DIR/EnsEMBL/Variation") or die "ERROR: Could not move directory\n".$!;
  }
  elsif($module eq 'ensembl-funcgen') {
    move("$DEST_DIR/tmp/$module\-release-$API_VERSION/modules/Bio/EnsEMBL/Funcgen", "$DEST_DIR/EnsEMBL/Funcgen") or die "ERROR: Could not move directory\n".$!;
  }
  
  rmtree("$DEST_DIR/tmp/$module\-release-$API_VERSION") or die "ERROR: Failed to remove directory $DEST_DIR/tmp/$module\-release-$API_VERSION\n";
}



# BIOPERL
#########

# now get BioPerl
print " - fetching BioPerl\n" unless $QUIET;

$bioperl_file = (split /\//, $BIOPERL_URL)[-1];

my $target_file = $DEST_DIR.'/tmp/'.$bioperl_file;

download_to_file($BIOPERL_URL, $target_file);

print " - unpacking $target_file\n" unless $QUIET;
unpack_arch("$DEST_DIR/tmp/$bioperl_file", "$DEST_DIR/tmp/");

print " - moving files\n" unless $QUIET;

$bioperl_file =~ /(bioperl.+?)\.tar\.gz/i;
my $bioperl_dir = $1;

opendir BIO, "$DEST_DIR/tmp/$bioperl_dir/Bio/";
move("$DEST_DIR/tmp/$bioperl_dir/Bio/$_", "$DEST_DIR/$_") for readdir BIO;
closedir BIO;

rmtree("$DEST_DIR/tmp") or die "ERROR: Failed to remove directory $DEST_DIR/tmp\n";



# TEST
######

print "\nTesting VEP script\n" unless $QUIET;
my $test_vep = `perl $dirname/variant_effect_predictor.pl --help 2>&1`;

$test_vep =~ /ENSEMBL VARIANT EFFECT PREDICTOR/ or die "ERROR: Testing VEP script failed with the following error\n$test_vep\n";

print " - OK!\n" unless $QUIET;

# CACHE FILES
#############

CACHE:

print "\nThe VEP can either connect to remote or local databases, or use local cache files.\nUsing local cache files is the fastest and most efficient way to run the VEP\n" unless $QUIET;
print "Cache files will be stored in $CACHE_DIR\n" unless $QUIET;

my $ok;

if($AUTO) {
  $ok = $AUTO =~ /c/i ? 'y' : 'n';
}
else {
  print "Do you want to install any cache files (y/n)? ";
  
  $ok = <>;
}

if($ok !~ /^y/i) {
  print "Skipping cache installation\n" unless $QUIET;
  goto FASTA;
}

# check cache dir exists
if(!(-e $CACHE_DIR)) {
  if(!$AUTO) {
    print "Cache directory $CACHE_DIR does not exists - do you want to create it (y/n)? ";
  
    my $ok = <>;
  
    if($ok !~ /^y/i) {
      print "Exiting\n";
      exit(0);
    }
  }
  
  mkdir($CACHE_DIR) or die "ERROR: Could not create directory $CACHE_DIR\n";
}

mkdir($CACHE_DIR.'/tmp') unless -e $CACHE_DIR.'/tmp';

# get list of species
print "\Getting list of available cache files\n" unless $QUIET;

my $num = 1;
my $species_list;
my @files;

if($CACHE_URL =~ /^ftp/i) {
  $CACHE_URL =~ m/(ftp:\/\/)?(.+?)\/(.+)/;
  my $ftp = Net::FTP->new($2) or die "ERROR: Could not connect to FTP host $2\n$@\n";
  $ftp->login($FTP_USER) or die "ERROR: Could not login as $FTP_USER\n$@\n";
  
  foreach my $sub(split /\//, $3) {
    $ftp->cwd($sub) or die "ERROR: Could not change directory to $sub\n$@\n";
  }
  
  push @files, sort grep {$_ =~ /tar.gz/} $ftp->ls;
}
else {
  eval q{
    use File::Listing qw(parse_dir);
  };
  push @files, sort map {$_->[0]} grep {$_->[0] =~ /tar.gz/} @{parse_dir(get($CACHE_URL))};
  
  if(!@files) {
    opendir DIR, $CACHE_URL;
    @files = grep {$_ =~ /tar.gz/} readdir DIR;
    closedir DIR;
  }
}

# if we don't have a species list, we'll have to guess
if(!scalar(@files)) {
  print "Could not get current species list - using predefined list instead\n";
  
  @files = (
    "bos_taurus_vep_$API_VERSION.tar.gz",
    "danio_rerio_vep_$API_VERSION.tar.gz",
    "homo_sapiens_vep_$API_VERSION.tar.gz",
    "mus_musculus_vep_$API_VERSION.tar.gz",
    "rattus_norvegicus_vep_$API_VERSION.tar.gz",
  );
}

foreach my $file(@files) {
  $species_list .= $num++." : ".$file."\n";
}

my @store_species;
my @indexes;

if($AUTO) {
  if($SPECIES->[0] eq 'all') {
    @indexes = (1..(scalar @files));
  }
  
  else {
    foreach my $sp(@$SPECIES) {
      for my $i(0..$#files) {
        if($sp =~ /refseq/i) {
          push @indexes, $i + 1 if $files[$i] =~ /$sp/i;
        }
        else {
          push @indexes, $i + 1 if $files[$i] =~ /$sp/i && $files[$i] !~ /refseq/i;
        }
      }
    }
  }
  
  die("ERROR: No matching species found") unless scalar @indexes;
  
  # uniquify and sort
  @indexes = sort {$a <=> $b} keys %{{map {$_ => 1} @indexes}};
}
else {
  print "The following species/files are available; which do you want (can specify multiple separated by spaces or 0 for all): \n$species_list\n? ";
  @indexes = split /\s+/, <>;
  
  # user wants all species found
  if(scalar @indexes == 1 && $indexes[0] == 0) {
    @indexes = [1..(scalar @files)];
  }
}

foreach my $file(@indexes) {
  die("ERROR: File number $file not valid\n") unless defined($file) && $file =~ /^[0-9]+$/ && defined($files[$file - 1]);
  
  my $file_path = $files[$file - 1];
  
  my $refseq = 0;
  my ($species, $file_name);
  
  if($file_path =~ /\//) {
    ($species, $file_name) = (split /\//, $file_path);
  }
  else {
    $file_name = $file_path;
    $file_name =~ m/^(\w+?)\_vep/;
    $species = $1;
  }
  
  push @store_species, $species;
  
  # check if user already has this species and version
  if(-e "$CACHE_DIR/$species/$API_VERSION") {
    
    my $ok;
    
    print "\nWARNING: It looks like you already have the cache for $species (v$API_VERSION) installed.\n" unless $QUIET;
    
    if($AUTO) {
      print "\nDelete the folder $CACHE_DIR/$species/$API_VERSION and re-run INSTALL.pl if you want to re-install\n";
    }
    else {
      print "If you continue the existing cache will be overwritten.\nAre you sure you want to continue (y/n)? ";
      
      $ok = <>;
    }
    
    if($ok !~ /^y/i) {
      print " - skipping $species\n" unless $QUIET;
      next;
    }
    
    rmtree("$CACHE_DIR/$species/$API_VERSION") or die "ERROR: Could not delete directory $CACHE_DIR/$species/$API_VERSION\n";
  }
  
  if($species =~ /refseq/i) {
    print "NB: Remember to use --refseq when running the VEP with this cache!\n" unless $QUIET;
  }
  
  my $target_file = "$CACHE_DIR/tmp/$file_name";
  
  if($CACHE_URL =~ /^ftp/) {
    print " - downloading $CACHE_URL/$file_path\n" unless $QUIET;
    download_to_file("$CACHE_URL/$file_path", $target_file) unless $TEST;
  }
  else {
    print " - copying $CACHE_URL/$file_path\n" unless $QUIET;
    copy("$CACHE_URL/$file_path", $target_file) unless $TEST;
  }
  
  print " - unpacking $file_name\n" unless $QUIET;
  
  
  unpack_arch($target_file, $CACHE_DIR.'/tmp/') unless $TEST;
  
  # does species dir exist?
  if(!-e "$CACHE_DIR/$species" && !$TEST) {
    mkdir("$CACHE_DIR/$species") or die "ERROR: Could not create directory $CACHE_DIR/$species\n";
  }
  
  # move files
  unless($TEST) {
    opendir CACHEDIR, "$CACHE_DIR/tmp/$species/";
    move("$CACHE_DIR/tmp/$species/$_", "$CACHE_DIR/$species/$_") for readdir CACHEDIR;
    closedir CACHEDIR;
  }
  
  # convert?
  if($CONVERT && !$TEST) {
    print " - converting cache\n" unless $QUIET;
    system("perl $dirname/convert_cache.pl --dir $CACHE_DIR --species $species --version $API_VERSION") == 0 or print STDERR "WARNING: Failed to run convert script\n";
  }
}



# FASTA FILES
#############

FASTA:

print "\nThe VEP can use FASTA files to retrieve sequence data for HGVS notations and reference sequence checks.\n" unless $QUIET;
print "FASTA files will be stored in $CACHE_DIR\n" unless $QUIET;

if($AUTO) {
  $ok = $AUTO =~ /f/i ? 'y' : 'n';
}
else {
  print "Do you want to install any FASTA files (y/n)? ";
  
  $ok = <>;
}

if($ok !~ /^y/i) {
  print "Skipping FASTA installation - Exiting\n";
  exit(0);
}

my @dirs = ();
my $ftp;

if($FASTA_URL =~ /^ftp/i) {
  $FASTA_URL =~ m/(ftp:\/\/)?(.+?)\/(.+)/;
  $ftp = Net::FTP->new($2) or die "ERROR: Could not connect to FTP host $2\n$@\n";
  $ftp->login($FTP_USER) or die "ERROR: Could not login as $FTP_USER\n$@\n";
  
  foreach my $sub(split /\//, $3) {
    $ftp->cwd($sub) or die "ERROR: Could not change directory to $sub\n$@\n";
  }
  
  push @dirs, sort $ftp->ls;
}
else {
  opendir DIR, $FASTA_URL;
  @dirs = grep {-d $FASTA_URL.'/'.$_ && $_ !~ /^\./} readdir DIR;
  closedir DIR;
}

$species_list = '';
$num = 1;
foreach my $dir(@dirs) {
  $species_list .= $num++." : ".$dir."\n";
}

my @species;
if($AUTO) {
  if($SPECIES->[0] eq 'all') {
    @species = scalar @store_species ? @store_species : @dirs;
  }
  else {
    @species = scalar @store_species ? @store_species : @$SPECIES;
  }
}
else {
  print "FASTA files for the following species are available; which do you want (can specify multiple separated by spaces, \"0\" to install for species specified for cache download): \n$species_list\n? ";
  
  my $input = <>;
  my @nums = split /\s+/, $input;
  
  @species = @store_species if grep {$_ eq '0'} @nums;
  push @species, $dirs[$_ - 1] for grep {$_ > 0} @nums;
}

foreach my $species(@species) {
  
  # remove refseq name
  my $orig_species = $species;
  $species =~ s/_refseq//;
  
  my @files;
  
  if($ftp) {
    $ftp->cwd($species) or die "ERROR: Could not change directory to $species\n$@\n";
    $ftp->cwd('dna') or die "ERROR: Could not change directory to dna\n$@\n";
    @files = grep {!/_(s|r)m\./} $ftp->ls;
  }
  else {
    if(!opendir DIR, "$FASTA_URL/$species/dna") {
      warn "WARNING: Could not read from directory $FASTA_URL/$species/dna\n$@\n";
      next;
    }
    @files = grep {$_ !~ /^\./} readdir DIR;
    closedir DIR;
  }
  
  my ($file) = grep {/primary_assembly.fa.gz$/} @files;
  ($file) = grep {/toplevel.fa.gz$/} @files if !defined($file);
  
  unless(defined($file)) {
    warn "WARNING: No download found for $species\n";
    next;
  }
  
  my $ex = "$CACHE_DIR/$orig_species/$API_VERSION/$file";
  $ex =~ s/\.gz$//;
  if(-e $ex) {
    print "Looks like you already have the FASTA file for $orig_species, skipping\n" unless $QUIET;
    
    if($ftp) {
      $ftp->cwd('../');
      $ftp->cwd('../');
    }
    next;
  }
  
  # create path
  mkdir($CACHE_DIR) unless -d $CACHE_DIR || $TEST;
  mkdir("$CACHE_DIR/$orig_species") unless -d "$CACHE_DIR/$orig_species" || $TEST;
  mkdir("$CACHE_DIR/$orig_species/$API_VERSION") unless -d "$CACHE_DIR/$orig_species/$API_VERSION" || $TEST;
  
  if($ftp) {
    print " - downloading $file\n" unless $QUIET;
    download_to_file("$FASTA_URL/$species/dna/$file", "$CACHE_DIR/$orig_species/$API_VERSION/$file") unless $TEST;
  }
  else {
    print " - copying $file\n" unless $QUIET;
    copy("$FASTA_URL/$species/dna/$file", "$CACHE_DIR/$orig_species/$API_VERSION/$file") unless $TEST;
  }
  
  print " - extracting data\n" unless $QUIET;
  unpack_arch("$CACHE_DIR/$orig_species/$API_VERSION/$file", "$CACHE_DIR/$orig_species/$API_VERSION/") unless $TEST;
  
  print " - attempting to index\n" unless $QUIET;
  eval q{
    use Bio::DB::Fasta;
  };
  if($@) {
    print "Indexing failed - VEP will attempt to index the file the first time you use it\n" unless $QUIET;
  }
  else {
    Bio::DB::Fasta->new($ex) unless $TEST;
    print " - indexing OK\n" unless $QUIET;
  }
  
  print "The FASTA file should be automatically detected by the VEP when using --cache or --offline. If it is not, use \"--fasta $ex\"\n\n" unless $QUIET;
  
  if($ftp) {
    $ftp->cwd('../');
    $ftp->cwd('../');
  }
}


# CLEANUP
#########
if(-d "$CACHE_DIR/tmp" && !$TEST) {
  rmtree("$CACHE_DIR/tmp") or die "ERROR: Could not delete directory $CACHE_DIR/tmp\n";
}

print "\nSuccess\n" unless $QUIET;


# SUBS
######

sub download_to_file {
  my ($url, $file) = @_;
  
  $url =~ s/([a-z])\//$1\:21\// if $url =~ /ftp/;
  
  if($use_curl) {
    my $output = `curl $url > $file`;
  }
  
  if(have_LWP()) {
    my $response = getstore($url, $file);
    
    unless($response == 200) {
      # try no proxy
      $ua->no_proxy('github.com');
      
      $response = getstore($url, $file);
      
      unless($response == 200) {
        warn "WARNING: Failed to fetch from $url\nError code: $response\n" unless $QUIET;
        print "Trying to fetch using curl\n" unless $QUIET;
        $use_curl = 1;
        download_to_file($url, $file);
      }
    }
  }
  else {
    my $response = HTTP::Tiny->new(no_proxy => 'github.com')->get($url);
    
    if($response->{success}) {
      open OUT, ">$file" or die "Could not write to file $file\n";
      binmode OUT;
      print OUT $response->{content};
      close OUT;
    }
    else {
      warn "WARNING: Failed to fetch from $url\nError code: $response->{reason}\nError content:\n$response->{content}\nTrying without no_proxy\n" unless $QUIET;
      $response = HTTP::Tiny->new->get($url);
      
      if($response->{success}) {
        open OUT, ">$file" or die "Could not write to file $file\n";
        binmode OUT;
        print OUT $response->{content};
        close OUT;
      }
      else {
        warn "WARNING: Failed to fetch from $url\nError code: $response->{reason}\nError content:\n$response->{content}\n" unless $QUIET;
        print "Trying to fetch using curl\n" unless $QUIET;
        $use_curl = 1;
        download_to_file($url, $file);
      } 
    }
  }
}

sub have_LWP {
  return $have_LWP if defined($have_LWP);
  
  eval q{
    use LWP::Simple qw($ua getstore get);
  };
  
  if($@) {
    $have_LWP = 0;
    warn("Using HTTP::Tiny - this may fail when downloading large files; install LWP::Simple to avoid this issue\n");
    
    eval q{
      use HTTP::Tiny;
    };
    
    if($@) {
      die("ERROR: No suitable package installed - this installer requires either HTTP::Tiny or LWP::Simple\n");
    }
  }
  else {
    $have_LWP = 1;
    
    # set up a user agent's proxy (excluding github)
    $ua->env_proxy;
    
    # enable progress
    eval q{
      $ua->show_progress(1);
    } unless $QUIET;
  }
}

# unpack a tarball
sub unpack_arch {
  my ($arch_file, $dir) = @_;
  
  my $ar = Archive::Extract->new(archive => $arch_file);
  my $ok = $ar->extract(to => $dir) or die $ae->error;
  unlink($arch_file);
}

# update or initiate progress bar
sub progress {
  my ($i, $total) = @_;
  
  my $width = 60;
  my $percent = int(($i/$total) * 100);
  my $numblobs = (($i/$total) * $width) - 2;

  return unless $numblobs != $prev_progress;
  $prev_progress = $numblobs;
  
  printf("\r% -${width}s% 1s% 10s", '['.('=' x $numblobs).($numblobs == $width - 2 ? '=' : '>'), ']', "[ " . $percent . "% ]");
}

sub usage {
    my $usage =<<END;
#---------------#
# VEP INSTALLER #
#---------------#

version $VERSION

By Will McLaren (wm2\@ebi.ac.uk)

http://www.ensembl.org/info/docs/variation/vep/vep_script.html#installer

Usage:
perl INSTALL.pl [arguments]

Options
=======

-h | --help        Display this message and quit

-d | --DESTDIR     Set destination directory for API install (default = './')
-v | --VERSION     Set API version to install (default = $VERSION)
-c | --CACHEDIR    Set destination directory for cache files (default = '$HOME/.vep/')

-n | --UPDATE      EXPERIMENTAL! Check for and download new VEP versions

-a | --AUTO        Run installer without user prompts. Use a (API), c (cache),
                   f (FASTA) to specify parts to install e.g. -a ac for API and
                   cache
-s | --SPECIES     Comma-separated list of species to install when using --AUTO
-q | --QUIET       Don't write any status output when using --AUTO
-p | --PREFER_BIN  Use this if the installer fails with out of memory errors

-t | --CONVERT     Convert downloaded caches to use tabix for retrieving
                   co-located variants (requires tabix)
                   
                   
-u | --CACHEURL    Override default cache URL; this may be a local directory or
                   a remote (e.g. FTP) address.
-f | --FASTAURL    Override default FASTA URL; this may be a local directory or
                   a remote (e.g. FTP) address. The FASTA URL/directory must have
                   gzipped FASTA files under the following structure:
                   [species]/[dna]/
END

    print $usage;
}
