use strict;
use warnings;
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'clustalUtils',
    license             => 'perl',
    dist_abstract       => 'Utilities for clustal file handling',
    dist_author         => 'Mark Kelley <meb223@cs.unh.edu>',
    build_requires => {
        'Test::More' => '0.10',
    },
);

$builder->create_build_script();