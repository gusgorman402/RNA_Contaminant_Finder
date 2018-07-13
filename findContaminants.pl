#!/usr/bin/perl -w

use Bio::Taxon;
use Bio::DB::Taxonomy;
use Bio::Tree::Tree;
use Bio::SearchIO;
use strict;


my $tree_functions = Bio::Tree::Tree->new();
my $db = Bio::DB::Taxonomy->new( -source => 'flatfile',
                                 -directory => '/home/mwachholtz/tmp',
                                 -nodesfile => '/home/mwachholtz/taxonomy/nodes.dmp',
                                 -namesfile => '/home/mwachholtz/taxonomy/names.dmp' );
my $report = Bio::SearchIO->new( -format => 'blastxml', 
                                 -file => $ARGV[0] );
my $numFile = ">".$ARGV[0].".contaminants";
my $taxonFile = ">".$ARGV[0].".taxons";
open( NUMS, $numFile );
open( TAXONS, $taxonFile );



my ( $taxonName, $isPlant, $plantHit, $counter );
my ( $numPlant, $numBug, $numOther, $numUnlisted, $numBacteria, $numFungi, $numArchaea, $numVirus );

my %gi_taxid_table = ();
open( ACC, "/home/mwachholtz/taxonomy/gi_taxid_prot.dmp");

my $count = 0;
print "Creating hash\n";
while( my $fline = <ACC> )
{
    $count++;
    if( $count % 100000 == 0 ){ print $count,"\n" }
    chomp $fline;
    $fline =~ /^(\d+)\s+(\d+)/;
    $gi_taxid_table{$1} = $2;
}
close( ACC );
print "Hash created\n";

while( my $result = $report->next_result )
{
    $counter = 0;
    $isPlant = 0;
    $plantHit = 0;
    $numPlant = 0;
    $numBug = 0;
    $numFungi = 0;
    $numBacteria = 0;
    $numArchaea = 0;
    $numVirus = 0;
    $numOther = 0;
    $numUnlisted = 0;

    print NUMS $result->query_name;

    while( my $hit = $result->next_hit )
    {
        $counter++;
       
        my $gi = $hit->ncbi_gi;
        #print "GI for hit = $gi\n";
        #my $gi_taxid_line = `grep "gi:$gi " ~/taxonomy/gi_taxid_list.txt`;
        $taxonName = "";
        
        if( defined( $gi_taxid_table{$gi} ) )
        {
            my $speciesID = $gi_taxid_table{$gi};
            my $taxid = $db->get_taxon( -taxonid => $speciesID );

            my @lineages = $tree_functions->get_lineage_nodes( $taxid );
            $taxonName = "";

            foreach my $node ( @lineages )
            {
                my $taxonomy = $db->get_taxon( -taxonid => $node->id() );
                $taxonName = $taxonName.$taxonomy->scientific_name.":";
            }

            if( $counter == 1 && ($taxonName !~ /Viridiplantae/) )
            {
                print TAXONS $result->query_name, " ", $taxonName, $taxid->scientific_name, "\n";
            }

        }
        else
        {
            $taxonName = "Not listed";
        }

        if( $taxonName =~ /Viridiplantae/ && $isPlant == 0 )
        {
            $plantHit = $counter;
            $isPlant = 1;
        }
        if( $plantHit == 1 ){ last }

        if( $taxonName =~ /Viridiplantae/ ){ $numPlant++ }
        if( $taxonName =~ /Arthropoda/ ){ $numBug++ }
        if( $taxonName =~ /Fungi/ ){ $numFungi++ }
        if( $taxonName =~ /Bacteria/ ){ $numBacteria++ }
        if( $taxonName =~ /Archaea/ ){ $numArchaea++ }
        if( $taxonName =~ /Viruses/ || $taxonName =~ /Viroids/ ){ $numVirus++ }
        if( $taxonName =~ /Not\s+listed/ ){ $numUnlisted++ }
        if( $taxonName !~ /Viridiplantae/ && $taxonName !~ /Arthropoda/ && $taxonName !~ /Fungi/ && $taxonName !~ /Bacteria/ && $taxonName !~ /Archaea/
            && $taxonName !~ /Viruses/ && $taxonName !~ /Viroids/ && $taxonName !~ /Not\s+listed/ ){ $numOther++ }
        
    }

    print NUMS ",$counter,$plantHit,$numPlant,$numBug,$numFungi,$numBacteria,$numArchaea,$numVirus,$numOther,$numUnlisted\n";
}

