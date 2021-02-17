#!/gsc/bin/perl

#  modDFS: Find whether a KEGG module is complete in a proteome
#  
#  Author: Rahul Tyagi
#  Date: 3/12/2013
#  
#  ===========================
#  
#  PROGRAM DESCRIPTION
#  -------------------
#  
#  Given module definitions and reaction information based on a KEGG database, modDFS does a depth first search to find whether a set of KOs are sufficient for completion of a module.
#  For cases where only 1 absent KO would've been enough for module completion, that module is reported to be complete(only lenient).
#  For cyclic modules, completion is defined as having all the module reaction steps present. Again, lenient completion cases are reported.
#  Also, all absent reactions are reported for an cyclic module not strictly complete, so that a manual inspection of completion can be done for modules that are partially cyclic.
#  
#  The user needs to provide following information:
#  1. Available KOs in the dataset, module's definition, and reaction data in input files (see below).
#  2. whether the module is (even partially) cyclic (default is assumed to be acyclic, and the results for such modules need to be manually checked).
#  3. Any "given compounds" in the module. i.e. compounds that are available to the module from external sources. these are generally the substrates of the first reaction step of an acyclic module (this is assumed to be the default), but do not necessarily exclude other compounds. For instance, there might be other reactions whose substrates aren't produced by any reaction in the module. these subsrates need to be pointed out by the user.
#  
#  There are some modules that have "forks" in their definition (e.g. module M00074). i.e. there is more than one product of the module, each of them produced by a different branch of the module. In these cases, it is recommended that the user manually generates module definition files that divide this module into smaller modules with one branch each.
#  
#  
#  ===========================
#  

use warnings;
use Getopt::Long;
use Storable 'dclone';

my ($keggdb,$modfile);
my ($module,$kofile,$rnfile);
my $help;
my ($forcerxn,$last_comp,$first_comp,$ends,$recursive,$cyc)=("","","",0,0,0);
my ($givfile,$given);

GetOptions(
	   "module=s"         => \$module,      #required
       "kofile=s"         => \$kofile,      #required
	   "rnfile=s"         => \$rnfile,      #required
	   "keggdb:s"         => \$keggdb,      #alternate required with modfile
	   "modfile:s"        => \$modfile,     #alternate required with keggdb
	   "givfile:s"        => \$givfile,     #alternate required with given
	   "given:s"          => \$given,       #alternate required with givfile
	   "forcerxn:s"       => \$forcerxn,    #optional
	   "ends"             => \$ends,        #optional 
	   "last_comp:s"      => \$last_comp,   #optional
	   "first_comp:s"     => \$first_comp,  #optional
	   "recursive"        => \$recursive,   #optional
	   "cyc"              => \$cyc,         #optional
	   "help"             => \$help,        # optional
	   );

# ----------------------------------------------------------------------
# C H E C K    I N P U T
# ----------------------------------------------------------------------
my $usage = "

Usage: modDFS.pl [Required Options] [Required alternate options] [optional options]

Order of the options doesn't matter. Abbreviation is allowed, if unambigous.

Please read the accompanying README file for important explanations about some of the options, input and output.


Required options:

 -module   (required) Name of the KEGG module

 -kofile   (required) Path to the species KO file. formatted as in example/example.ko

 -rnfile   (required) Path to the reactions data file. contains reversibility and KO mapping of reactions with one reaction information per line. reversibility is indicated by a 0 for irreversible,1 for reversible, and NA for not available. KOs are in the third column, separated by colon


Required alternate options (one of the two is required):

 -keggdb   (required if -modfile isn't provided) Path to the KEGG database directory
 -modfile  (required if -keggdb isn't provided)  File that contains information about given compounds and whether the module is cyclic

 -givfile  (required if -given isn't provided) File which has a 'c' in second column if the module is (even partially) cyclic and the 'given compounds' of the module in the third column, separated by colon (:). the columns are separated by whitespace. this information can also directly be provided on command line using -cyc and -given
 -given    (required if -givfile isn't provided)  Given compounds separated by colon (:)


Optional options:

 -forcerxn (optional) Rxn id that is to be forced to be present. this reaction will be considered present even if its corresponding KOs aren't in the kofile

 -ends     (optional) Signifies that the bounding substrates (last_comp and first_comp) are provided at command line. this can be used to check completion of a smaller part of the whole module without having to manipulating the module definition file.
 -last_comp   (optional) Command line enforced last_comp that needs to be obtained from first_comp. see -ends
 -first_comp  (optional) Command line enforced first_comp from which last_comp needs to be obtained. see -ends

 -recursive   (optional) Option internal to the program, signifies that the progrma is being called from an instance of itself rather than the user

 -cyc      (optional)    Option indicating that the module has to be treated as cyclic. i.e. complete only if all steps are complete and output contains all incomplete steps. This should be provided even if only a part of the module is cyclic.

 -help     (optional)  View this help information.

";
my $usage1 = "

ERROR!

first_comp and last_comp should be provided if and only if -ends is set

";

die $usage unless ((defined $module) and (defined $kofile) and (defined $rnfile) and ((defined $keggdb) or (defined $modfile)) and ((defined $givfile)or (defined $given)));
die $usage1 unless (($ends==0 && $last_comp eq "" && $first_comp eq "")||($ends==1 && $last_comp ne "" && $first_comp ne ""));


my $warning1 ="

WARNING: both -cyc and -givfile are provided. -cyc will prevail over givfile information.

";
my $warning2 ="

WARNING: both -given and -givfile are provided. -given will prevail over givfile information.

";
my $warning3 ="

WARNING: both -keggdb and -modfile are provided. -keggdb will prevail over modfile information.

";
my $warning4 ="

WARNING: neither -cyc nor -givfile are provided. treating module $module as acyclic.

";

if($recursive!=1){ #warnings are not given when the program calls itself. only resulting value is printed.
	if($cyc==1 && defined $givfile){print $warning1;}
	if($cyc==0 && !(defined $givfile)){print $warning4;}
	if(defined $given && defined $givfile){print $warning2;}
	if(defined $keggdb && defined $modfile){print $warning3;}
}

if(defined $keggdb){$modfile="$keggdb/module/module";}


my $usage2 = "

ERROR!

File '$rnfile' isn't in correct format.  see file example/rn_data. 

reversibility is indicated by a 0 for irreversible,1 for reversible, and NA for not available.
KOs are in the third column, separated by colon

";
my $usage3 = "

ERROR!

File '$kofile' isn't in correct format.  see file example/example.ko 

";

my $givhash; #hash to store "given" compounds of the module. i.e. compounds that are assumed to be present regardless of any reaction presence.
if(defined $givfile){
	open(GIV,$givfile) or die("couldn't open the given data file $givfile\n\n");
	while(<GIV>){
		chomp;
		my($mod,$yc,$given1)=split(/\s+/);
		if($mod eq $module){
			if($cyc!=1){$cyc=($yc eq 'c')?1:0;}
			if(defined $given){$given1=$given;}
			if($given1!~/NA/){
				foreach my $giv(split(/:/,$given1)){
					$givhash->{$giv}=1;
				}
			}
		}
	}
	close GIV;
}else{
	foreach my $giv(split(/:/,$given)){
		$givhash->{$giv}=1;
	}
}

my $ko_hash; #hash to save the information about ko presence . 

open(KO,$kofile) or die("ERROR! couldn't open the ko file $kofile\n\n");
while(<KO>){
	if(/^K\d+/){
		my($ko)=(/^(K\d+)/);
		$ko_hash->{$ko}=1;
	}elsif(/\S/){die $usage3;}
}
close KO;


open(RXN, $rnfile) or die("ERROR! couldn't open the reaction data file $rnfile\n");

my $rn_data; #hash to store information about reactions. includes reversibility, KOs, number of KOs, spontaneity, presence in species, and the step number in the module.

while(<RXN>){
	chomp;
	my @data=split(/\s+/);
	if(scalar @data !=3 || $data[1]!~/(1|0|(NA))/){die $usage2;}
	if($data[1] eq 'NA'){
		$rn_data->{$data[0]}->{'rev'}=0; #if unknown, consider it irreversible
	}else{
		$rn_data->{$data[0]}->{'rev'}=$data[1];
	}
	$rn_data->{$data[0]}->{'present'}=0;
	my @ko=split(':',$data[2]);
	if($ko[0] ne "NA"){
		foreach my $ko(@ko){
			$rn_data->{$data[0]}->{'ko'}->{$ko}=1; #rn_data mapping to ko
			if(exists $ko_hash->{$ko}){
				$rn_data->{$data[0]}->{'present'}=1; #ko is present, so rxn is present
			}
		}
		$rn_data->{$data[0]}->{'numko'}=scalar(@ko);
	}else{
		$rn_data->{$data[0]}->{'numko'}=0;
	}
}
close RXN;

if($forcerxn ne ""){$rn_data->{$forcerxn}->{'present'}=1;} #only one reaction can be "forced". i.e. present even if corresponding KOs might not be present in the ko file

my @spont; #to keep track of the rare 'spontaneous' reactions
my $compsprod;#compsprod is a hash that stores information about the substrates and reactions that lead to the index compound. first index is the product. second is the string with all subsrates (C1+C2 etc.) and the value is the reaction group needed for this reaction (i.e. R1 or R1+R2 etc.)
my @comp;#array that stores all compounds of the module definition
my $comprxn; #hash to store the information about the reactions that the compound is a part of. 1st index is the compound. second is either 'sub' or 'prod' depending on whether the reaction has this compound as substrate or product. 3rd is either 'rxns' (storing the reaction ids) or 'rxnnum' (storing the number of such reactions). 4th (only for 'rxns', not 'rxnnum') is the reaction id. the value is the corresponding prodgrp (for 'sub'->'rxns') or subgrp (for 'prod'->'rxns')
my $rxncomp; #hash  to store substrates of every reaction. first index is reaction (or reaction group R1+R2 etc.) and second is the whole substrate group for that reaction(C1+C2 etc.)
my $comps; #hash indexed by all compounds in the network and storing the compounds that lead to it via a reaction that is present .
my $depends; #hash to save the number and identities of dependencies (i.e. compounds that depend on the indexed compound) 
my @tstepr; #array to store the lines of reaction

my $flag=0; #flag used in module data parsing, keeps track of whether the REACTION of COMPOUND block is being parsed
my $found=0; #flag for the module of interest. it is 1 only when that module record is found

open(IN, $modfile) or die("ERROR! couldn't open the module file $modfile\n\n");
while(<IN>){
	if(/^ENTRY\s+$module\s+Pathway\s+Module/){$found=1;}
	if(/^\/\/\// && $found==1){last;}
	
	if($found==1){
		chomp;
		my $line=$_;
		if(/^REACTION/ && $flag==0){
			$flag=1;
			$line=~s/^REACTION//;

			#ascertaining first_comp from the last substrate in the first reaction of the module definition. if there are more than one reactions, separated by comma, use the last reaction's last substrate.
			my $temp=$line;
			if($first_comp eq ""){
				my ($subs)=($temp=~/^\s+R\S+\s+(.*)\s+<?->\s*.*\s*$/);
				$subs=~s/\s//g;
		
				my @subgrp=split(/\,/,$subs); #some substrates have a comma between them. this means that these are alternative substrates either for same reaction or for different ones, if there are more than one reactions separated by comma on the same step. 
				my $temp1=$subgrp[@subgrp-1];
				my @temp=split(/\+/,$temp1);
				$first_comp=$temp[@temp-1];
			}
		}
		if($flag==1){ #i.e. if inside REACTION block
			if($line=~/^\s/){
				my ($rxns,$subs,$prods)=($line=~/^\s+(R\S+)\s+(.*)\s+<?->\s*(.*)\s*$/);
				$rxns=~s/\s//g;
				push(@tstepr,$rxns); #the whole reaction step, including comma separated and + separated reactions.
				$subs=~s/\s//g;
				$prods=~s/\s//g;
	
				my @rxngrp=split(/\,/,$rxns);	#comma separated reactions form different, alternative reaction groups that can complete this step in module
				for(my $j=0;$j<@rxngrp;$j++){
					my @rxns=split(/\+/,$rxngrp[$j]); #reactions separated by a + are both needed to be complete for the step to be complete. i.e. they are "reaction pairs" that have to happen sequentially in this step. 
					for(my $k=0;$k<@rxns;$k++){
						die ("\nERROR!\ninformation for reaction $rxns[$k] missing from the file $rnfile\n") unless (exists $rn_data->{$rxns[$k]});
						if($prods=~/spontaneous/){
							$rn_data->{$rxns[$k]}->{'spont'}=1;
							$rn_data->{$rxns[$k]}->{'rev'}=0; #a spontaneous reaction isn't considered reversible. 
						}else{$rn_data->{$rxns[$k]}->{'spont'}=0;}
					}
				}
				
				my @subgrp=split(/\,/,$subs); #some substrates have a comma between them. this means that these are alternative substrates either for same reaction or for different ones, if there are more than one reactions separated by comma on the same step. 
				my @prods=split(/\+/,$prods); 

				#for loop to use fillcomps function to fill the hashes defined above. i.e. rxncomp, comprxn, comps, depends etc. 
				for(my $j=0;$j<@prods;$j++){
					$prods[$j]=~s/^\d+//; #stripping the stoichiometric information from products, if any exists. 
					if($j==@prods-1 && $ends==0){$last_comp=$prods[$j]; } #if last_comp isn't provided by user, then last compound of last reaction in the definition is considered to be last_comp. 
					$prods[$j]=~s/\(spontaneous\)//;
					for(my $k=0;$k<@subgrp;$k++){
						
						#if same number of substrate and reaction groups (separated by commas) are there, then nth reaction group corresponds to nth substrate group
						if(@subgrp==@rxngrp){
							$compsprod->{$prods[$j]}->{$subgrp[$k]}=$rxngrp[$k];
							fillcomps($subgrp[$k],$rxngrp[$k],$prods[$j],$prods);
						}else{
							$compsprod->{$prods[$j]}->{$subgrp[0]}=$rxns;#we are assuming that @subgrp>@rxngrp is unlikely (i.e. same reaction with different, alternative substrates. Also, we assume that @subgrp<@rxngrp essentially means R1,R2 C1+C2+C3->P1+P2 etc.. i.e. different reactions mediating the same substrate/product group)
							if($k>0){die("\n\nERROR!! modDFS doesn't handle reaction steps with more than 1 comma-separated substrates involving more than 1 comma-separetd reactions, unless they are equal in number. $subgrp[$k] is part of one such step\n\n");} 
							else{
								for(my $m=0;$m<@rxngrp;$m++){ #if there is just one substrate group and many reaction groups, then all of them have the same substrates.
									fillcomps($subgrp[0],$rxngrp[$m],$prods[$j],$prods);
								}
							}
						}
					}
				}

				#following if-else is to take care of filling compsprod hash for reversible reactions. i.e. treating substrates as products and vice versa.
				if(@subgrp==@rxngrp){
					for(my $j=0;$j<@rxngrp;$j++){
						if($rn_data->{$rxngrp[$j]}->{'rev'} eq '1'){
							my @subs=split(/\+/,$subgrp[$j]);
							for(my $k=0;$k<@subs;$k++){
								$subs[$k]=~s/^\d+//;
								my $prodclean=$prods;
								$prodclean=~s/\(spontaneous\)//;#if the reaction is spontaneous, it shouldn't be reversible, unless something is annotated wrongly. still, we keep this here.
								$compsprod->{$subs[$k]}->{$prodclean}=$rxngrp[$j];
							}
						}
					}
				}else{
					for(my $j=0;$j<@rxngrp;$j++){
						if($rn_data->{$rxngrp[$j]}->{'rev'} eq '1'){
							my @subs=split(/\+/,$subgrp[0]);
							$subs[0]=~s/^\d+//;
							my $prodclean=$prods;
							$prodclean=~s/\(spontaneous\)//;
							$compsprod->{$subs[0]}->{$prodclean}=$rxngrp[$j];
						}
					}
				}
			}elsif($line=~/^COMPOUND/){
				$flag=2; #compound 	block is being parsed
				$line=~s/^COMPOUND//;
			}else{$flag=0;} #for some module definitions where no COMPOUND is defined. though this should anyway break the script. 
		}
		if($flag==2){
			if($line=~/^\s+\S\d\d/){
				my ($comp,$trash)=($line=~/^\s+([A-Z]\d+)/);
				if(!exists $compsprod->{$comp}){$comps->{$comp}->{'present'}=1;} #this means that $comp is a compound that only exists as a substrate in the reaction, and has to be considered 'given' in this module.
				if(exists $givhash->{$comp}){$comps->{$comp}->{'present'}=1;}
				push(@comp,$comp);
			}else{$flag=0;}
		}
	}
}
close IN;

if($found==0){
	print "\n\nno pathway module record for $module was found in the file $modfile!\n\n";
	exit;
}

if($cyc==1){
	my @incsteps; #incomplete steps
	my @incrxn; #incomplete rxn
	foreach my $rxns(@tstepr){
		my $flag1=0;
		my @incrxn1;
		my @rxngrp=split(/\,/,$rxns);	
		for(my $j=0;$j<@rxngrp;$j++){
			my $flag2=1;
			my @incrxn2;
			my @rxns=split(/\+/,$rxngrp[$j]);
			for(my $k=0;$k<@rxns;$k++){
				if(!exists $rn_data->{$rxns[$k]}->{'present'} || $rn_data->{$rxns[$k]}->{'present'}!=1){
					$flag2=0; 
					push(@incrxn2,$rxns[$k]);
				}
			}
			$flag1+=$flag2;
			if(scalar @incrxn2 >1){
				@incrxn2=();
			}elsif(scalar @incrxn2 ==1){push @incrxn1,$incrxn2[0];}
		}
		if($flag1==0){
			push(@incsteps,$rxns);
			push(@incrxn,@incrxn1);
		}else{@incrxn1=();}
	}
	if(scalar @incsteps ==0){
		print "OUTPUT:\n\n$module(cyclic)\tcomplete(strict)\n";
	}elsif(scalar @incsteps>1){
		print "OUTPUT:\n\n$module(cyclic)\tincomplete\n\nincomplete steps:\t";
		print join(";",@incsteps);
		print "\n";
	}else{
		print "OUTPUT:\n\n$module(cyclic)\tcomplete(only lenient)\n\n";
		print "incomplete steps:\t$incsteps[0]\nneeded for completion:\t";
		foreach my $rxn(@incrxn){
			print "\t$rxn;";
			print join(",",(keys %{$rn_data->{$rxn}->{'ko'}}));
		}
		print "\n";
	}
	exit;
}


my $metasp; #metabolites present in the species
foreach my $comp(@comp){
	if(exists $comps->{$comp}){$metasp->{$comp}=1;}
}

#the next do-until loop essentially makes sure that any substrates of a reaction that is forced to go in only one direction (because at least one of the substrates of that reversible reaction cannot be produced by any other reaction in the network and is not a 'given' compound) are not considered a substrate for that reaction in the pathway search. The do-until loop wrapping the foreach loop makes sure that this is repeated for all compounds until the impact of changing one reaction to irreversible percolates throughout the network.
my $del;
my $comprxn1=dclone $comprxn;#temp copy of comprxn
do{
	$del=0;
	foreach my $comp(keys %{$metasp}){
		unless(exists $comps->{$comp}->{'present'}){ #the compound is 'given'
			my $total=0;
			my $onlyrxn;
			foreach my $rxn1(keys %{$comprxn1->{$comp}->{'prod'}->{'rxns'}}){#no. of reactions that have $comp as their product
				if($rn_data->{$rxn1}->{'present'}==1){
					$total++;
					$onlyrxn=$rxn1;
				}
			}
			if($total==1 && $comprxn->{$comp}->{'prod'}->{'rxnnum'}>1){
				foreach my $subgrp1(keys %{$rxncomp->{$onlyrxn}}){
					if($subgrp1=~/$comp/){
						foreach my $subs1(split(/\+/,$subgrp1)){#co-products of $comp in $onlyrxn
							if(exists $comprxn->{$subs1}->{'prod'}->{'rxns'}->{$onlyrxn}){
								foreach my $subs2(split(/\+/,$comprxn->{$subs1}->{'prod'}->{'rxns'}->{$onlyrxn})){#substrates that take part in $onlyrxn to result in $subs1 and its co substrates
									if(exists $comps->{$subs2}->{'altnum'} && $comps->{$subs2}->{'altnum'}>0){
										foreach my $alt(keys %{$comps->{$subs2}}){
											if(exists $comps->{$subs2}->{$alt}->{$subs1} && $comps->{$subs2}->{$alt}->{$subs1}=~/$onlyrxn/){
												delete $comps->{$subs2}->{$alt};
												$comps->{$subs2}->{'altnum'}--;
												if(exists $comprxn1->{$subs2}->{'prod'}->{'rxns'}->{$onlyrxn}){
													delete $comprxn1->{$subs2}->{'prod'}->{'rxns'}->{$onlyrxn};
													$comprxn1->{$subs2}->{'prod'}->{'rxnnum'}--;
												}
												$del++;
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}until($del==0);

do {
	$del=0;
	foreach my $comp(keys %{$metasp}){
		my $flag=1; #flag to remember if the entry needs to be deleted. 1 for keep, 0 for delete
		if(exists $comps->{$comp}->{'present'}){$flag=1;} #$comps is 'given' in the module
		elsif($comps->{$comp}->{'altnum'}==0){$flag=0;}#no alternative subgrps are remaining for this species 
		else{
			$flag=0;
			foreach my $alt(keys %{$comps->{$comp}}){ #to see if all the substrates of all the alternatives that lead to $comp are deleted from this species
				if($alt ne 'altnum'){
					my $flag1=1;
					foreach my $sub(keys %{$comps->{$comp}->{$alt}}){
						if(!exists $metasp->{$sub}){$flag1=0;}
					}
					$flag+=$flag1;
				}
			}
			if($flag>0){$flag=1;}
		}
		if($flag==0){
			$del++;
			delete $metasp->{$comp};
			foreach my $depen(keys %{$depends->{$comp}->{'deps'}}){
				my @alts=split(':',$depends->{$comp}->{'deps'}->{$depen});
				foreach my $alt(@alts){
					if(exists $comps->{$depen}->{$alt}){
						delete $comps->{$depen}->{$alt};
						$comps->{$depen}->{'altnum'}--;
					}
				}
			}
		}
	}
}until($del==0);

#open(PRE,">>modules_new4/$ARGV[0].present");
#open(ABS,">>modules_new4/$ARGV[0].absent");

my($value,@connected)=connected($last_comp,$first_comp);
if($recursive==1){print $value;}
elsif($value==1){
	print "OUTPUT:\n\n$module\tcomplete(strict)\n";
}else{
	my @cmprxn=lenient(); #reactions that'd complete this module
	if(scalar @cmprxn == 0){
		print "OUTPUT:\n\n$module\tincomplete\n";
	}else{
		print "OUTPUT:\n\n$module\tcomplete(only lenient)\n\nneeded for completion:";
		foreach my $rxn(@cmprxn){
			print "\t$rxn;";
			print join(",",(keys %{$rn_data->{$rxn}->{'ko'}}));
		}
		print "\n";
	}
}

exit;

sub connected{
	my($Cp, $Cs)=@_;
	my @connected;
	if((!exists $metasp->{$Cp})||(!exists $metasp->{$Cs})){return (0,@connected);}
	elsif($Cp eq $Cs){return (1,@connected);}
	delete $metasp->{$Cp};
	my ($remain,@check,$conn); #hash that stores remainig;stack that stores to be checked and hash to store connected comps
	my $tot=0;
	foreach my $comp(keys %{$metasp}){
		$remain->{$comp}=0;
		$tot++;
	}
	push(@check,$Cp);
#	delete $remain->{$Cp};

	while(@check>0){
		my $cur=pop @check;
		push(@connected,$cur);
		if($cur eq $Cs){return (1,@connected);}
		$conn->{$cur}=0;
		if($comps->{$cur}->{'altnum'}>0){
			foreach my $alt(keys %{$comps->{$cur}}){
				foreach my $sub(keys %{$comps->{$cur}->{$alt}}){
					#next foreach loop is to count the number of different reactions (i.e. rxngrps) that are part of different alternates of this substrate. This will be used later to make sure that any substrate that has only one reaction left after knocking out any absent reactions in that species, is only considered as a product in any pathway, and not as a substrate, since it can't be present in this species in this module as all other enzymes are absent.
#					my $altrxns;
#					$altrxns->{'num'}=0;
#					foreach my $alt1(keys %{$comps->{$sub}->{$sp}}){
#						foreach my $sub1(keys %{$comps->{$sub}->{$sp}->{$alt1}}){
#							my ($trash,$rxnleft)=split(':',$comps->{$sub}->{$sp}->{$alt1}->{$sub1});
#							unless(exists $altrxns->{$rxnleft}){
#								$altrxns->{$rxnleft}=1;
#								$altrxns->{'num'}++;
#							}
#						}
#					}
					#the part after && is the check that was mentioned immediately above 
					if(exists $remain->{$sub}){
						push(@check,$sub); #this is where the problem lies. if we just put the substrate on the check stack, it means that it is already considered to be connected, even though it might need a re-traversal of one of its reactions to be made. 
						delete $remain->{$sub};
						if(exists $comps->{$sub}->{'altnum'}){
							if($comps->{$sub}->{'altnum'}>0){

								#next foreach loop is to find out which substrates need to be disconnected from $sub because the reaction being considered in this alternative can't be traversed in the opposite direction anymore.
								my $delsub;
								foreach my $alt1(keys %{$comps->{$sub}}){
									if(exists $comps->{$sub}->{$alt1}->{$cur}){
										my ($subgrp1,$rxngrp1)=split(':',$comps->{$sub}->{$alt1}->{$cur});
										foreach my $sub1(split (/\+/,$subgrp1)){
											if($sub1 ne $cur){
												$delsub->{$sub1}->{$rxngrp1}=1;
											}
										}
									}
								}
								#actual deletion of the connection talked about above after checking that we are deleting only the substrates of the same rxngrp. 
								foreach my $alt1(keys %{$comps->{$sub}}){
									foreach my $sub1(keys %{$delsub}){
										if (exists $comps->{$sub}->{$alt1}->{$sub1}){
											my ($subgrp2,$rxngrp2)=split(':',$comps->{$sub}->{$alt1}->{$sub1});
											foreach my $rxngrp1(keys %{$delsub->{$sub1}}){
												if($rxngrp2 eq $rxngrp1){
													delete $comps->{$sub}->{$alt1};
													$comps->{$sub}->{'altnum'}--;
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	if(exists $conn->{$Cs}){return (1,@connected);}
	else{return (0,@connected);}
}

sub lenient{
	my @rtrn;
	if(@tstepr>2){
		for(my $rxnstep=0;$rxnstep<@tstepr;$rxnstep++){
			my @rxngrp=split(/\,/,$tstepr[$rxnstep]);
			my $needsprxn;
			my $rxnstepflag=0;
			my $maxneeded=0;
			my $neededrxn;
			for(my $j=0;$j<@rxngrp;$j++){
				my @rxns=split(/\+/,$rxngrp[$j]);
				my $rxnflag=1;
				my $needed=0;
				my $needrxn="";
				for(my $k=0;$k<@rxns;$k++){
					#check for each remaining species whether this rxnstep is absent. if so, check each rxngrp whether forcing the presence of exactly 1 reaction can make that rxngroup complete. if so, add that species to this reaction's checkgroup and then send all species in this checkgroup along with this reaction to the sub program. after this, remove those species that are now present from checksp, print it in the .present file and run the same thing again for remaining species until all rxngrps are completed in this step. then go to the next step, assuming any species is left in checksp. after the last step, any species still left in checksp are printed in the .absent file
					$rxnflag*=$rn_data->{$rxns[$k]}->{'present'};
					if($needed==0){$needrxn=$rxns[$k]};
					$needed+=(($rn_data->{$rxns[$k]}->{'present'}-1)*(-1)); #changes 1 to 0 and 0 to 1
				}
				$rxnstepflag+=$rxnflag;
				if($needed==1){$neededrxn->{$needrxn}=0;}
			}
			if($rxnstepflag==0){
				foreach my $needrxn(keys %{$neededrxn}){
					if(exists $needsprxn->{$needrxn}->{'num'}){
						$needsprxn->{$needrxn}->{'num'}++;
					}else{$needsprxn->{$needrxn}->{'num'}=1;}
				}
			}
			if(defined $needsprxn){
				my $max=-1;
				my $maxrxn="";
				foreach my $needrxn(keys %{$needsprxn}){
					if($max<$needsprxn->{$needrxn}->{'num'}){
						$max=$needsprxn->{$needrxn}->{'num'};
						$maxrxn=$needrxn;
					}
				}
				$needsprxn->{$maxrxn}->{'num'}--;
				my $givecmd=(defined $givfile)?"-givfile $givfile":"-given $given";
				my $present=`modDFS.pl -module $module -kofile $kofile -rnfile $rnfile -modfile $modfile $givecmd -forcerxn $maxrxn -ends -last_comp $last_comp -first_comp $first_comp -recursive`;
				if($present eq '1'){push(@rtrn,$maxrxn);}
			}
		}
	}
	return @rtrn;
}


#subroutine fillcomps fills the substrate/reaction/product related hashes when given single subgroups/reactiongroups/productgroups (i.e. not separated by commas, only +'s)
sub fillcomps{
	my ($subgrp,$rxngrp,$prod,$prodgrp)=@_;
	
	my @subs=split(/\+/,$subgrp);
	my @rxns=split(/\+/,$rxngrp);
	my @prods=split(/\+/,$prodgrp);


	#filling rxncomp with substrates of reactions or reaction groups. 
	my $revfl=1; #flag that the whole reaction is reversible. it stays 1 after this for loop only if all the reactions in the +-separated group are reversible. 
	for(my $m=0;$m<@rxns;$m++){
		$revfl*=$rn_data->{$rxns[$m]}->{'rev'};
		$rxncomp->{$rxns[$m]}->{$subgrp}=1;
		if($rn_data->{$rxns[$m]}->{'rev'}==1){
			$rxncomp->{$rxns[$m]}->{$prodgrp}=1;
		}
	}
	if($rxngrp=~/\+/){ #also saving for the whole reaction group (R1+R2 etc.), though I probably never use this information. 
		$rxncomp->{$rxngrp}->{$subgrp}=1;
		if($revfl==1){$rxncomp->{$rxngrp}->{$prodgrp}=1;}
	}

	#next for loop is for collecting the information about the different reactions that involve the compounds in the subgrp (and also prods, for reversible reaction). this information is stored in the hash comprxn
	for(my $m=0;$m<@rxns;$m++){
		for my $sub(@subs){
			#products for reactions that have $sub as 'sub'
			unless(exists $comprxn->{$sub}->{'sub'}->{'rxns'}->{$rxns[$m]}){
				if(exists $comprxn->{$sub}->{'sub'}){
					$comprxn->{$sub}->{'sub'}->{'rxnnum'}++;
				}else{
					$comprxn->{$sub}->{'sub'}->{'rxnnum'}=1;
				}
				$comprxn->{$sub}->{'sub'}->{'rxns'}->{$rxns[$m]}=$prodgrp;
			}

			#if reversible, then these products are stored as substrates of these reactions which also have $sub as 'prod'
			if($rn_data->{$rxns[$m]}->{'rev'}==1){
				unless(exists $comprxn->{$sub}->{'prod'}->{'rxns'}->{$rxns[$m]}){
					if(exists $comprxn->{$sub}->{'prod'}){
						$comprxn->{$sub}->{'prod'}->{'rxnnum'}++;
					}else{
						$comprxn->{$sub}->{'prod'}->{'rxnnum'}=1;
					}
					$comprxn->{$sub}->{'prod'}->{'rxns'}->{$rxns[$m]}=$prodgrp;
				}
			}
		}

		#same thing for products now. 'prod' for normal. and 'sub' for reversible
		for my $sub(@prods){
			if($rn_data->{$rxns[$m]}->{'rev'}==1){
				unless(exists $comprxn->{$sub}->{'sub'}->{'rxns'}->{$rxns[$m]}){
					if(exists $comprxn->{$sub}->{'sub'}){
						$comprxn->{$sub}->{'sub'}->{'rxnnum'}++;
					}else{
						$comprxn->{$sub}->{'sub'}->{'rxnnum'}=1;
					}
					$comprxn->{$sub}->{'sub'}->{'rxns'}->{$rxns[$m]}=$subgrp;
				}
			}
			unless(exists $comprxn->{$sub}->{'prod'}->{'rxns'}->{$rxns[$m]}){
				if(exists $comprxn->{$sub}->{'prod'}){
					$comprxn->{$sub}->{'prod'}->{'rxnnum'}++;
				}else{
					$comprxn->{$sub}->{'prod'}->{'rxnnum'}=1;
				}
				$comprxn->{$sub}->{'prod'}->{'rxns'}->{$rxns[$m]}=$subgrp;
			}
		}
	}


	my $prsnt=1;
	for(my $m=0;$m<@rxns;$m++){
		$prsnt*=$rn_data->{$rxns[$m]}->{'present'};
	}
	if($prsnt!=0){
		for(my $m=0;$m<@subs;$m++){
			if(exists $comps->{$prod}){$comps->{$prod}->{'altnum'}++;}
			else{$comps->{$prod}->{'altnum'}=1;}
			$comps->{$prod}->{$comps->{$prod}->{'altnum'}}->{$subs[$m]}="$subgrp:$rxngrp";
			if(exists $depends->{$subs[$m]}){
				$depends->{$subs[$m]}->{'depnum'}++;
				if(exists $depends->{$subs[$m]}->{'deps'}->{$prod}){
					$depends->{$subs[$m]}->{'deps'}->{$prod}.=":$comps->{$prod}->{'altnum'}";
				}else{
					$depends->{$subs[$m]}->{'deps'}->{$prod}="$comps->{$prod}->{'altnum'}";
				}
			}else{
				$depends->{$subs[$m]}->{'depnum'}=1;
				if(exists $depends->{$subs[$m]}->{'deps'}->{$prod}){
					$depends->{$subs[$m]}->{'deps'}->{$prod}.=":$comps->{$prod}->{'altnum'}";
				}else{
					$depends->{$subs[$m]}->{'deps'}->{$prod}="$comps->{$prod}->{'altnum'}";
				}
			}
			if($revfl==1){
				if(exists $comps->{$subs[$m]}){$comps->{$subs[$m]}->{'altnum'}++;}
				else{$comps->{$subs[$m]}->{'altnum'}=1;}
				$comps->{$subs[$m]}->{$comps->{$subs[$m]}->{'altnum'}}->{$prod}="$prodgrp:$rxngrp";
				if(exists $depends->{$prod}){
					$depends->{$prod}->{'depnum'}++;
					if(exists $depends->{$prod}->{'deps'}->{$subs[$m]}){
						$depends->{$prod}->{'deps'}->{$subs[$m]}.=":$comps->{$subs[$m]}->{'altnum'}";
					}else{
						$depends->{$prod}->{'deps'}->{$subs[$m]}="$comps->{$subs[$m]}->{'altnum'}";
					}
				}else{
					$depends->{$prod}->{'depnum'}=1;
					if(exists $depends->{$prod}->{'deps'}->{$subs[$m]}){
						$depends->{$prod}->{'deps'}->{$subs[$m]}.=":$comps->{$subs[$m]}->{'altnum'}";
					}else{
						$depends->{$prod}->{'deps'}->{$subs[$m]}="$comps->{$subs[$m]}->{'altnum'}";
					}
				}
			}
		}
	}
}

