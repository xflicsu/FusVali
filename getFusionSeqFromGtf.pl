use strict;
unless(@ARGV==3){
	print "perl $0 <fusion.result> <dir> <name>\n";
	exit;
}

my $outdir=$ARGV[1];
my $name=$ARGV[2];
open GTF,"/histor/sun/wujinyu/soft/fusioncatcher/data/current/organism.gtf"||die"$!";
my %trans;
my %strand;
my %exonSite;
while(<GTF>){
	next if /^#/;
	s/"//g;
	chomp;
	my @e=split(/\t/,$_);
	next if $e[2] ne "exon";
	$e[8]=~/exon_id (ENSE\d+?);/;
	my $exon=$1;
	my @infor=split(/\s+/,$e[8]);
	$e[3]-=1;
	$strand{$exon}=$e[6];
	$exonSite{$exon}="chr$e[0]\t$e[3]\t$e[4]\t$exon\t.\t$e[6]";
	##ENSG00000227232 ENST00000488147 1 region
}

my %fa;
open FA,"/histor/sun/wujinyu/soft/fusioncatcher/data/current/transcripts.fa"||die"$!";
$/ = ">"; <FA>; $/ = "\n";
while (<FA>){
	chomp;
	my $tr = (split ( /;/, $_))[0];
	$/ = ">";
	my $seq = <FA>;
	$/ = "\n";
	chomp($seq);
	$seq =~ s/\s+//g;
	$seq =~ s/>$//;
	$seq = uc($seq);
	$fa{$tr}=$seq;
}

open FUS,"$ARGV[0]"||die"$!";
<FUS>;
open OUT,">$outdir/$name.fusion.for.PCR.sequence.txt"||die"$!";
print OUT "Gene_1_symbol(5end_fusion_partner)\tGene_2_symbol(3end_fusion_partner)\tPredicted_effect\tFusion_description\tPredicted_trancripts_or_exons\tFusion_sequence\tFusion_sequence_extended\n";
while(<FUS>){
	chomp;
	my @e=split(/\t/,$_);
	# 5'end gene,3'end gene,5'end exon,3'end exon,Predicted_fused_transcripts
	$e[8];
	$e[9];
	$e[10];
	$e[11];
	$e[12];
	$e[13];
	# transcripts
	if($e[16] ne ""){
		my @t=split(/;/,$e[16]);
		foreach my $ts(@t){
			$ts=~/^(ENST\d+?):(\d+?)\/(ENST\d+?):(\d+)/;
			my $t5=$1;
			my $t5len=$2;
			my $t3=$3;
			my $t3len=$4;
			print "xxx\t$t3len\n";
			my $t5reg=substr($fa{$t5},0,$t5len);
			my $t3reg=substr($fa{$t3},$t3len-1);
			my $str="$t5reg*$t3reg";
			print OUT "$e[0]\t$e[1]\t$e[15]\t$e[3]\t$ts\t$e[14]\t$str\n";
		}
	}else{
		my $reg5;
		my $reg3;
		### 5 end
		if($strand{$e[12]} eq "+"){
			my @break=split(/:/,$e[8]);
			my @lastExon=split(/\t/,$exonSite{$e[12]});
			print "ssssssssss\t$break[1]\t$lastExon[1]\n";
			if($break[1]<$lastExon[2]){
				print "$break[1]<$lastExon[2]\n";
			}
			$lastExon[2]=$break[1];
			my $r=join "\t",@lastExon;
			$reg5="$r\n";
		}else{
			my @break=split(/:/,$e[8]);
			my @lastExon=split(/\t/,$exonSite{$e[12]});
			print "ooooooo\t$break[1]\t$lastExon[1]\n";
			if($break[1]>$lastExon[1]){
				print "$break[1]<$lastExon[1]\n";
			}
			$lastExon[1]=$break[1]-1;
			my $r=join "\t",@lastExon;
			$reg5="$r\n";
		}
		open O,">$outdir/$e[12].BED"||die"$!";
		print O "$reg5";
		close(O);
		`fastaFromBed -fi /histor/sun/wujinyu/DB/SOAPfuse/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa -bed $outdir/$e[12].BED -s -tab -fo $outdir/$e[12].fa`;
		### 3 end
		if($strand{$e[13]} eq "+"){
			my @break=split(/:/,$e[9]);
			my @lastExon=split(/\t/,$exonSite{$e[13]});
			if($break[1]>$lastExon[2]){
				print "$break[1]>$lastExon[2]\n";
			}
			$lastExon[1]=$break[1]-1;
			my $r=join "\t",@lastExon;
			$reg3="$r\n";
		}else{
			my @break=split(/:/,$e[9]);
			my @lastExon=split(/\t/,$exonSite{$e[13]});
			if($break[1]<$lastExon[2]){
				print "$break[1]<$lastExon[2]\n";
			}
			$lastExon[2]=$break[1]-1;
			my $r=join "\t",@lastExon;
			$reg3="$r\n";
		}
		open O,">$outdir/$e[13].BED"||die"$!";
		print O "$reg3";
		close(O);
		`fastaFromBed -fi /histor/sun/wujinyu/DB/SOAPfuse/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa -bed $outdir/$e[13].BED -s -tab -fo $outdir/$e[13].fa`;
		open S5,"$outdir/$e[12].fa"||die"$!";
		my $s5=<S5>;
		chomp($s5);
		my @s5e=split(/\t/,$s5);
		open S3,"$outdir/$e[13].fa"||die"$!";
		my $s3=<S3>;
		chomp($s3);
		my @s3e=split(/\t/,$s3);
		my $str="$s5e[1]*$s3e[1]";
		print OUT "$e[0]\t$e[1]\t$e[15]\t$e[3]\t$e[12]/$e[13]\t$e[14]\t$str\n";
	}
}
