use strict;
unless(@ARGV==2){
	print "perl $0 <list> <dir>\n";
	exit;
}
my $list=$ARGV[0];
my $dir=$ARGV[1];
my %genes;
foreach my $db("ucsc","refseq","gencode"){
	open IN,"/histor/sun/wujinyu/soft/fusioncatcher/data/current/$db\_genes.txt"||die"$!";
	while(<IN>){
		chomp;
		my @e=split(/\t/,$_);
		$e[2]=$e[2]-500;
		$e[1]=$e[1]+500;
		$genes{$e[0]}="$e[4]\t$e[2]\t$e[1]";
	}
}

open LI, "$list"||die"$!";
while(<LI>){
	chomp;
	my @e=split(/\t/,$_);
	$e[0]=$_;
	my $fusionN="$dir/$e[0]/$e[0]N/final-list_candidate-fusion-genes.txt";
	my $fusionT="$dir/$e[0]/$e[0]T/final-list_candidate-fusion-genes.txt";
	my $RN="
library(chimeraviz)

fusions=importFusioncatcher(filename=\"$fusionN\",genomeVersion=\"hg38\")
#plotCircle(fusions)
edb <- ensembldb::EnsDb(\"/histor/sun/wujinyu/DB/SOAPfuse/Homo_sapiens.GRCh38.84.3.sqlite\")
";
	my $RT="
library(chimeraviz)

fusions=importFusioncatcher(filename=\"$fusionT\",genomeVersion=\"hg38\")
#plotCircle(fusions)
edb <- ensembldb::EnsDb(\"/histor/sun/wujinyu/DB/SOAPfuse/Homo_sapiens.GRCh38.84.3.sqlite\")
";
	
	`cd $dir/$e[0]/$e[0]N/;unzip -o $dir/$e[0]/$e[0]N/supporting-reads_gene-fusions_BOWTIE.zip`;
	print "$dir/$e[0]/$e[0]N\n";
	`cd $dir/$e[0]/$e[0]T/;unzip -o $dir/$e[0]/$e[0]T/supporting-reads_gene-fusions_BOWTIE.zip`;
	print "$dir/$e[0]/$e[0]T\n";
	my @fqN=`ls $dir/$e[0]/$e[0]N/*.fq`;
	my @fqT=`ls $dir/$e[0]/$e[0]T/*.fq`;
	foreach my $f(@fqN,@fqT){
		chomp($f);
		open F,"$f"||die"$!";
		my $name=$f;
		$name=~s/.fq//g;
		open O1,">$name.1.fq"||die"$!";
		open O2,">$name.2.fq"||die"$!";
		while(<F>){
			my $fq1=$_;
			$fq1.=<F>;
			$fq1.=<F>;
			$fq1.=<F>;
			print O1 "$fq1";
			my $fq2=<F>;
			$fq2.=<F>;
			$fq2.=<F>;
			$fq2.=<F>;
			print O2 "$fq2";
		}
		close(F);
	}
	open N,"$fusionN"||die"$!";
	<N>;
	my $n=0;
	my $commond;
	while(<N>){
		chomp;
		$n++;
		my @g=split(/\t/,$_);
		my @s1=split(/:/,$g[8]);
		my @s2=split(/:/,$g[9]);
		open TREG,">$dir/$e[0]/$e[0]N/$g[0]--$g[1]__$s1[1]--$s2[1].bed"||die"$!";
		print TREG "$genes{$g[0]}\n$genes{$g[1]}\n";
		close(TREG);
		print "
CrossMap.py bed ~/lxf/magic/liftOver/hg38ToHg19.over.chain.gz $dir/$e[0]/$e[0]N/$g[0]--$g[1]__$s1[1]--$s2[1].bed $dir/$e[0]/$e[0]N/$g[0]--$g[1]__$s1[1]--$s2[1].bed.hg19;
intersectBed -a $dir/../tophat/$e[0]N/accepted_hits.bam -b $dir/$e[0]/$e[0]N/$g[0]--$g[1]__$s1[1]--$s2[1].bed.hg19 -ubam >$dir/$e[0]/$e[0]N/$g[0]--$g[1]__$s1[1]--$s2[1].hg19.bam;
CrossMap.py bam ~/lxf/magic/liftOver/hg19ToHg38.over.chain.gz $dir/$e[0]/$e[0]N/$g[0]--$g[1]__$s1[1]--$s2[1].hg19.bam $dir/$e[0]/$e[0]N/$g[0]--$g[1]__$s1[1]--$s2[1].bam;\n";
		$commond.="CrossMap.py bed ~/lxf/magic/liftOver/hg38ToHg19.over.chain.gz $dir/$e[0]/$e[0]N/$g[0]--$g[1]__$s1[1]--$s2[1].bed $dir/$e[0]/$e[0]N/$g[0]--$g[1]__$s1[1]--$s2[1].bed.hg19;
intersectBed -a $dir/../tophat/$e[0]N/accepted_hits.bam -b $dir/$e[0]/$e[0]N/$g[0]--$g[1]__$s1[1]--$s2[1].bed.hg19 -ubam >$dir/$e[0]/$e[0]N/$g[0]--$g[1]__$s1[1]--$s2[1].hg19.bam;
CrossMap.py bam ~/lxf/magic/liftOver/hg19ToHg38.over.chain.gz $dir/$e[0]/$e[0]N/$g[0]--$g[1]__$s1[1]--$s2[1].hg19.bam $dir/$e[0]/$e[0]N/$g[0]--$g[1]__$s1[1]--$s2[1]
samtools index $dir/$e[0]/$e[0]N/$g[0]--$g[1]__$s1[1]--$s2[1].bam
";
		## 1. get coordinate of two fusion genes from database.
		## 2. bedtools get region reads
		## 3. CrossMap.py hg19 To hg38

		$RN.="
pdf(\"$dir/$e[0]/$e[0]N/$g[0]--$g[1]__$s1[1]--$s2[1].pdf\")
fusion=getFusionById(fusions,$n)
referenceFilename=\"$dir/$e[0]/$e[0]N/$g[0]--$g[1]__$s1[1]--$s2[1].fa\"
writeFusionReference(fusion = fusion, filename = referenceFilename)
#source(\"/panfs/home/sun/wujinyu/R/chimeraviz/scripts/bowtie.R\")
#bowtieIndex(bowtieBuildLocation =\"/panfs/home/sun/wujinyu/lxf/soft/bowtie2-2.1.0/bowtie2-build\",referenceFasta = referenceFilename)
#outputbam = \"$dir/$e[0]/$e[0]N/$g[0]--$g[1]__$s1[1]--$s2[1]\"
#fastq1=\"$dir/$e[0]/$e[0]N/$g[0]--$g[1]__$s1[1]--$s2[1]_reads.1.fq\"
#fastq2=\"$dir/$e[0]/$e[0]N/$g[0]--$g[1]__$s1[1]--$s2[1]_reads.2.fq\"
#bowtieAlign(bowtieLocation = \"/panfs/home/sun/wujinyu/lxf/soft/bowtie2-2.1.0/bowtie2\",referenceName = referenceFilename,fastq1 = fastq1,fastq2 = fastq2,outputBamFilename=outputbam)
fusion <- addFusionReadsAlignment(fusion, \"$dir/$e[0]/$e[0]N/$g[0]--$g[1]__$s1[1]--$s2[1].bam\")
#plotFusion( fusion = fusion, bamfile = \"$dir/$e[0]/$e[0]N/$g[0]--$g[1]__$s1[1]--$s2[1].bam\", edb = edb, nonUCSC = TRUE)
plotFusion( fusion = fusion, bamfile = \"$dir/$e[0]/$e[0]N/$g[0]--$g[1]__$s1[1]--$s2[1].bam\", edb = edb, nonUCSC = FALSE, reduceTranscripts = TRUE)
";
	}
	open R,">$dir/$e[0]/plot.N.R"||die"$!";
	print R "$RN
pdf(\"$dir/$e[0]/$e[0]N/fusion.pdf\")
plotCircle(fusions)";
	close R;
	close N;
	my $n=0;
	open T,"$fusionT"||die"$!";
	<T>;
	while(<T>){
        chomp;
        $n++;
        my @g=split(/\t/,$_);
		my @s1=split(/:/,$g[8]);
		my @s2=split(/:/,$g[9]);
		open TREG,">$dir/$e[0]/$e[0]T/$g[0]--$g[1]__$s1[1]--$s2[1].bed"||die"$!";
		print TREG "$genes{$g[0]}\n$genes{$g[1]}\n";
		close(TREG);
		$commond.="
CrossMap.py bed ~/lxf/magic/liftOver/hg38ToHg19.over.chain.gz $dir/$e[0]/$e[0]T/$g[0]--$g[1]__$s1[1]--$s2[1].bed $dir/$e[0]/$e[0]T/$g[0]--$g[1]__$s1[1]--$s2[1].bed.hg19;
intersectBed -a $dir/../tophat/$e[0]T/accepted_hits.bam -b $dir/$e[0]/$e[0]T/$g[0]--$g[1]__$s1[1]--$s2[1].bed.hg19 -ubam >$dir/$e[0]/$e[0]T/$g[0]--$g[1]__$s1[1]--$s2[1].hg19.bam;
CrossMap.py bam ~/lxf/magic/liftOver/hg19ToHg38.over.chain.gz $dir/$e[0]/$e[0]T/$g[0]--$g[1]__$s1[1]--$s2[1].hg19.bam $dir/$e[0]/$e[0]T/$g[0]--$g[1]__$s1[1]--$s2[1]
samtools index $dir/$e[0]/$e[0]T/$g[0]--$g[1]__$s1[1]--$s2[1].bam
";
        $RT.="
pdf(\"$dir/$e[0]/$e[0]T/$g[0]--$g[1]__$s1[1]--$s2[1].pdf\")
fusion=getFusionById(fusions,$n)
referenceFilename=\"$dir/$e[0]/$e[0]T/$g[0]--$g[1]__$s1[1]--$s2[1].fa\"
writeFusionReference(fusion = fusion, filename = referenceFilename)
#source(\"/panfs/home/sun/wujinyu/R/chimeraviz/scripts/bowtie.R\")
#bowtieIndex(bowtieBuildLocation =\"/panfs/home/sun/wujinyu/lxf/soft/bowtie2-2.1.0/bowtie2-build\",referenceFasta = referenceFilename)
#outputbam = \"$dir/$e[0]/$e[0]T/$g[0]--$g[1]__$s1[1]--$s2[1]\"
#fastq1=\"$dir/$e[0]/$e[0]T/$g[0]--$g[1]__$s1[1]--$s2[1]_reads.1.fq\"
#fastq2=\"$dir/$e[0]/$e[0]T/$g[0]--$g[1]__$s1[1]--$s2[1]_reads.2.fq\"
#bowtieAlign(bowtieLocation = \"/panfs/home/sun/wujinyu/lxf/soft/bowtie2-2.1.0/bowtie2\",referenceName = referenceFilename,fastq1 = fastq1,fastq2 = fastq2,outputBamFilename=outputbam)
fusion <- addFusionReadsAlignment(fusion, \"$dir/$e[0]/$e[0]T/$g[0]--$g[1]__$s1[1]--$s2[1].bam\")
#plotFusion( fusion = fusion, bamfile = \"$dir/$e[0]/$e[0]T/$g[0]--$g[1]__$s1[1]--$s2[1].bam\", edb = edb, nonUCSC = TRUE)
plotFusion( fusion = fusion, bamfile = \"$dir/$e[0]/$e[0]T/$g[0]--$g[1]__$s1[1]--$s2[1].bam\", edb = edb, nonUCSC = FALSE, reduceTranscripts = TRUE)
";
	}
	close T;
	open R,">$dir/$e[0]/plot.T.R"||die"$!";
	print R "$RT
pdf(\"$dir/$e[0]/$e[0]T/fusion.pdf\")
plotCircle(fusions)";
	close(R);
	`echo "#!/bin/bash
#PBS -N fusion
#PBS -l nodes=node77:ppn=1
#PBS -j oe
#PBS -q superfat
cd $dir/$e[0]/
$commond
R CMD BATCH $dir/$e[0]/plot.N.R
R CMD BATCH $dir/$e[0]/plot.T.R
" >$dir/$e[0]/$e[0].sh;qsub $dir/$e[0]/$e[0].sh`;
}
