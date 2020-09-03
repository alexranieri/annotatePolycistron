#!/usr/bin/env bash

# This script annotates policistrons.

# Variables are better than hard-coded values or magic numbers.
E_WRONGARG=1      # error return code when the script is called without one argument
INPUT=${1}        # input GFF filename
CHROM=""          # chromosome on the current line being processed
TYPE=""           # type on the current line being processed
STRAND=""         # strand on the current line being processed
PC_CHROM=""       # chromosome of the current policistron being annotated
PC_TYPE=""        # type of the current policistron being annotated
PC_START=""       # chromStart of the current policistron being annotated
PC_END=""         # chromEnd of the current policistron being annotated
PC_STRAND=""      # strand of the current policistron being annotated
PC_PHASE=""       # phase of the current policistron being annotated
PC_DESCRIPTION="" # description of the current policistron being annotated
NUM_PC_CDS="0"    # unique identifier of policistron of CDS
NUM_PC_NCRNA="0"  # unique identifier of policistron of ncRNA
NUM_PC_RRNA="0"   # unique identifier of policistron of rRNA
NUM_PC_SNORNA="0" # unique identifier of policistron of snoRNA
NUM_PC_TRNA="0"   # unique identifier of policistron of tRNA


# This function prints the current chromosome.
function getChromosome {
   echo "${LINE}" | cut -f1
}


# This function prints the current type.
function getType {
   echo "${LINE}" | cut -f3
}


# This function prints the current chromStart.
function getChromStart {
   echo "${LINE}" | cut -f4
}


# This function prints the current chromEnd.
function getChromEnd {
   echo "${LINE}" | cut -f5
}


# This function prints the current strand.
function getStrand {
   echo "${LINE}" | cut -f7
}


# This function prints the current phase.
function getPhase {
   echo "${LINE}" | cut -f8
}


# This function verifies if the type of the current line is one of the target types.
function processEventualPolicistron {
   # Printing the current policistron, if there is one to be printed.
   if [ "${PC_CHROM}" != "" ]
   then
      echo -e "${PC_CHROM}\tannotatePolicistron\tpolicistron-${PC_TYPE}\t${PC_START}\t${PC_END}\t.\t${PC_STRAND}\t${PC_PHASE}\t${PC_DESCRIPTION}"
   fi

   # Check eventual policistron start.
   TYPE=`getType`
   case "${TYPE}" in
      "CDS" | "ncRNA" | "rRNA" | "snoRNA" | "tRNA")
         updatePCData
         ;;
      *)
         # Reset policistron data.
         PC_CHROM=""
         ;;
   esac
}


# This function updates the data related to a policistron.
function updatePCData {
   PC_CHROM=`getChromosome`
   PC_TYPE=`getType`
   PC_START=`getChromStart`
   PC_END=`getChromEnd`
   PC_STRAND=`getStrand`
   PC_PHASE=`getPhase`

   # Updating the description and increasing by one the number of policistron of the identified type.
   case "${PC_TYPE}" in
      "CDS")
         PC_DESCRIPTION="ID=policistron-${PC_TYPE}_${NUM_PC_CDS}"
         let "NUM_PC_CDS = NUM_PC_CDS + 1"
         ;;
      "ncRNA")
         PC_DESCRIPTION="ID=policistron-${PC_TYPE}_${NUM_PC_NCRNA}"
         let "NUM_PC_NCRNA = NUM_PC_NCRNA + 1"
         ;;
      "rRNA")
         PC_DESCRIPTION="ID=policistron-${PC_TYPE}_${NUM_PC_RRNA}"
         let "NUM_PC_RRNA = NUM_PC_RRNA + 1"
         ;;
      "snoRNA")
         PC_DESCRIPTION="ID=policistron-${PC_TYPE}_${NUM_PC_SNORNA}"
         let "NUM_PC_SNORNA = NUM_PC_SNORNA + 1"
         ;;
      "tRNA")
         PC_DESCRIPTION="ID=policistron-${PC_TYPE}_${NUM_PC_TRNA}"
         let "NUM_PC_TRNA = NUM_PC_TRNA + 1"
         ;;
      *)
         ;;
   esac
}


# Checking the number of arguments.
if [[ $# -ne 1 ]]
then
   echo "You should pass a GFF file as argument."
   echo "Usage: `basename ${0}` <file.gff>"
   echo "Aborting."
   exit ${E_WRONGARG}
fi


# Sort text using case-sensitive comparisons.
LC_ALL=C


# Sort by chromosome (or scaffold, or transcript, column 1), then sort numerically by chromStart
# (column 4) and finally sort numerically by chromEnd (column 5).
SORTED_INPUT=`tempfile`
< ${INPUT} sort -k1,1 -k4n,4 -k5n,5 > ${SORTED_INPUT}


# In order to process each line as a unit, we have to set the input field separator as a '\n'.
IFS='
'


# Start processing the sorted file.
for LINE in `cat "${SORTED_INPUT}"`
do
   # We process a line only if it is not a comment (i. e., does not start with the character '#').
   if [ ${LINE:0:1} != "#" ]
   then
      # Read the chromosome.
      CHROM=`getChromosome`

      # If the chromosome is different, we have an eventual policistron boundary.
      if [ "${CHROM}" != "${PC_CHROM}" ]
      then
         processEventualPolicistron
      else
         # Checking if the type is the same of the current policistron.
         TYPE=`getType`
         if [ "${TYPE}" != "${PC_TYPE}" ] && [ "${TYPE}" != "exon" ] && [ "${TYPE}" != "gene" ] && [ "${TYPE}" != "mRNA" ]
         then
            processEventualPolicistron
         fi

         # Read the strand.
         STRAND=`getStrand`

         if [ "${TYPE}" = "${PC_TYPE}" ]
         then
            # If the strand is different, we also have an eventual policistron boundary.
            if [ "${STRAND}" != "${PC_STRAND}" ]
            then
               processEventualPolicistron
            # At this point we know that the current chromosome, type and strand are the same of the
            # policistron being annotated. This means that a new feature of the same type composes
            # the policistron and we have to update the chromEnd in order to reflect the extension.
            else
               PC_END=`getChromEnd`
            fi
         fi
      fi
   fi

   # Printing the current line. This way, all input lines are printed. If this line indicates the
   # end of a policistron, the later has already been printed by the above code.
   echo "${LINE}"
done

# TODO: Como tratar o último policistron?

rm "${SORTED_INPUT}"

# Test if it is a target feature.
# If it is a target feature, verify if the next feature is of the same
# category.
# If the category is the same, verify if the strand is the same.
# Any negative to the above conditions implies the end of a policistron. At this moment, we should
# annotate the policistron, the SSR (which can be a DSSR (divergent) or a CSSR (convergent) and the
# TSS (or TTS, which can be unique, depending of the distance).
# If the policistron is a  monocistron, which fields should have at ninth column to identify this
# issue?
# At policistron annotation, we should list its components.
# At SSR annotation, we should explicit the TSS or TTS that are part of it. Maybe it would be
# interesting to also identify the neighbor policistrons, in the case they are monocistrons.

exit 0

my @array_1;
my @array_2;
my $direction;
my $type;
my $start; 
my $end;
my $contig_id;
my $size;

# Lê o arquivo linha por linha
while ($linha) {
	chomp($linha);
	# Armazena os atributos da primeira linha	
	@array_1 = split(/\t/,$linha); 
	$direction = $array_1[6];
	$type = $array_1[2];
	$start = $array_1[3];
	$contig_id = $array_1[0];
	@features = ();

	# Os blocos abaixo recuperam o id da feature
	my @annotations = split(";",$array_1[8]);
	foreach my $annotation (@annotations) {
		if ($annotation =~ /\bid=\b/){    
		$feature_id = substr($annotation, 3, length($annotation)-3)."\n";
		push @features, $feature_id;
		last;
		}
	}
	# Fim de recuperação de id da feature

	$end = $array_1[4];
	$size = $end-$start;

	# Agora vamos ler a próxima linha e continuar lendo enquanto for o mesmo tipo, direção e contig id
	my $same = "yes"; 
	while ($linha && $same eq "yes") {
		$linha = readline($input);
		chomp($linha);
		my @array_2 = split(/\t/,$linha); 
		# Checa as condições
		if (($direction eq $array_2[6]) && ($type eq $array_2[2]) && ($contig_id eq $array_2[0])) {
			$end = $array_2[4];
			# Os blocos abaixo recuperam o id da feature
			my @annotations = split(";",$array_2[8]);
			foreach my $annotation (@annotations) {
				if ($annotation =~ /\bid=\b/){    
				$feature_id = substr($annotation, 3, length($annotation)-3)."\n";
				push @features, $feature_id;
				last;
				}
			}
			# Fim de recuperação de id da feature	
			$size = $end-$start;		
		} 
		else {
			#print("Entrou parada de policistron\n");
			$same = "no"			
	  	last; # Aqui paro o loop caso as features das linhas não sejam iguais			
		}
	}

	if ($type eq "CDS") {
		$CDS_count++;
		# Escrevendo no arquivo
		chomp @features;
		if ($direction eq "positive") {
			$direction = "+"; }
		elsif ($direction eq "negative") {
			$direction = "-";}
		else {
			$direction = ".";}
		my $data = $contig_id."\t"."write_policistron.pl"."\t"."Policistron CDS"."\t".$start."\t".$end."\t"."."."\t".$direction."\t"."."."\t"."Policistron CDS_".$CDS_count."\n";
		print $output $data;
		#$linha = readline($input);
		#chomp($linha);	
		#print("entrou CDS\n");	
	}

	elsif ($type eq "tRNA") {
		$tRNA_count++;
		# Escrevendo no arquivo
		chomp @features;
		if ($direction eq "positive") {
			$direction = "+"; }
		elsif ($direction eq "negative") {
			$direction = "-";}
		else {
			$direction = ".";}
		my $data = $contig_id."\t"."write_policistron.pl"."\t"."Policistron tRNA"."\t".$start."\t".$end."\t"."."."\t".$direction."\t"."."."\t"."Policistron tRNA_".$tRNA_count."\n";	
		print $output $data;
		#$linha = readline($input);
		#chomp($linha);
		#print("entrou tRNA\n");
	}

	elsif ($type eq "rRNA") {
		$rRNA_count++;
		# Escrevendo no arquivo
		chomp @features;
		if ($direction eq "positive") {
			$direction = "+"; }
		elsif ($direction eq "negative") {
			$direction = "-";}
		else {
			$direction = ".";}
		my $data = $contig_id."\t"."write_policistron.pl"."\t"."Policistron rRNA"."\t".$start."\t".$end."\t"."."."\t".$direction."\t"."."."\t"."Policistron rRNA_".$rRNA_count."\n";	
		print $output $data;
		#$linha = readline($input);
		#chomp($linha);
		#print("entrou rRNA\n");
	}

	else {
		$ncRNA_count++;
		# Escrevendo no arquivo
		chomp @features;
		if ($direction eq "positive") {
			$direction = "+"; }
		elsif ($direction eq "negative") {
			$direction = "-";}
		else {
			$direction = ".";}
		my $data = $contig_id."\t"."write_policistron.pl"."\t"."Policistron ncRNA"."\t".$start."\t".$end."\t"."."."\t".$direction."\t"."."."\t"."Policistron ncRNA_".$ncRNA_count."\n";	
		print $output $data;
		#$linha = readline($input);
		#chomp($linha);
		#print("entrou ncRNA\n");
	} 


}

my $soma = $CDS_count + $tRNA_count + $rRNA_count + $ncRNA_count;
print("A quantidade de policistrons obtidos foi: ".$soma."\n");
print("O total para cada categoria foi de: \n");
print("CDS: ".$CDS_count."\n");
print("tRNA: ".$tRNA_count."\n");
print("rRNA: ".$rRNA_count."\n");
print("ncRNA: ".$ncRNA_count."\n");
