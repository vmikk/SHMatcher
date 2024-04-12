## Script to linearize FASTA file

NR==1 {print ; next} {printf /^>/ ? "\n"$0"\n" : $1} END {printf "\n"}
