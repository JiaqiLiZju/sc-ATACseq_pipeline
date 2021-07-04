cat hg38.Refseq.bed |sed "s/chr/HUMAN_/g" >Mix.Human.Refseq.bed
cat mm10.Refseq.bed |sed "s/chr/MOUSE_/g" >Mix.Mouse.Refseq.bed