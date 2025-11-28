
awk '$1 ~ /^#/ || ($5=="<DUP>" && $7=="PASS")' gnomad.v4.1.sv.sites.vcf > gnomad.v4.1.sv.pass.DUP.vcf
awk '$1 ~ /^#/ || ($5=="<DEL>" && $7=="PASS")' gnomad.v4.1.sv.sites.vcf > gnomad.v4.1.sv.pass.DEL.vcf

awk -F'\t' '
BEGIN {
    OFS="\t";
    print "#CHROM","POS","ID","QUAL","AF","SVTYPE","CHR2","END","PREDICTED_LOF"
}
!/^#/ {
    AF="."; SVTYPE="."; CHR2="."; END_POS="."; LOF=".";
    n=split($8,info,";");
    for(i=1;i<=n;i++){
        split(info[i],x,"=");
        if(x[1]=="AF") AF=x[2];
        else if(x[1]=="SVTYPE") SVTYPE=x[2];
        else if(x[1]=="CHR2") CHR2=x[2];
        else if(x[1]=="END") END_POS=x[2];
        else if(x[1]=="PREDICTED_LOF") LOF=x[2];
    }
    print $1,$2,$3,$6,AF,SVTYPE,CHR2,END_POS,LOF
}' gnomad.v4.1.sv.pass.DUP.vcf > gnomad.v4.1.sv.pass.DUP.tsv
