prefetch --option-file /data/cephfs-1/scratch/groups/dubois/users/rauertc_c/scRNAseq_LPS/SRR_Acc_List.txt --output-directory /data/cephfs-1/scratch/groups/dubois/users/rauertc_c/scRNAseq_LPS/SRA
for dir in /data/cephfs-1/scratch/groups/dubois/users/rauertc_c/scRNAseq_LPS/SRA/*/
do
    fasterq-dump $dir
done
