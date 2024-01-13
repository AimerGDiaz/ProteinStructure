#bash metrics.sh P3_polimerization/P3/p3-tetramer/ 3  
#
module load GCCcore/10.2.0 GCCcore/9.3.0 Python/3.8.6
module load intel-compilers/2021.4.0  impi/2021.4.0 matplotlib/3.4.3
field=$2
Folder_name=`echo $1 | awk -F '/' '{print $'"$field"'}'`
mkdir $1/$Folder_name"_results" 2>/dev/null
best_model=`grep -A 1 order $1/ranking_debug.json | tail -n 1 | tr -d '"| |,'`
models=`echo $best_model  | awk -F'_multimer' '{print $1}'`
version=`echo $best_model  | awk -F'v[0-9]_|.pdb' '{print $2}'`
cp $1/*$models*pdb $1/$Folder_name"_results"
cp $1/features.pkl $1/result_"$models"_multimer_v3_"$version".pkl $1/$Folder_name"_results"
python visualize_alphafold_results.py --input_dir $1/$Folder_name"_results" --output_dir $1/$Folder_name"_results" --name PKL

gzip $1/$Folder_name"_results"/result_"$models"_multimer_v3_"$version".pkl

##iptm+ptm
iptm=`grep $best_model $1/ranking_debug.json  | head -n1 |  awk -F' ' '{ gsub(",","",$2); printf "%.0f\n", $2*100}'`
plddt=`awk -v i=1 '{sum+=$11; i=i+1 }END{printf ("%d\n", sum/i)}'  $1/ranked_0.pdb`
plddt_ur=`awk -v i=1 '{sum+=$11; i=i+1 }END{printf ("%d\n", sum/i)}'  $1/unrelaxed_"$models"_multimer_v3_"$version".pdb`

File_name=`echo $Folder_name |  tr '[:upper:]' '[:lower:]'`
ms=`echo $models | awk -F'_' '{print $2}' `
vs=`echo $version | awk -F'_' '{print $2}' `
mv $1/$Folder_name"_results"/unrelaxed_"$models"_multimer_v3_"$version".pdb $1/$Folder_name"_results"/$File_name"_unrelax_m"$ms"_p"$vs"_iptm-"$iptm"_plddt-"$plddt_ur".pdb" 
mv $1/$Folder_name"_results"/relaxed_"$models"_multimer_v3_"$version".pdb $1/$Folder_name"_results"/$File_name"_relax_m"$ms"_p"$vs"_iptm-"$iptm"_plddt-"$plddt".pdb" 

tar  -czf $Folder_name"_results".tar.gz $1/$Folder_name"_results" $1/msas  

#for f  in CelegansHexa/SpastinCelegansHexa/*model_4*pdb CelegansHexa/SpastinCelegansHexa/ranked_0.pdb; do echo -ne "$f \t"; awk -v i=1 '{sum+=$11; i=i+1 }END{print ("%.0%.0fn",sum/i)}' $f ; done
#relaxed_model_4_multimer_v3_pred_0.pdb 	65.1175
#unrelaxed_model_4_multimer_v3_pred_0.pdb 	64.9047
#unrelaxed_model_4_multimer_v3_pred_1.pdb 	65.418
#unrelaxed_model_4_multimer_v3_pred_2.pdb 	64.9367
#unrelaxed_model_4_multimer_v3_pred_3.pdb 	65.224
#unrelaxed_model_4_multimer_v3_pred_4.pdb 	65.1054
