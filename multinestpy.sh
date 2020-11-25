#!/bin/bash

fittingcode=/home/astro/variu/phd/voids/chengscodes/BAOfit/voidnw_FFTlog_myVer/baofit.py
plotfile=/home/astro/variu/phd/voids/chengscodes/BAOfit/voidnw_FFTlog_myVer/analyse/plot.py
statsfile=/home/astro/variu/phd/voids/chengscodes/BAOfit/voidnw_FFTlog_myVer/analyse/stats.py

inpath_mocks_path=/scratch/variu/clustering/patchy_cmass_subset/box1/real/overlapping/vv2pcf/
inpath_fit_path=/scratch/variu/clustering/patchy_cmass_subset/box1/real/overlapping/vv2pcf/
outpath=/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/16R/void_G1024CIC_60_150_mocks/
templatepath_arg=/scratch/variu/clustering/linearnw_CIC/patchy_cmass_like/fix_ampl/box1/G1024/

# model=parab or pvoid
model=pvoid

#avgfile=True or False
avgfile_bool=True

#RminARR=(4 6 8  10 12 14 16 18 20 22 16)
#RmaxARR=(6 8 10 12 14 16 18 20 22 30 NONE)

RminARR=(16)
RmaxARR=(NONE)
N=1

#mode=COMBINED or NONCOMBINED
mode=NONCOMBINED

###################### Check which model is used.
if [[ "${model}" == pvoid ]]; then
    templatepath=${templatepath_arg}
else 
    templatepath=NONE
fi

##### Test for Code path, Plot, Stats, InPath and OutPath
if [[ -f "$fittingcode" ]]; then
    echo "Using the fitting file from ${fittingcode}"
else
    echo "ERROR: The path to the code does not exist"
    echo "${fittingcode}"
    exit 1
fi

if [[ -f "$plotfile" ]]; then
    echo "Using the plot file from ${plotfile}"
else
    echo "ERROR: The path to the plot file does not exist"
    echo "${plotfile}"
    exit 1
fi

if [[ -f "$statsfile" ]]; then
    echo "Using the plot file from ${statsfile}"
else
    echo "ERROR: The path to the stats file does not exist"
    echo "${statsfile}"
    exit 1
fi

if [[ -d "$outpath" ]]; then
    echo "The path to the output folder is ${outpath}"
else
    echo "ERROR: The path to the output folder does not exist"
    echo "${outpath}"
    exit 1
fi

if [[ -d "$inpath_mocks_path" ]]; then
    echo "The path to the list of mocks is ${inpath_mocks_path}"
else
    echo "ERROR: The path to the list of mocks does not exist"
    echo "${inpath_mocks_path}"
    exit 1
fi

if [[ -d "$inpath_fit_path" ]]; then
    echo "The path to the fitting file is ${inpath_fit_path}"
else
    echo "ERROR: The path to the fitting file does not exist"
    echo "${inpath_fit_path}"
    exit 1
fi

if [[ "${model}" == pvoid ]]; then
    if [[ -d "$templatepath" ]]; then
        echo "The path to the template folder is ${templatepath}"
    else
        echo "ERROR: The path to the template folder does not exist"
        echo "${templatepath}"
        exit 1
    fi
else
    echo "The code uses the PARABOLIC or Galaxy model"
fi

#####

for (( i=0; i<N; i++ ))
do
    Rmin=${RminARR[$i]}
    Rmax=${RmaxARR[$i]}
  
    if [[ "${Rmin}" == NONE && "${Rmax}" == NONE ]]; then
        limits=""
        echo "Write the output files directly in ${outpath}/${limits}"
    else

        if [[ "${Rmin}" == NONE ]]; then            
            echo R${Rmax}
            limits=R${Rmax}
        else
            if [[ "${Rmax}" == NONE ]]; then
                echo ${Rmin}R
                limits=${Rmin}R
            else
                echo ${Rmin}R${Rmax}
                limits=${Rmin}R${Rmax}
            fi
        fi
    fi
    
	if [[ ! -d "${inpath_mocks_path}/${limits}" ]]; then
        echo "ERROR: The INPUT folder ${inpath_mocks_path}/${limits} does not exist."
        exit 1
    fi

    if [[ ! -d "${inpath_fit_path}/${limits}" ]]; then
        echo "ERROR: The INPUT folder ${inpath_fit_path}/${limits} does not exist."
        exit 1
    fi

    if [[ "${model}" == pvoid ]]; then
        if [[ ! -d "${templatepath}/${limits}" ]]; then
            echo "ERROR: The TEMPLATE folder ${templatepath}/${limits} does not exist."
            exit 1
        fi
    fi
    

	if [[ "${mode}" == COMBINED ]]; then
        echo "Not Implemented"
        exit 1
		# for indir in ${inpath}/${limits}/w*
		# do
		# 	if [[ ! -d "${indir}" ]]; then
		# 		echo "ERROR: The INPUT folder ${indir} does not exist."
		# 		exit 1
		# 	fi
		# 	if [[ ! -f "${indir}/avg_${limits}_$(basename $indir).2pcf" ]]; then
		# 		echo "ERROR: The file $indir/avg_${limits}_$(basename $indir).2pcf does not exist."
		# 		exit 1
		# 	fi
        #     mockpaths=${indir}/500_${limits}_$(basename $indir).dat
		# 	if [[ ! -f "${mockpaths}" ]]; then
		# 		echo "ERROR: The file ${mockpaths} does not exist."
		# 		exit 1
		# 	fi
			
		# 	echo python ${fittingcode} ${indir}/avg_${limits}_$(basename $indir).2pcf ${mockpaths} >> ../temp_multinest_fit.txt
        #     echo python ${plotfile} ${outpath}/BAOfit_avg_${limits}_$(basename $indir).2pcf_ >> ../temp_plot_fit.txt
		# done
	else
		if [[ "${mode}" == NONCOMBINED ]]; then
			indir=${inpath_fit_path}/${limits}
            if [[ "${avgfile_bool}" == True ]]; then
                if [[ ! -f "${indir}/avg_${limits}.2pcf" ]]; then
                    echo "ERROR: The file $indir/avg_${limits}.2pcf does not exist."
                    exit 1
                else
                ### Possible modification
                    fitfile=${indir}/avg_${limits}.2pcf
                fi
            else
                if [[ ! -f "${indir}/voids_Box_z0.465600_${limits}.2pcf" ]]; then
                    echo "ERROR: The file $indir/voids_Box_z0.465600_${limits}.2pcf does not exist."
                    exit 1
                else
                ### Possible modification
                    fitfile=${indir}/voids_Box_z0.465600_${limits}.2pcf
                fi
            fi 
            ### Possible modification
            mockpaths=${inpath_mocks_path}/${limits}/500_${limits}.dat
            if [[ ! -f "${mockpaths}" ]]; then
				echo "ERROR: The file ${mockpaths} does not exist."
				exit 1
			fi

            if [[ "${model}" == pvoid ]]; then
                if [[ ! -f "${templatepath}/${limits}/avg_${limits}_extended.pspec" ]]; then
                    echo "ERROR: The TEMPLATE file ${templatepath}/${limits}/avg_${limits}_extended.pspec does not exist."
                    exit 1
                else
                    templatefile=${templatepath}/${limits}/avg_${limits}_extended.pspec
                fi
            else
                templatefile=NONE
            fi
        
            # for infmock in ${inpath}/${limits}/CA*2pcf
            # do
            #     #echo python ${fittingcode} ${outpath} ${infmock} ${mockpaths} ${templatepath}/${limits}/avg_${limits}_extended.pspec >> ../temp_multinest_fit.txt
            #     echo python ${fittingcode} ${outpath} ${infmock} ${mockpaths} NONE >> ../temp_multinest_fit.txt
			# done
            #echo python ${fittingcode} ${outpath} ${inpathMD}/${limits}/voids_Box_z0.465600_${limits}.2pcf ${mockpaths} NONE >> ../temp_multinest_fit.txt
            #echo python ${fittingcode} ${outpath} ${indir}/avg_${limits}.2pcf ${mockpaths} NONE >> ../temp_multinest_fit.txt
			
            echo python ${fittingcode} ${outpath} ${fitfile} ${mockpaths} ${templatefile} >> ./temp_multinest_fit.txt
			
            #echo python ${plotfile} ${outpath}/BAOfit_voids_Box_z0.465600_${limits}.2pcf_ >> ../temp_plot_fit.txt
            #echo python ${plotfile} ${outpath}/BAOfit_avg_${limits}.2pcf_ >> ../temp_plot_fit.txt
		else
			echo "ERROR: Wrong mode. You should put either COMBINED or NONCOMBINED"
			exit 1
		fi
	fi
done
echo "Successfully Done!"
