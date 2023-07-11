
rule all5:
   input: "/home/alexclear/HIC_amr_plasmids_project/sputum_analysis/report/amr_pie_SD8.pdf"



rule report:
    output:  
         percent_mapped = "/home/alexclear/HIC_amr_plasmids_project/sputum_analysis/report/mapping.pdf",
         amr_piechart1 = "/home/alexclear/HIC_amr_plasmids_project/sputum_analysis/report/amr_pie_SD2.pdf",
         amr_piechart2 = "/home/alexclear/HIC_amr_plasmids_project/sputum_analysis/report/amr_pie_SD8.pdf"
    conda:"tidyverse_env"
    shell: '''
    Rscript /home/alexclear/HIC_amr_plasmids_project/scripts/sputum_plots.R -f "SD002" -s "SD008" -m {output.percent_mapped} -p {output.amr_piechart1} -o {output.amr_piechart2}
    '''