import csv
from getconservativeness import conservativeness_estimator


genes = ['HSPB8']
gene = 'HSPB8'
window_size = 50
transitions = (0.5, 0.5)


rule estimate_conservativeness:
    input: genefile=expand("{gene}.fa", gene=genes),
        variants_file= expand("{gene}_variants.tsv", gene=genes)
    output: expand("output/{gene}_final_tab.tsv", gene=genes)
    run:
        link_to_genefile = f'{input.genefile}'
        link_to_variants_file = f'{input.variants_file}'
        gene_name, coords_list, u_list, n_list, mean_allele_number_list, path = conservativeness_estimator(genefilename=link_to_genefile,
                                                  variantsfilename=link_to_variants_file,
                                                   window_size=window_size,transitions=transitions)
        with open(f'{output}','w') as out_f:
            writer = csv.writer(out_f,delimiter='\t')
            writer.writerow(['Gene', 'Start', 'Stop', 'sumFreqLoF', 'sumAC',
                             'meanAN', 'State'])
            for i in range(len(coords_list)):
                writer.writerow([gene, coords_list[i][0], coords_list[i][1], u_list[i], n_list[i],
                                 mean_allele_number_list[i], path[i]])

rule generate_plofffinder_input:
    input: "data/variants_coords_ac_an_name.tsv",
    output: f"{gene}_variants.tsv"
    shell: 'grep -w "{gene}" {input} > {output}'

rule generate_freqcounter_input:
    input: "data/GRCh37.p13.genome.fa", f"snake_{gene}_cds.gff"
    output: f'{gene}.fa'
    shell: "bedtools getfasta -fi {input[0]} -bed {input[1]} -name > {output}"

rule generate_cds_gff:
    output: temp(f"snake_{gene}_cds.gff")
    shell: "sh get_cds.sh {gene} > {output}"
