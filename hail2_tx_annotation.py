import hail as hl
import hail.expr.aggregators as agg
from gnomad_hail import *
hl.init()

# MTs of interest
gnomad_vepped_mt_path = "gs://gnomad-berylc/tx-annotation/hail2/gnomad.exomes.r2.0.2.sites.split.vep.030818.mt"
exac_vepped_vds_path = "gs://gnomad/raw/hail-0.2/mt/exac/exac.r1.sites.vep.mt"
clinvar_vepped_mt_path = "gs://gnomad-resources/clinvar/hail-0.2/clinvar_alleles.single.b37.vep.mt"

# GTEx files
gtex_v6_rsem_path = "gs://gnomad-berylc/tx-annotation/reheadered.GTEx_Analysis_2016-09-07_RSEMv1.2.22_transcript_tpm.txt.bgz"
gtex_v7_rsem_path = "gs://gnomad-berylc/tx-annotation/reheadered.031216.GTEx_Analysis_2016-09-07_RSEMv1.2.22_transcript_tpm.txt.bgz"
gtex_v6_tx_summary_mt_path = "gs://gnomad-berylc/tx-annotation/hail2/GTEx.V6.tx_medians.030818.mt"
gtex_v7_tx_summary_mt_path = "gs://gnomad-berylc/tx-annotation/hail2/GTEx.V7.tx_medians.030818.mt"

# Gene lists
curated_haploinsufficient_genes = "gs://gnomad-berylc/tx-annotation/HI_genes_100417.tsv"
exac_pli_genes = "gs://gnomad-berylc/tx-annotation/pLI/pLI_data.tsv.gz"
clinvar_recessive_disease_genes = "gs://gnomad-berylc/tx-annotation/ClinVar/clinvar.recessive.disease.genes.tsv"

# CSQ terms
lof_csqs = ["stop_gained","splice_donor_variant", "splice_acceptor_variant","frameshift_variant"]
missense_csqs = ["missense_variant"]
syn_csqs = ["synonymous_variant"]


def make_clinvar_hail2(clinvar_vcf_path, clinvar_variants_table, clinvar_mt_out_path):
    """
    Import ClinVar vcf file, and turn it into a usable Hail2 mt

    :param str clinvar_vcf_path: Example : "gs://gnomad-berylc/tx-annotation/hail2/clinvar_alleles_single.b37.vcf.bgz"
    :param str clinvar_variants_table: Example : "gs://gnomad-berylc/tx-annotation/hail2/clinvar_alleles_single.b37.variants_table.tsv"
    :param bool repartition:
    :param int n_partitions: Number of partitions if repartition = True
    :param str clinvar_mt_out_path: "gs://gnomad-resources/clinvar/hail-0.2/clinvar_alleles.single.b37.hail2.vepped.mt"
    :return: split and VEP'd MT
    :rtype: MatrixTable
    """
    clinvar_mt = hl.import_vcf(clinvar_vcf_path)
    variants_table = hl.import_table(clinvar_variants_table, impute=True)
    variants_table = variants_table.annotate(v = hl.parse_variant(variants_table.v))
    variants_table = (variants_table.annotate(locus=variants_table.v.locus, alleles=variants_table.v.alleles).
                      key_by('locus', 'alleles'))


    clinvar_mt = clinvar_mt.annotate_rows(va=variants_table[clinvar_mt.locus, clinvar_mt.alleles])

    clinvar_mt = split_multi_dynamic(clinvar_mt, left_aligned=False)
    clinvar_mt = clinvar_mt.repartition(100)
    clinvar_vep = hl.vep(clinvar_mt, vep_config)
    clinvar_vep.write(clinvar_mt_out_path, overwrite=True)



    t = hl.read_matrix_table(clinvar_mt_out_path)
    t.rows().show()

def make_exac_release_hail2(vcf_path, mt_out):
    """
    From Konrad, who had already did this. I didn't actually run the code.

    :param str vcf_path: Example   "gs://gnomad/raw/source/ExAC.r1.sites.vep.vcf.gz"
    :param str mt_out: Example: "gs://gnomad/raw/hail-0.2/vds/exac/exac.r1.sites.vep.vds" should be mt but whatevs
    :param bool repartition:
    :param int n_partitions: Number of partitions if repartition = True
    :return: Writes out VEP'd Hail0.2 MatrixTable
    :rtype: None
    """

    mt = hl.import_vcf(vcf_path, force_bgz=True, min_partitions=1000)
    mt = split_multi_dynamic(mt)
    mt = hl.vep(mt, vep_config)
    mt.write(mt_out)


def make_gnomad_release_hail2(vcf_path, mt_out):
    """
    Used to import, filter and VEP existing "bootleg" gnomAD VCF (01.26.2018) and write out as a Hail 0.2 MatrixTable

    :param str vcf_path:
    Example: "gs://gnomad-public/release/2.0.2/vcf/exomes/gnomad.exomes.r2.0.2.sites.vcf.bgz"
    :param str mt_out:
    Example: "gs://gnomad-berylc/tx-annotation/hail2/gnomad.exomes.r2.0.2.sites.split.vep.030818.mt"
    :return: Writes out VEP'd Hail0.2 MatrixTable
    :rtype: None
    """

    release_mt = hl.import_vcf(vcf_path, min_partitions=8000)
    release_mt = split_multi_dynamic(release_mt)

    release_mt = release_mt.annotate_rows(
        as_pass=(release_mt.info.AS_FilterStatus[release_mt.a_index - 1] == "PASS") & (
                release_mt.filters.length() == 0))

    release_mt = release_mt.filter_rows(release_mt.as_pass)
    release_mt = hl.vep(release_mt, vep_config)
    release_mt.write(mt_out)
    # Confirmed 13,443,237 variants


def get_gtex_summary(gtex_rsem_path, gtex_tx_summary_out_path, get_medians=True):
    """
    Get GTEx RSEM table with ENSTs and ENSGs as rows and GTEx samples as columns (e.g. Muscle-Skeletal.12,
    Adipose.27 etc.) and write out a table with same rows, and tissues as columns (Muscle-Skeletal, Adipose etc.)
    with cells representing summary expression of transcripts across tissues (ie. mean or median).

    :param str gtex_rsem_path: Output of RSEM quantifications from GTEx
    Example: "gs://gnomad-berylc/reheadered.GTEx_Analysis_2016-09-07_RSEMv1.2.22_transcript_tpm.txt.bgz"
    :param str gtex_tx_summary_out_path: Path to write out.
    Example: "gs://gnomad-berylc/tx-annotation/hail2/GTEx.V7.tx_medians.030818.mt"
    :param bool get_medians: Default True. If False, returns mean transcript expression per tissue
    :return: Writes out summarized GTEx transcript expression as Table.
    :rtype: None
    """

    gtex = hl.import_matrix_table(gtex_rsem_path, row_key = 'transcript_id',
                               row_fields = {'transcript_id' : hl.tstr, 'gene_id' : hl.tstr}, entry_type= hl.tfloat64)

    gtex = gtex.annotate_cols(tissue=gtex.col_id.split("\\.")[0])

    if get_medians:
        gtex = gtex.group_cols_by(gtex.tissue).aggregate(median_tx_expr=hl.median(agg.collect(gtex.x)))
    else:
        gtex = gtex.group_cols_by(gtex.tissue).aggregate(mean_tx_expr = hl.mean(agg.collect(gtex.x)))

    # Make a new column as an array of the values across tissues (per transcript)
    gtex = gtex.annotate_rows(agg_expression=agg.collect(gtex.median_tx_expr))

    # Modify the gtex table to remove version numbers
    gtex = gtex.annotate_rows(transcript_id=gtex.transcript_id.split("\\.")[0])
    gtex = gtex.annotate_rows(gene_id=gtex.gene_id.split("\\.")[0])

    gtex.write(gtex_tx_summary_out_path, overwrite = True)


def import_gene_list(gene_list_path, gene_column, ensg=False, pLI_threshold=False, peek=False):
    """
    Imports a gene list tsv and returns a set of ENSG or gene symbols

    :param str gene_list_path: Path to TSV file with gene list of interest
    :param str or None gene_column: Column in TSV file that specifies gene symbol or ENSG id.
    This column will be turned into a set.
    :param str or bool ENSG: If there are no ENSGs with version numbers in the file, specify False (Default)
    If there are ENSGs with version numbers in the file, specify column containing the ENSGs.
    :param float or bool pLI_threshold: If the file does not contain pLI scores to filter, specify False (Default)
    If the file contains pLI scores, specify threshold to filter files.
    e.g. pLI threshold = 0.95
    :param bool peek: Default False.
    If you want to peek at the gene list to get the parameters
    Print out the first few lines of the gene list tsv, returns None
    :return: Set of genes of interest
    :rtype: set or None
    """

    genes = hl.import_table(gene_list_path, impute=True)
    if peek:
        genes.show(width = 200)
        return None

    if pLI_threshold:
        genes = genes.filter(genes.pLI > pLI_threshold)

    if ensg:
        genes = genes.annotate(ensg=genes[gene_column].split("\\.")[0])
        gene_column = "ensg"

    genes = genes.aggregate(agg.collect_as_set(genes[gene_column]))

    return genes


def filter_table_to_gene_list(mt_kt, genes, gene_column_in_mt_kt):
    """Take a matrix table and return a table filtered down to a set of genes

    :param Table mt_kt:
    :param list of str or set of str genes: Genes of interest to which to filter table
    :param str gene_column_in_mt_kt: Column in matrix table that contains gene information within
    vep.transcript_consequences. often ["gene_id", "gene_symbol"]
    :return: Filtered table
    :rtype: Table
    """
    gene_names = hl.literal(genes)

    mt_kt = mt_kt.annotate(
        in_gene_of_interest=gene_names.find(lambda x: mt_kt.vep.transcript_consequences[gene_column_in_mt_kt] == x))

    mt_kt = mt_kt.filter(mt_kt.in_gene_of_interest != "NA")

    return mt_kt

def filter_clinvar_to_gene_list(mt_kt, genes, gene_column_in_mt_kt):

    gene_names = hl.literal(genes)

    mt_kt = mt_kt.annotate(
        in_gene_of_interest=gene_names.find(lambda x: mt_kt[gene_column_in_mt_kt] == x))

    mt_kt = mt_kt.filter(mt_kt.in_gene_of_interest != "NA")

    return mt_kt



def filter_table_to_csqs(mt_kt, csqs):
    """Take a matrix table and return a table filtered down to a set of CSQs

    :param Table mt_kt:
    :param list of str or set of str csqs: CSQs of interest to which to filter table
    :return: Filtered matrix table
    :rtype: Table
    """
    csqs = hl.literal(csqs)
    mt_kt = mt_kt.annotate(
        in_csq_of_interest=csqs.find(lambda x: mt_kt.worst_csq == x))

    mt_kt = mt_kt.filter(mt_kt.in_csq_of_interest != "NA")

    return mt_kt


def simplify_worst_csq(mt):
    """
    Takes in a matrix table that has passed through mt.explode_rows(mt.vep.transcript_consequences)
    Creates a new annotation 'worst_csq' which is the most deleterious of consequence terms per transcript
    e.g. if transcript_consequences.consequence_terms is ["frameshift_variant", "missense_variant", "intron_variant"]
    worst_csq is "frameshift_variant"

    :param MatrixTable mt:
    :return: Return MatrixTable with addtional worst_csq annotation
    :rtype: MatrixTable
    """
    csqs = hl.literal(CSQ_ORDER)
    mt = mt.annotate_rows(
        worst_csq=csqs.find(lambda x: mt.vep.transcript_consequences.consequence_terms.contains(x)))
    return mt


def tx_annotate_mt(mt, gtex, filter_to_genes=None, gene_column_in_mt=None, filter_to_csqs=None, out_tx_annotation_tsv=None,
                    out_tx_annotation_kt=None, filter_to_homs=False):

    """
    Annotate variants in the input MatrixTable with transcript-based expression values accross GTEx. Returns Table.

    :param MatrixTable mt:
    :param MatrixTable gtex: Input GTEx summary MatrixTable, must have transcript_id column to key by
    :param None or set filter_to_genes: Default None. If you'd like to filter the mt before annotating
    (decreases time) feed in a list or set of genes
    :param str gene_column_in_mt: Must be set if filter_to_genes != None.
    Column in matrix table that contains gene information within vep.transcript_consequences.
    often ["gene_id", "gene_symbol"]
    :param Nonr or list filter_to_csqs: Default None. If you'd like to filter the mt before annotating
    (decreases time) feed in a list or set of consequence terms.
    Example = ["stop_gained","splice_donor_variant", "splice_acceptor_variant","frameshift_variant"]
    :param None or str out_tx_annotation_tsv: Default None.
    If you'd like to write out the results table as a tsv, provide a tsv path
    :param None or str out_tx_annotation_kt: Default None.
    If you'd like to write out the results table as a Hail 0.2 table, provide a .kt path
    :param bool filter_to_homs: Default False
    If True, filter to variants with at least one homozygote in dataset
    :return: Table with columns: variant, worst_csq, ensg, LOFTEE LOF, LOFTEE LOF Flag, transcript-aware expression
    by GTEx Tissue
    :rtype: Table with variants annotated with transcript-aware tissue expression
    """

    # Create a Table copy of GTEx, key'd by transcript_id
    gtex_table = gtex.rows().key_by("transcript_id")

    # Explode the mt for the transcript consequences to be able to key by transcript ID
    mt = mt.explode_rows(mt.vep.transcript_consequences)

    # Add worst csq to the mt
    mt = simplify_worst_csq(mt)

    mt_kt = mt.rows()

    if filter_to_genes:
        print("Filtering to genes of interest")
        mt_kt = filter_table_to_gene_list(mt_kt, filter_to_genes, gene_column_in_mt)

    if filter_to_csqs:
        print("Filtering to csqs in %s" %(",".join(filter_to_csqs)))
        mt_kt = filter_table_to_csqs(mt_kt, filter_to_csqs)

    if filter_to_homs:
        print("Filtering to variants with at least 1 homozygote sample in dataset")
        mt_kt = mt_kt.filter(mt_kt.info.Hom[mt_kt.a_index - 1] > 0)


    # Annotate mt with the gtex values (ie. join them)
    mt_kt = mt_kt.annotate(tx_data=gtex_table[mt_kt.vep.transcript_consequences.transcript_id])

    # Group by gene, worst_csq and variant, and do a pairwise-sum
    grouped_table = (
        mt_kt.group_by(worst_csq=mt_kt.worst_csq,
                       ensg=mt_kt.vep.transcript_consequences.gene_id,
                       locus=mt_kt.locus,
                       alleles = mt_kt.alleles,
                       lof=mt_kt.vep.transcript_consequences.lof,
                       lof_flag=mt_kt.vep.transcript_consequences.lof_flags).aggregate(tx_annotation=agg.array_sum(mt_kt.tx_data.agg_expression)))

    # Expand the columns from the arrays and add tissues as headers
    tissue_ids = gtex.aggregate_cols(agg.collect(gtex.tissue))
    d = {tiss: i for i, tiss in enumerate(tissue_ids)}

    # This is currently a hack but Tim will fix it at which point I can remove "replace"s
    tx_annotation_table = grouped_table.annotate(
        **{tissue_id.replace("-", "_").replace(" ", "_").replace("(", "_").replace(")", "_"):
               grouped_table.tx_annotation[d[tissue_id]] for tissue_id in tissue_ids})

    if out_tx_annotation_tsv:
        print("Writing tsv file to %s" %out_tx_annotation_tsv)
        tx_annotation_table.export(out_tx_annotation_tsv)

    if out_tx_annotation_kt:
        print("Writing Table to %s" % out_tx_annotation_kt)
        tx_annotation_table.write(out_tx_annotation_kt)

    return tx_annotation_table



def read_tx_annotation_tables(mt_path, gtex_tx_summary_path):
    # Read in MT and filter to PASS variants
    mt = hl.read_matrix_table(mt_path)

    # Read in GTEx summary, turn into a table
    gtex = hl.read_matrix_table(gtex_tx_summary_path)

    return mt, gtex


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack('@berylc', main, args)
    else:
        main(args)



# #This is the workflow for gnomad

mt, gtex = read_tx_annotation_tables(gnomad_vepped_mt_path, gtex_v7_tx_summary_mt_path)

# #gnomAD HI genes
genes_of_interest = import_gene_list(curated_haploinsufficient_genes, gene_column = "ENSGID", ensg=True)

tx_table_gnomad_hi = tx_annotate_mt(mt, gtex, filter_to_genes=genes_of_interest, gene_column_in_mt = "gene_id",
                                    out_tx_annotation_tsv = "gs://gnomad-berylc/tx-annotation/hail2/gene_list_comparisons/gtex_v7/" \
                                    "v7.gnomad.exomes.r2.0.2.tx_annotated.higenes.030818.tsv.bgz")


# #gnomAD autosomal recessive genes
recessive_genes = import_gene_list(clinvar_recessive_disease_genes, gene_column="symbol")

tx_table_gnomad_recessive = tx_annotate_mt(mt, gtex, filter_to_genes=recessive_genes, gene_column_in_mt = "gene_symbol",
                                            filter_to_homs = True, out_tx_annotation_tsv =
                                    "gs://gnomad-berylc/tx-annotation/hail2/gene_list_comparisons/gtex_v7/"
                                    "v7.gnomad.exomes.r2.0.2.tx_annotated.recessive_genes.030818.tsv.bgz")


# #gnomAD pLI
pli_genes = import_gene_list(exac_pli_genes, gene_column="gene", pLI_threshold=0.95)

tx_table_gnomad_pli = tx_annotate_mt(mt, gtex, filter_to_genes=pli_genes, gene_column_in_mt="gene_symbol",
                                      filter_to_csqs = lof_csqs, out_tx_annotation_tsv=
                                      "gs://gnomad-berylc/tx-annotation/hail2/gene_list_comparisons/gtex_v7/"
                                      "v7.gnomad.exomes.r2.0.2.tx_annotated.pli_genes.030818.tsv.bgz")

tx_table_gnomad_pli_syn = tx_annotate_mt(mt, gtex, filter_to_genes=pli_genes, gene_column_in_mt="gene_symbol",
                                          filter_to_csqs=syn_csqs, out_tx_annotation_tsv=
                                          "gs://gnomad-berylc/tx-annotation/hail2/gene_list_comparisons/gtex_v7"
                                          "v7.gnomad.exomes.r2.0.2.tx_annotated.pli_genes_synonymous.030818.tsv.bgz")


#ClinVar
clinvar_mt, gtex = read_tx_annotation_tables(clinvar_vepped_mt_path, gtex_v7_tx_summary_mt_path)
tx_table_clinvar = tx_annotate_mt(clinvar_mt, gtex)


# #Now joining to add other information into the file
clinvar_kt = clinvar_mt.rows()
clinvar_kt = clinvar_kt.select(clinvar_kt.locus, clinvar_kt.alleles, clinvar_kt.va.symbol, clinvar_kt.va.clinical_significance,
                               clinvar_kt.va.pathogenic, clinvar_kt.va.conflicted, clinvar_kt.va.gold_stars)

tx_table_clinvar = tx_table_clinvar.key_by('locus', 'alleles')

tx_table_clinvar = tx_table_clinvar.join(clinvar_kt)

tx_table_clinvar.export("gs://gnomad-berylc/tx-annotation/hail2/gene_list_comparisons/"
                          "clinvar.alleles.single.b37.tx_annotated.allgenes.030818.tsv.bgz")

tx_table_clinvar.write("gs://gnomad-berylc/tx-annotation/hail2/gene_list_comparisons/gtex_v7/"
                       "v7.clinvar.alleles.single.b37.tx_annotated.allgenes.030818.kt")

# #Now filtering clinvar to gene lists for quick analysis in R
clinvar_full = hl.read_table("gs://gnomad-berylc/tx-annotation/hail2/gene_list_comparisons/gtex_v7/"
                             "v7.clinvar.alleles.single.b37.tx_annotated.allgenes.030818.kt")

recessive_genes = import_gene_list(clinvar_recessive_disease_genes, gene_column="symbol")

clinvar_rec = filter_clinvar_to_gene_list(clinvar_full, genes=recessive_genes, gene_column_in_mt_kt="symbol")

clinvar_rec.export("gs://gnomad-berylc/tx-annotation/hail2/gene_list_comparisons/gtex_v7/"
                   "v7.clinvar.alleles.single.b37.tx_annotated.recessive_genes.030818.tsv.bgz")

hi_genes = import_gene_list(curated_haploinsufficient_genes, gene_column="ENSGID", ensg=True)

clinvar_hi = filter_clinvar_to_gene_list(clinvar_full, genes=hi_genes, gene_column_in_mt_kt="ensg")

clinvar_hi.export("gs://gnomad-berylc/tx-annotation/hail2/gene_list_comparisons/gtex_v7/"
                   "v7.clinvar.alleles.single.b37.tx_annotated.hi_genes.030818.tsv.bgz")
