import hailtop.batch as hb
import argparse
import numpy as np
import hailtop.fs as hfs
from fraser2_rscirpts import count_split_reads_single_sample_r, \
    count_reads_all_samples_r, run_fraser_r
from metadata_utils import read_from_airtable, \
    write_to_airtable, RNA_SEQ_BASE_ID, \
    DATA_PATHS_TABLE_ID, \
    DATA_PATHS_VIEW_ID, get_gtex_metadata, switch_to_gmail_account

REGION = ["us-central1"]
# DEFAULT_BATCH_NAME = "all"
# SAVE_PREFIX = "with_volcano"
DEFAULT_CPU = 2 ** 4
DEFAULT_MEMORY = "highmem"
PER_BAM_SIZE = 50
DEFAULT_STORAGE = f"{PER_BAM_SIZE}G"
ANNOTATION_DELIMITER = "\t"  # delimiter of the annotation input file
SAMPLE_ID_COL = "sample_id"
BAM_PATH_COL = "star_bam"

DOCKER_IMAGE = "gcr.io/cmg-analysis/fraser2@sha256:38e7e777a08886b5d4789b4c06f5433af60953542774134165281b7c39d35eeb"
GENE_MODELS_GFF = "gs://tgg-rnaseq/ref/MANE.GRCh38.v1.0.ensembl_genomic.without_ensg_versions.gff.gz"


def get_ids_and_bam_paths():
    data_paths_table = read_from_airtable(RNA_SEQ_BASE_ID,
                                          DATA_PATHS_TABLE_ID,
                                          DATA_PATHS_VIEW_ID)
    tissue_table = data_paths_table.loc[data_paths_table["imputed_tissue"] ==
                                        args.tissue]
    tissue_table = tissue_table[~(tissue_table["exclude"] == "yes")]
    sample_ids = np.array(tissue_table["sample_id"])
    bam_paths = np.array(tissue_table["star_bam"])
    return sample_ids, bam_paths


# def create_symbolic_links(batch, job, path, link_path):
#     localized_path = batch.read_input(path)
#     job.command(f"ln -s {localized_path} {link_path}")

def create_symbolic_links(batch, job, path, link_path, is_gtex):
    if is_gtex:
        job.command(
            f"gsutil -u {args.requester_pays_project} -m cp {path} {link_path}")
    else:
        localized_path = batch.read_input(path)
        job.command(f"ln -s {localized_path} {link_path}")


def load_bam_and_index_file(batch, cur_job, sample_id, bam_path, is_gtex):
    bam_index_path = f"{bam_path}.bai"
    link_bam_path = f"{sample_id}.bam"
    link_bam_index_path = f"{sample_id}.bam.bai"
    create_symbolic_links(batch, cur_job, bam_path, link_bam_path, is_gtex)
    create_symbolic_links(batch, cur_job, bam_index_path, link_bam_index_path, is_gtex)


def count_and_cache_split_reads(batch, sample_id, bam_path):
    # count and cache split reads for one sample
    cur_count_job = batch.new_job(f"count_{sample_id}")
    bam_size = hfs.ls(bam_path)[0].size
    # print(bam_size)
    cur_count_job.storage(bam_size + 1e11)  # set the job storage to be the file size
    cur_count_job.command("cd /io")  # enter the dir where storage is mounted

    # switch email account for GTEx samples
    if args.with_gtex and "GTEX" in sample_id:
        switch_to_gmail_account(cur_count_job)
        # create local bam file paths for R
        load_bam_and_index_file(batch, cur_count_job, sample_id, bam_path, True)
    else:
        load_bam_and_index_file(batch, cur_count_job, sample_id, bam_path, False)

    # count split reads
    cur_count_job.command(f"""xvfb-run Rscript -e '
    {count_split_reads_single_sample_r(DEFAULT_CPU, sample_id, sample_id + ".bam")}
    '
    """)

    # save the split read counts in cache folder
    cached_filename = f"count_split_reads_{sample_id}.tar.gz"
    cur_count_job.command("ls -lh .")
    cur_count_job.command(f"tar czf {cached_filename} cache")
    cur_count_job.command(f"cp {cached_filename} {cur_count_job.ofile}")
    batch.write_output(cur_count_job.ofile, f"{fraser_dir}/{cached_filename}")

    return cur_count_job


def read_from_cloud(path):
    dat = []
    with hfs.open(path) as f:
        for line in f:
            dat.append(line.strip())
    return dat


def save_iter_to_cloud(job, iterable_obj):
    for item in iterable_obj:
        job.command(f"echo {item} >> {job.ofile}")


def save_file_to_cloud(batch, cur_job, file, id):
    cur_job.command(f"cp {file} {cur_job.ofile}")
    batch.write_output(cur_job.ofile, f"{fraser_dir}/{id}_{file}")


def copy_split_read_counts_files(batch, job, sample_ids):
    for cur_id in sample_ids:
        path = f"{fraser_dir}/count_split_reads_{cur_id}.tar.gz"  # cloud path
        link_path = f"count_split_reads_{cur_id}.tar.gz"  # soft link path
        create_symbolic_links(batch, job, path, link_path, False)


def get_split_reads(batch, sample_ids, bam_paths):
    count_jobs = []
    print(len(sample_ids))

    # if the split read counts of a sample are not cached, save the sample id
    for i in range(len(sample_ids)):
        cur_id = sample_ids[i]
        cur_bam_path = bam_paths[i]
        cur_saved_bam_path = f"{fraser_dir}/count_split_reads_{cur_id}.tar.gz"
        # if sample_id_cached_set is None or cur_id not in sample_id_cached_set or not hfs.is_file(cur_saved_bam_path):  # count split read counts for new samples
        if not hfs.is_file(cur_saved_bam_path):
            print(cur_id)
            # print(cur_saved_bam_path)
            cur_count_job = count_and_cache_split_reads(batch, cur_id, cur_bam_path)
            count_jobs.append(cur_count_job)

    return count_jobs


def get_all_reads(batch, cur_job, sample_ids, bam_paths):
    saved_fds_path = f"{fraser_dir}/{args.job_name}_savedObjects"
    if hfs.is_file(saved_fds_path):  # if the fds exists
        return None

    # cur_job = batch.new_job(f"get_all_reads_{type}")
    cur_job.storage(f"{15 * to_use_ids.shape[0] + 100}G")
    if args.with_gtex:
        switch_to_gmail_account(cur_job)
        load_split_reads_and_bam_files(batch, cur_job, sample_ids, bam_paths, True)
    else:
        load_split_reads_and_bam_files(batch, cur_job, sample_ids, bam_paths, False)
    env_var_id, env_var_path = get_env_vars(sample_ids)

    print(env_var_id)
    print(env_var_path)
    cur_job.command(f"""xvfb-run Rscript -e '
            {count_reads_all_samples_r(DEFAULT_CPU)}
            ' {env_var_id} {env_var_path}
            """)

    cur_job.command("ls -lh .")
    cur_job.command(f"tar czf savedObjects.tar.gz savedObjects")
    cur_job.command(f"cp savedObjects.tar.gz {cur_job.ofile}")
    batch.write_output(cur_job.ofile, saved_fds_path)

    return cur_job


def load_split_reads_and_bam_files(batch, cur_job, sample_ids, bam_paths, is_gtex):
    cur_job.command("cd /io")

    # copy all split read counts files to the current directory
    copy_split_read_counts_files(batch, cur_job, sample_ids)
    cur_job.command("ls -lh")
    # decompress .tar.gz files and rebuild the cache folder
    cur_job.command("for i in count_split_reads*.tar.gz; do tar xzf $i; done")
    cur_job.command("ls -lh cache")

    # load bam and bam index file to the current directory
    for i in range(len(sample_ids)):
        cur_id = sample_ids[i]
        cur_bam_path = bam_paths[i]
        if is_gtex and "GTEX" in cur_id:
            load_bam_and_index_file(batch, cur_job, cur_id, cur_bam_path, True)
        else:
            load_bam_and_index_file(batch, cur_job, cur_id, cur_bam_path, False)


def get_env_vars(sample_ids):
    env_var_id = ":".join(sample_ids)
    link_bam_paths = [f"{cur_id}.bam" for cur_id in sample_ids]
    env_var_path = ":".join(link_bam_paths)
    return env_var_id, env_var_path


def run_fraser(batch, cur_job, sample_ids, bam_paths, cur_type):
    if args.with_gtex:
        switch_to_gmail_account(cur_job)

    # localize saved fds to the container
    saved_fds_path = f"{fraser_dir}/{args.job_name}_savedObjects"
    create_symbolic_links(batch, cur_job, saved_fds_path, "savedObjects", False)
    cur_job.command("ls -lh")
    cur_job.command(f"tar xzf savedObjects")

    gene_models_gff_path = "MANE.GRCh38.v1.0.ensembl_genomic.without_ensg_versions.gff.gz"  # link path
    create_symbolic_links(batch, cur_job, GENE_MODELS_GFF, gene_models_gff_path, False)

    # run fraser2
    env_var_id, env_var_path = get_env_vars(sample_ids)
    result_table = f"{cur_type}_results.txt"
    heatmap_before_ae = f"{cur_type}_before_ae_heatmap.png"
    heatmap_after_ae = f"{cur_type}_after_ae_heatmap.png"
    enc_dim_auc = f"{cur_type}_enc_dim_auc.png"
    enc_dim_loss = f"{cur_type}_enc_dim_loss.png"
    aberrant = f"{cur_type}_aberrant.png"

    delta_psi_threshold = 0.1
    padj_threshold = 0.3
    min_reads = 2

    zip_dat = f"{args.job_name}_{cur_type}_zip_dat.tar.gz"

    cur_job.command(f"""xvfb-run Rscript -e '
    {run_fraser_r(cur_type, DEFAULT_CPU, result_table, heatmap_before_ae, heatmap_after_ae, enc_dim_auc, enc_dim_loss,
                  aberrant, delta_psi_threshold, padj_threshold, min_reads, gene_models_gff_path)}
    ' {env_var_id} {env_var_path}
    """)
    cur_job.command("ls -lh .")

    # save results to cloud path
    # cur_job.command(
    #     f"tar czf {zip_dat}.tar.gz {result_table} {heatmap_before_ae} {heatmap_after_ae} "
    #     f"{enc_dim_auc} {enc_dim_loss}")
    cur_job.command(f"tar czf {zip_dat}.tar.gz {result_table} filtered_{result_table} "
                    f"volcano_*")
    cur_job.command(f"cp {zip_dat}.tar.gz {cur_job.ofile}")
    batch.write_output(cur_job.ofile, f"{fraser_dir}/{zip_dat}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--billing-project", type=str, help="Project to bill under.",
                        default="tgg-rare-disease")
    parser.add_argument("--requester-pays-project", type=str,
                        help="Requester pays project to bill under.",
                        default="cmg-analysis")
    parser.add_argument("--file-dir", type=str,
                        help="The directory to store results table.",
                        default="gs://jialan-storage")
    parser.add_argument("--tissue", type=str,
                        help="Tissue type of the samples to be processed. Should be "
                             "one of [\"muscle\", \"fibroblasts\", \"lymphocytes\","
                             "\"whole_blood\"].",
                        required=True)
    parser.add_argument("--job-name", type=str,
                        required=True)
    parser.add_argument("--with-gtex", action="store_true",
                        help="Add the top 100 GTEx samples to RDG samples.")
    parser.add_argument("-s", "--sample-ids", type=str,
                        help="A text file containing IDs of samples with the "
                             "specified tissue type to be processed. "
                             "If none provided, all samples of the specified tissue "
                             "type in the airtable will be processed.")
    args = parser.parse_args()

    backend = hb.ServiceBackend(billing_project=args.billing_project,
                                remote_tmpdir=args.file_dir,
                                regions=REGION)

    batch = hb.Batch(backend=backend, name=args.job_name,
                     requester_pays_project=args.requester_pays_project,
                     default_image=DOCKER_IMAGE,
                     default_cpu=DEFAULT_CPU,
                     default_memory=DEFAULT_MEMORY,
                     default_storage=DEFAULT_STORAGE)
    fraser_dir = f"{args.file_dir}/fraser2"

    to_use_ids, to_use_bam_paths = get_ids_and_bam_paths()

    if args.with_gtex:
        gtex_table = get_gtex_metadata()
        gtex_for_tissue = gtex_table[gtex_table["tissue"] == args.tissue]

        gtex_ids = np.array(gtex_for_tissue["sample_id"])
        gtex_bam_paths = np.array(gtex_for_tissue["bam_path"])
        to_use_ids = np.concatenate([to_use_ids, gtex_ids], axis=0)
        to_use_bam_paths = np.concatenate([to_use_bam_paths, gtex_bam_paths], axis=0)

    psi_types = ["jaccard"]
    print("The number of samples used in this batch run is: ", len(to_use_ids))

    count_jobs = get_split_reads(batch, to_use_ids, to_use_bam_paths)

    all_reads_job = batch.new_job(f"get_all_reads")
    if len(count_jobs) > 0:
        all_reads_job.depends_on(*count_jobs)
    get_all_reads(batch, all_reads_job, to_use_ids, to_use_bam_paths)

    for cur_type in psi_types:
        fraser_job = batch.new_job(f"fraser2_{cur_type}")
        fraser_job.depends_on(all_reads_job)
        fraser_job.storage(f"100G")
        run_fraser(batch, fraser_job, to_use_ids, to_use_bam_paths, cur_type)

    batch.run()