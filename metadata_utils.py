import re

from collections import defaultdict
import gspread
import os
import pandas as pd
import subprocess
import yaml
from pyairtable import Api

from google.oauth2.service_account import Credentials

VALUE_RENDER_OPTION__FORMATTED_VALUE = "FORMATTED_VALUE"
VALUE_RENDER_OPTION__UNFORMATTED_VALUE = "UNFORMATTED_VALUE"
VALUE_RENDER_OPTION__FORMULA = "FORMULA"

_GSPREAD_CLIENT = None

# Airtable IDs.
RNA_SEQ_BASE_ID = "app034nYTJwpo9xxw"
DATA_PATHS_TABLE_ID = "tblmT2O7VQasTKSDx"
DATA_PATHS_VIEW_ID = "viwUzuAxagca7zF0l"

METADATA_TABLE_ID = "tblcluIOjz23rZDya"
METADATA_VIEW_ID = "viwTjCsk3SccsOJIz"

DOWNSTREAM_ANALYSIS_TABLE_ID = "tblyJvv0ozZjXIutA"
DOWNSTREAM_ANALYSIS_VIEW_ID = "viwkCT48ws9yDTDmf"

# Parse YAML config file
config_file_path = os.path.expanduser("~/.tgg_rnaseq_pipelines")
if not os.path.isfile(config_file_path):
    raise Exception(f"Config file not found: {config_file_path}")

with open(config_file_path, "r") as f:
    config = yaml.safe_load(f)

credentials = ["service-account", "airtable-token"]
for credential in credentials:
    if credential not in config:
        raise Exception(
            f"Config file {config_file_path} doesn't contain a \"{credential}\" key")

SERVICE_ACCOUNT_CREDENTIALS_JSON_PATH = os.path.expanduser(
    str(config["service-account"]))
AIRTABLE_TOKEN = str(config["airtable-token"])


# %%
def normalize_rnaseq_sample_id(sample_id):
    sample_id = sample_id.strip()
    sample_id = sample_id.replace(".", "-")
    return sample_id


def get_spreadsheet(spreadsheet_name):
    global _GSPREAD_CLIENT
    if _GSPREAD_CLIENT is None:
        creds = Credentials.from_service_account_file(
            SERVICE_ACCOUNT_CREDENTIALS_JSON_PATH,
            scopes=[
                'https://www.googleapis.com/auth/spreadsheets',
                'https://www.googleapis.com/auth/drive.file',
                'https://www.googleapis.com/auth/drive',
            ]
        )

        _GSPREAD_CLIENT = gspread.authorize(creds)

    spreadsheet = _GSPREAD_CLIENT.open(spreadsheet_name)

    return spreadsheet


## Spreadsheet must be Shared with 733952080251-compute@developer.gserviceaccount.com
_RNASEQ_DOWNSTREAM_ANALYSIS_METADATA_WORKSHEET = None
_RNASEQ_METADATA_SPREADSHEET = None
_RNASEQ_METADATA_WORKSHEET = None
_DATA_PATHS_WORKSHEET = None
_IMPUTED_METADATA_WORKSHEET = None
_BERYLS_WORKSHEET = None
_BERYLS_WORKSHEET_2 = None
_BERYLS_WORKSHEET_3 = None

_GTEX_METADATA_SPREADSHEET = None
_GTEX_RNASEQ_SAMPLE_METADATA_WORKSHEET = None
_GTEX_WES_SAMPLE_METADATA_WORKSHEET = None
_GTEX_WGS_SAMPLE_METADATA_WORKSHEET = None
_GTEX_INDIVIDUAL_METADATA_WORKSHEET = None


def get_table(base_id, table_id, view_id):
    api = Api(AIRTABLE_TOKEN)
    all_tables = api.table(base_id, table_id)
    cur_table = all_tables.all(view=view_id)
    return cur_table


def read_from_airtable(base_id, table_id, view_id):
    cur_table = get_table(base_id, table_id, view_id)
    for i in range(len(cur_table)):
        record = cur_table[i]
        cur_table[i] = record["fields"]

    df = pd.DataFrame(cur_table)
    return df


def map_airtable_and_sample_id(base_id, table_id, view_id):
    cur_table = get_table(base_id, table_id, view_id)
    map_dict = {}
    for i in range(len(cur_table)):
        record = cur_table[i]
        record_id = record["id"]
        sample_id = record["fields"]["sample_id"]
        map_dict[sample_id] = record_id
    return map_dict


def get_missing_records(base_id, table_id, view_id, column):
    cur_table = get_table(base_id, table_id, view_id)
    records_with_na_in_column = []
    ids_with_na_in_column = []
    for i in range(len(cur_table)):
        record = cur_table[i]
        if column and column not in record["fields"].keys():
            records_with_na_in_column.append(record["fields"])
            ids_with_na_in_column.append(record["id"])

    records_with_na_in_column_df = pd.DataFrame(records_with_na_in_column)
    return records_with_na_in_column_df, ids_with_na_in_column


def write_to_airtable(base_id, table_id, new_records, to_replace=False,
                      fields_to_match=["sample_id"]):
    api = Api(AIRTABLE_TOKEN)
    all_tables = api.table(base_id, table_id)
    all_tables.batch_upsert(new_records,
                            key_fields=fields_to_match,
                            replace=to_replace)


def inner_join_two_tables_into_df(base_id1, table_id1, view_id1,
                                  base_id2, table_id2, view_id2):
    df1 = read_from_airtable(base_id1,
                             table_id1,
                             view_id1)
    df2 = read_from_airtable(base_id2,
                             table_id2,
                             view_id2)
    df_res = pd.merge(df1, df2,
                      how='inner',
                      left_index=True, right_index=True,
                      suffixes=('', '_drop'))
    df_res.drop(
        [col for col in df_res.columns if 'drop' in col],
        axis=1,
        inplace=True)
    return df_res


def get_gtex_v8_metadata_spreadsheet():
    global _GTEX_METADATA_SPREADSHEET
    _GTEX_METADATA_SPREADSHEET = get_spreadsheet("GTEx v8 metadata")
    return _GTEX_METADATA_SPREADSHEET


def get_gtex_rnaseq_sample_metadata_worksheet():
    global _GTEX_RNASEQ_SAMPLE_METADATA_WORKSHEET
    _GTEX_RNASEQ_SAMPLE_METADATA_WORKSHEET = get_gtex_v8_metadata_spreadsheet().worksheet(
        "RNA-seq sample metadata (auto)")
    return _GTEX_RNASEQ_SAMPLE_METADATA_WORKSHEET


def get_gtex_wes_sample_metadata_worksheet():
    global _GTEX_WES_SAMPLE_METADATA_WORKSHEET
    _GTEX_WES_SAMPLE_METADATA_WORKSHEET = get_gtex_v8_metadata_spreadsheet().worksheet(
        "WES sample metadata (auto)")
    return _GTEX_WES_SAMPLE_METADATA_WORKSHEET


def get_gtex_wgs_sample_metadata_worksheet():
    global _GTEX_WGS_SAMPLE_METADATA_WORKSHEET
    _GTEX_WGS_SAMPLE_METADATA_WORKSHEET = get_gtex_v8_metadata_spreadsheet().worksheet(
        "WGS sample metadata (auto)")
    return _GTEX_WGS_SAMPLE_METADATA_WORKSHEET


def get_gtex_individual_metadata_worksheet():
    global _GTEX_INDIVIDUAL_METADATA_WORKSHEET
    _GTEX_INDIVIDUAL_METADATA_WORKSHEET = get_gtex_v8_metadata_spreadsheet().worksheet(
        "individual metadata (auto)")
    return _GTEX_INDIVIDUAL_METADATA_WORKSHEET


def get_rnaseq_metadata_spreadsheet():
    global _RNASEQ_METADATA_SPREADSHEET
    _RNASEQ_METADATA_SPREADSHEET = get_spreadsheet("RNA-seq metadata")
    return _RNASEQ_METADATA_SPREADSHEET


def get_rnaseq_metadata_worksheet():
    global _RNASEQ_METADATA_WORKSHEET
    _RNASEQ_METADATA_WORKSHEET = get_rnaseq_metadata_spreadsheet().worksheet(
        "sample_metadata")
    return _RNASEQ_METADATA_WORKSHEET


def get_data_paths_worksheet():
    global _DATA_PATHS_WORKSHEET
    _DATA_PATHS_WORKSHEET = get_rnaseq_metadata_spreadsheet().worksheet("data_paths")
    return _DATA_PATHS_WORKSHEET


def get_rnaseq_downstream_analysis_metadata_worksheet():
    global _RNASEQ_DOWNSTREAM_ANALYSIS_METADATA_WORKSHEET
    _RNASEQ_DOWNSTREAM_ANALYSIS_METADATA_WORKSHEET = get_rnaseq_metadata_spreadsheet().worksheet(
        "downstream_analysis")
    return _RNASEQ_DOWNSTREAM_ANALYSIS_METADATA_WORKSHEET


def get_beryls_supplementary_table_worksheet():
    global _BERYLS_WORKSHEET
    _BERYLS_WORKSHEET = get_rnaseq_metadata_spreadsheet().worksheet(
        "Beryl's Supplementary Table 1")
    return _BERYLS_WORKSHEET


def get_beryls_rnaseq_probands_worksheet():
    global _BERYLS_WORKSHEET_2
    _BERYLS_WORKSHEET_2 = get_rnaseq_metadata_spreadsheet().worksheet(
        "Copy of Beryl's RNAseq Probands")
    return _BERYLS_WORKSHEET_2


def get_beryls_seqr_data_worksheet():
    global _BERYLS_WORKSHEET_3
    _BERYLS_WORKSHEET_3 = get_rnaseq_metadata_spreadsheet().worksheet(
        "Copy of Beryl's Seqr-data")
    return _BERYLS_WORKSHEET_3


def get_rnaseq_metadata_df(value_render_option=VALUE_RENDER_OPTION__FORMATTED_VALUE):
    rows = get_rnaseq_metadata_worksheet().get(value_render_option=value_render_option)
    return pd.DataFrame(data=rows[1:], columns=rows[0])


def get_rnaseq_downstream_analysis_metadata_df(
        value_render_option=VALUE_RENDER_OPTION__FORMATTED_VALUE):
    rows = get_rnaseq_downstream_analysis_metadata_worksheet().get(
        value_render_option=value_render_option)
    return pd.DataFrame(data=rows[1:], columns=rows[0])


def get_rnaseq_data_paths_df(value_render_option=VALUE_RENDER_OPTION__FORMATTED_VALUE):
    rows = get_data_paths_worksheet().get(value_render_option=value_render_option)
    return pd.DataFrame(data=rows[1:], columns=rows[0])


def get_beryls_supplementary_table_df(
        value_render_option=VALUE_RENDER_OPTION__FORMATTED_VALUE):
    rows = get_beryls_supplementary_table_worksheet().get(
        value_render_option=value_render_option)
    return pd.DataFrame(data=rows[1:], columns=rows[0])


def get_beryls_rnaseq_probands_df(
        value_render_option=VALUE_RENDER_OPTION__FORMATTED_VALUE):
    rows = get_beryls_rnaseq_probands_worksheet().get(
        value_render_option=value_render_option)
    return pd.DataFrame(data=rows[1:], columns=rows[0])


def get_beryls_seqr_data_df(value_render_option=VALUE_RENDER_OPTION__FORMATTED_VALUE):
    rows = get_beryls_seqr_data_worksheet().get(value_render_option=value_render_option)
    return pd.DataFrame(data=rows[1:], columns=rows[0])


def get_rnaseq_metadata_joined_with_paths_df(
        value_render_option=VALUE_RENDER_OPTION__FORMATTED_VALUE):
    df1 = get_rnaseq_data_paths_df(value_render_option=value_render_option)
    df2 = get_rnaseq_metadata_df(value_render_option=value_render_option)
    df2 = df2[[c for c in df2.columns if c not in ("star_pipeline_batch",
                                                   "batch_date_from_hg19_bam_header")]]  # remove columns that exist in both tables
    return df1.merge(df2, on="sample_id", how="left").set_index("sample_id", drop=False)


def get_rnaseqc_metrics(rnaseqc_metrics_file_path):
    output = subprocess.check_output("gsutil cat %s" % rnaseqc_metrics_file_path,
                                     shell=True, encoding="UTF-8")
    metrics_dict = {}
    for i, line in enumerate(output.rstrip().split("\n")):
        key, value = line.split("\t")
        metrics_dict[key] = value

    return metrics_dict


def get_gtex_rnaseq_sample_metadata_df(
        value_render_option=VALUE_RENDER_OPTION__FORMATTED_VALUE):
    rows = get_gtex_rnaseq_sample_metadata_worksheet().get(
        value_render_option=value_render_option)
    return pd.DataFrame(data=rows[1:], columns=rows[0]).set_index("SAMPID", drop=False)


def get_gtex_wes_sample_metadata_df(
        value_render_option=VALUE_RENDER_OPTION__FORMATTED_VALUE):
    rows = get_gtex_wes_sample_metadata_worksheet().get(
        value_render_option=value_render_option)
    return pd.DataFrame(data=rows[1:], columns=rows[0]).set_index("SAMPID", drop=False)


def get_gtex_wgs_sample_metadata_df(
        value_render_option=VALUE_RENDER_OPTION__FORMATTED_VALUE):
    rows = get_gtex_wgs_sample_metadata_worksheet().get(
        value_render_option=value_render_option)
    return pd.DataFrame(data=rows[1:], columns=rows[0]).set_index("SAMPID", drop=False)


def get_analysis_batches():
    df = get_rnaseq_downstream_analysis_metadata_df(
        value_render_option=VALUE_RENDER_OPTION__FORMATTED_VALUE)
    analysis_batch_to_tissue = defaultdict(set)
    analysis_batch_to_sex = defaultdict(set)

    analysis_batches = {}
    for _, r in df.iterrows():
        analysis_batch = r["tissue"]
        if not analysis_batch:
            continue
        analysis_batch = analysis_batch.strip()
        if not analysis_batch or analysis_batch == "x":
            continue

        analysis_batch_to_tissue[analysis_batch].add(r["tissue"])
        analysis_batch_to_sex[analysis_batch].add(r["sex"])

    for analysis_batch, tissue in analysis_batch_to_tissue.items():
        if len(tissue) != 1:
            raise ValueError(f"Expected 1 tissue for {analysis_batch}. Found: {tissue}")
        tissue = next(iter(tissue))
        sex = analysis_batch_to_tissue[analysis_batch]
        if len(sex) > 1:
            sex = "both"
        else:
            sex = next(iter(sex))

        analysis_batches[analysis_batch] = {
            "tissue": tissue,
            "sex": sex,
            "samples": list(df[df["tissue"] == analysis_batch].sample_id)
        }

    # TODO fix empty values in spreadsheet "analysis batch" column
    return analysis_batches