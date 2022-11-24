from Bio import SeqIO
import pandas as pd
import glob, os
from tqdm import tqdm

from ont_fast5_api.fast5_file import Fast5FileTypeError
from ont_fast5_api.fast5_interface import check_file_type, MULTI_READ
from ont_fast5_api.multi_fast5 import MultiFast5File
from ont_fast5_api.conversion_tools.multi_to_single_fast5 import create_single_f5

def get_barcodes(directory):
    df_all = None
    fastQ_files = glob.glob(directory + '/pass/**/*.fastq', recursive=True)
    for fastQ_file in tqdm(fastQ_files, desc='Getting barcodes'):
        try:
            barcode = [i for i in fastQ_file.split('/') if 'barcode' in i][0]
            records = []
            for record in SeqIO.parse(fastQ_file, "fastq"):
                records.append(record.id)
            df = pd.DataFrame([barcode] * len(records), index=records, columns=['barcode'])

            if df_all is None:
                df_all = df
            else:
                df_all = pd.concat((df_all, df))

        except IndexError:
            pass
    return (df_all)


def try_multi_to_single_conversion_split_barcodes(directory, df):
    fast5_files = glob.glob(directory + '/workspace/*.fast5', recursive=False)

    output_files = []
    not_found = 0

    for fast5_file in tqdm(fast5_files, desc='Splitting fast5 files'):
        with MultiFast5File(fast5_file, 'r') as multi_f5:
            file_type = check_file_type(multi_f5)
            if file_type != MULTI_READ:
                raise Fast5FileTypeError("Could not convert Multi->Single for file type '{}' with path '{}'"
                                         "".format(file_type, fast5_file))
            for read in multi_f5.get_reads():
                try:
                    subfolder = r'workspace/' + df.at[read.read_id, 'barcode']
                    output_file = os.path.join(directory, subfolder, "{}.fast5".format(read.read_id))
                    create_single_f5(output_file, read)
                    output_files.append(os.path.basename(output_file))
                except Exception as e:
                    not_found += 1
    print(f'Converted {len(output_files)} reads; Skipped {not_found} reads.')
    return output_files

if __name__ == '__main__':
    directory = r'/home/kuijntjes/Desktop/2022-11-02_WholeCellExtractGalLocusCindy2_DEPC/no_sample/20221102_1621_MN30914_ajg585_e7b9be34/fast5_skip_test/guppy'
    df = get_barcodes(directory)
    tmp = try_multi_to_single_conversion_split_barcodes(directory, df)
