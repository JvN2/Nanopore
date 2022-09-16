from Bio import SeqIO
import numpy as np

def phred_to_quality(phred):
    return np.asarray([ord(I) - 33 for I in phred])


def parse_file(filename):
    with open(filename, 'r') as f:
        content = f.readlines()
        records = []
        for line in range(len(content)):
            if content[line][0] == '@':
                record = {}
                items = content[line].split(' ')
                record['ID'] = items[0][1:]
                for item in items[1:]:
                    record[item.split('=')[0]] = item.split('=')[1].split('\n')[0]
                record['seq'] = content[line+1][:-1]
                record['phred'] = content[line+3][:-1]
#                record['quality'] = phred_to_quality(record['phred'])
            records.append(record)

    return records

if __name__ == '__main__':
    filename = r'/media/noort/Data/users/noort/20220816_1950_MN30914_AJF795_9344cc69/guppy/pass/fastq_runid_c89ca53bdde84e9f9c833a7a4d1a0dc7fd2f53a4_0_0.fastq'
    data = parse_file(filename)
    entry_nr = 9
    for key in data[entry_nr].keys():
        print(key, data[entry_nr][key])
    print('length', len(data[entry_nr]['seq']))
