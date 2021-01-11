from Bio import Entrez
import csv
import sys


def listIdAndValidateInput(gse_id):
    Entrez.email = "ranshoransho@gmail.com"
    if gse_id[:3] != 'GSE' or len(gse_id) < 4:
        print("id is not valid\n")
        sys.exit()

    handle = Entrez.esearch(db='gds', term=gse_id, idtypr='acc', retmode='text')
    record = Entrez.read(handle)
    handle.close()
    id_list = record['IdList']
    return id_list

def microarray(gse_id):

    id_list = listIdAndValidateInput(gse_id)

    handle = Entrez.esummary(db="gds", id=",".join(id_list), retmode="text")
    records = Entrez.parse(handle)
    records_list = []

    try:
        for record in records:
            record_dict = {}
            for key, val in record.items():  # Search for the desired fields
                if key == 'Accession':
                    record_dict['accession_id'] = val
                if key == 'GPL':
                    record_dict['gpls'] = val
                if key == 'suppFile':
                    record_dict['suppfile'] = val
                if key == 'FTPLink':
                    record_dict['ftplink'] = val

            records_list.append(record_dict)  # add record to list
        handle.close()

        with open('microarrayas.csv', 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=records_list[0].keys())
            writer.writeheader()  # write columns names
            for record_dict in records_list:
                writer.writerow(record_dict)
        return
    except RuntimeError or StopIteration:
        print("Experiment not found\n")
        return


def rnaseq(gse_id):
    id_list = listIdAndValidateInput(gse_id)

    handle = Entrez.elink(db="sra", dbfrom="gds", idtypr='acc', id=",".join(id_list), retmode="xml")
    records = Entrez.read(handle)
    rec_list = []
    parse_rec = []

    # Preparation after unification
    try:
        for record in records:
            parse_rec.append(record['LinkSetDb'][0]['Link'])
        for id in parse_rec:
            for i in id:
                rec_list.append(i['Id'])
        handle.close()

    except RuntimeError or StopIteration:
        print("sequence not found\n")
        return

    handle = Entrez.esummary(db="sra", id=",".join(rec_list), retmode="xml")
    records = Entrez.parse(handle)
    records_list = []

    try:
        for record in records:
            record_dict = {}
            for key, val in record.items():
                if key == 'Runs':
                    record_dict['SRRid'] = val[val.find("SRR"):val.find("total_spots") - 2]
                if key == 'ExpXml':
                    tmp_rsp = val[val.find('RNA-Seq"/><Study acc='):val.find('name="RNAseq') - 2]
                    record_dict['SRPid'] = tmp_rsp[tmp_rsp.find('SRP'):]
                    tmp_rsp = val[val.find('<Bioproject>'):val.find("</Bioproject>") - 1]
                    record_dict['Bioproject'] = tmp_rsp[tmp_rsp.find('PRJ'):]
                if key == 'ExtLinks':
                    record_dict['SRAlink'] = val
            records_list.append(record_dict)
        handle.close()

        with open('rnaseq.csv', 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=records_list[0].keys())
            writer.writeheader()  # write columns names
            for record_dict in records_list:
                writer.writerow(record_dict)
        return
    except RuntimeError or StopIteration:
        print("sequence not found\n")
        return

if __name__ == '__main__':
    microarray('GSE59847')
    rnaseq('GSE89408')
