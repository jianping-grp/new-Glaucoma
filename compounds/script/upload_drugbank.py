import os
os.environ.setdefault("DJANGO_SETTINGS_MODULE", 'glaucoma.settings')
import django
django.setup()
import pandas as pd

from compounds.models import *

def upload_drugbank(line):
    drugbanks = line['Drugbank'].split('\n') if line['Drugbank'] else []
    entry_name = line['Entry_name']

    target = Target.objects.get(entry_name=entry_name)

    for db in drugbanks:
        drugbank, created = DrugBankID.objects.get_or_create(drugbank_id=db)
        target.drugbankids.add(drugbank)
        target.save()

if __name__ == '__main__':
    drugbank_files = '/home/jianping/Desktop/weiyu/drugs/Drug_bank.xlsx'
    df = pd.read_excel(drugbank_files, sheet_name=0)
    for idx, line in df.iterrows():
        print(idx)
        upload_drugbank(line)

