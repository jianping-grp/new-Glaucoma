import os
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'glaucoma.settings')
import django
django.setup()

from rdkit.Chem import Draw
from compounds.models import *
from django.core.files import File

base_dir = '/home/jianping/django_test/glaucoma'
IMAGE_DIR = os.path.join(base_dir, 'media', 'mol_images')
MOL_DIR = os.path.join(base_dir, 'media', 'mol_files')

def mol_plotter():
    if not os.path.isdir(base_dir):
        os.mkdir(base_dir)
    if not os.path.isdir(IMAGE_DIR):
        os.makedirs(IMAGE_DIR)
    if not os.path.join(MOL_DIR):
        os.makedirs(MOL_DIR)

    drugs = Drug.objects.exclude(mol=None)
    for drug in drugs:
        image_name = '%s.svg' % drug.pk
        image_file = os.path.join(IMAGE_DIR, image_name)
        print(image_file)
        if not os.path.exists(image_file):
            Draw.MolToFile(drug.mol, image_file)
        img = File(open(image_file))
        drug.mol_image.save(image_name, img)

        mol_file_name = '%s.mol' % drug.pk
        mol_file = os.path.join(MOL_DIR, mol_file_name)
        print(mol_file)
        if not os.path.exists(mol_file):
            f = open(mol_file, 'w')
            f.write(drug.mol_block)
            f.close()
        mol_f = File(open(mol_file))
        drug.mol_file.save(mol_file_name, mol_f)
        drug.save()

if __name__ == '__main__':
    mol_plotter()

