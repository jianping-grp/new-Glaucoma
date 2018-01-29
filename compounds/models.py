# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.encoding import python_2_unicode_compatible
from django.db import models
from django.utils.translation import ugettext_lazy as _
from django_rdkit.models.fields import *
from django_rdkit.models import *
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.Crippen import MolLogP
from rdkit.Chem.Lipinski import NumHAcceptors, NumHDonors, NumRotatableBonds
from django_rdkit.config import config


# Create your models here.

class CompoundManager(models.Manager):
    def structure_search(self, smiles, similarity):
        search_mfp2 = MORGANBV_FP(Value(smiles))
        config.tanimoto_threshold = similarity
        queryset = super(CompoundManager, self).get_queryset().filter(bfp__tanimoto=search_mfp2)
        queryset = queryset.annotate(similarity=TANIMOTO_SML('bfp', search_mfp2))
        queryset = queryset.order_by('-similarity')
        return queryset

@python_2_unicode_compatible
class DrugBankID(models.Model):
    drugbank_id = models.CharField(max_length=200, blank=True, null=True)
    url = models.URLField(max_length=1024, blank=True, null=True)


    def __str__(self):
        return self.drugbank_id

    def save(self, force_insert=False, force_update=False, using=None,
             update_fields=None, *args, **kwargs):
        self.url = 'https://www.drugbank.ca/drugs/{}'.format(self.drugbank_id)
        super(DrugBankID, self).save(*args, **kwargs)




@python_2_unicode_compatible
class Drug(models.Model):

    objects = CompoundManager()

    drug_name = models.CharField(max_length=1024, null=True, blank=True)
    drug_state = models.CharField(max_length=1024, null=True, blank=True)
    drug_state_bool = models.IntegerField(null=True, blank=True)
    cas = models.CharField(max_length=1024, blank=True, null=True)
    IUPAC = models.CharField(max_length=2048, blank=True, null=True)
    drugbank_id = models.CharField(max_length=1024, null=True, blank=True)
    drugbank_url = models.URLField(max_length=1024, blank=True, null=True)
    cid = models.CharField(max_length=200, null=True, blank=True)
    first_created = models.DateTimeField(auto_now_add=True)
    last_modified = models.DateTimeField(auto_now=True)
    smiles = models.CharField(max_length=2048, null=True, blank=True)
    formula = models.CharField(max_length=1024, blank=True, null=True)

    mol = MolField(null=True, blank=True)
    mol_block = models.TextField(blank=True)
    bfp = BfpField(blank=True, null=True)
    mol_file = models.FileField(upload_to='mol_files', blank=True, null=True)
    mol_image = models.ImageField(upload_to='mol_images', blank=True, null=True)
    mol_weight = models.DecimalField(max_digits=9, decimal_places=2, blank=True, null=True)

    alogp = models.DecimalField(max_digits=9, decimal_places=2, verbose_name=_('Calculated AlogP'), null=True)
    hba = models.SmallIntegerField(verbose_name=_('Number of hydrogen bond acceptor'), blank=True, null=True)
    hbd = models.SmallIntegerField(verbose_name=_('Number of hydrogen bond donor'), blank=True, null=True)
    psa = models.DecimalField(max_digits=9, decimal_places=2, verbose_name=_('Polar surface area'), blank=True, null=True)
    rtb = models.SmallIntegerField(verbose_name=_('Number of rotatable bonds'), blank=True, null=True)
    targets = models.ManyToManyField('Target', blank=True)

    def save(self, force_insert=False, force_update=False, using=None, update_fields=None, *args, **kwargs):
        smiles = self.smiles
        if smiles:
            try:
                self.mol = Chem.MolFromSmiles(smiles)
                self.mol_block = Chem.MolToMolBlock(self.mol)
                self.mol_weight = Descriptors.ExactMolWt(self.mol)
                self.alogp = MolLogP(self.mol)
                self.hba = NumHAcceptors(self.mol)
                self.hbd = NumHDonors(self.mol)
                self.psa = Chem.MolSurf.TPSA(self.mol)
                self.rtb = NumRotatableBonds(self.mol)
                super(Drug, self).save(*args, **kwargs)
                self.formula = Chem.rdMolDescriptors.CalcMolFormula(self.mol)
                self.bfp = MORGANBV_FP(Value(smiles))
            except (ValueError, TypeError):
                print('Error when storing mol object')
                pass
        super(Drug, self).save(*args, **kwargs)

    def __str__(self):
        return self.drug_name

@python_2_unicode_compatible
class Target(models.Model):
    entry_name = models.CharField(max_length=200, blank=True, null=True)
    uniprot_accession = models.CharField(max_length=200, blank=True, null=True)
    gene = models.CharField(max_length=200, blank=True, null=True)
    protein_description = models.CharField(max_length=2048, blank=True, null=True)
    keggid = models.CharField(max_length=200, blank=True, null=True)
    kegg_url = models.URLField(max_length=2048, blank=True, null=True)
    chemblid = models.CharField(max_length=200, blank=True, null=True)
    chembl_url = models.URLField(max_length=2048, blank=True, null=True)
    type = models.CharField(max_length=2048, blank=True, null=True)
    pdbid = models.CharField(max_length=4096, blank=True, null=True)

    drugbankids = models.ManyToManyField('DrugBankID', blank=True)
    #chembl_small_molecules = models.ManyToManyField('ChEMBLSmallMolecule', blank=True)
    chembl_small_molecules_all_infos = models.ManyToManyField('ChEMBL_small_molecule_all_info', blank=True)


    def __str__(self):
        return self.entry_name
    def save(self, force_insert=False, force_update=False, using=None,
             update_fields=None, *args, **kwargs):
        self.kegg_url = 'http://www.kegg.jp/dbget-bin/www_bget?{}'.format(self.keggid)
        self.chembl_url = 'https://www.ebi.ac.uk/chembl/target/inspect/{}'.format(self.chemblid)
        super(Target, self).save(*args, **kwargs)

@python_2_unicode_compatible
class Pathway(models.Model):
    pathway_name = models.CharField(max_length=200, null=True, blank=True)
    url = models.URLField(max_length=1024, blank=True, null=True)
    descripor = models.CharField(max_length=2048, null=True, blank=True)
    drugs = models.ManyToManyField(Drug, blank=True)
    targets = models.ManyToManyField(Target, blank=True)

    def __str__(self):
        return self.pathway_name

    def save(self, force_insert=False, force_update=False, using=None,
             update_fields=None, *args, **kwargs):
        self.url = 'http://www.kegg.jp/kegg-bin/show_pathway?{}'.format(self.pathway_name)
        super(Pathway, self).save(*args, **kwargs)


class ChEMBL_small_molecule_all_info(models.Model):
    molecule_chembl_id = models.CharField(max_length=200, null=True, blank=True)
    molecule_chembl_id_url = models.URLField(max_length=1024, blank=True, null=True)
    molecule_smile = models.CharField(max_length=2048, null=True, blank=True)

    activity_pchembl_value = models.CharField(max_length=50, blank=True, null=True)
    activity_standard_value = models.CharField(max_length=50, null=True, blank=True)
    activity_standard_type = models.CharField(max_length=50, blank=True, null=True)
    activity_standard_units = models.CharField(max_length=50, blank=True, null=True)
    activity_doc_id = models.CharField(max_length=50, blank=True, null=True)
    activity_assay_id = models.CharField(max_length=50, null=True, blank=True)

    assay_chembl_id = models.CharField(max_length=50, blank=True, null=True)
    assay_chembl_id_url = models.CharField(max_length=1024, blank=True, null=True)
    assay_description = models.CharField(max_length=2048, blank=True, null=True)
    assay_assay_id = models.CharField(max_length=50, null=True, blank=True)

    doc_chembl_id = models.CharField(max_length=50, null=True, blank=True)
    doc_chembl_id_url = models.URLField(max_length=200, blank=True, null=True)
    doc_doi = models.CharField(max_length=200, null=True, blank=True)
    docs_title = models.CharField(max_length=1024, blank=True, null=True)
    doc_pubmed_id = models.CharField(max_length=50, null=True, blank=True)
    doc_pubmed_id_url = models.URLField(max_length=1024, null=True, blank=True)


    def save(self, force_insert=False, force_update=False, using=None,
             update_fields=None):
        self.molecule_chembl_id_url = 'https://www.ebi.ac.uk/chembl/compound/inspect/{}'.format(self.molecule_chembl_id)
        self.assay_chembl_id_url = 'https://www.ebi.ac.uk/chembl/assay/inspect/{}'.format(self.assay_chembl_id)
        self.doc_chembl_id_url = 'https://www.ebi.ac.uk/chembl/doc/inspect/{}'.format(self.doc_chembl_id)
        self.doc_pubmed_id_url = 'https://www.ncbi.nlm.nih.gov/pubmed/?term={}'.format(self.doc_pubmed_id)
        super(ChEMBL_small_molecule_all_info, self).save()

class Feedback(models.Model):
    username = models.CharField(max_length=256)
    email = models.EmailField()
    phone = models.CharField(max_length=64, blank=True, null=True)
    message = models.TextField()
    ip = models.GenericIPAddressField(blank=True, null=True)
    create_at = models.DateTimeField(auto_now=True)



# @python_2_unicode_compatible
# class ChEMBLSmallMolecule(models.Model):
#     chembl_id = models.CharField(max_length=100, null=True, blank=True)
#     chembl_url = models.URLField(max_length=200, null=True, blank=True)
#     # description = models.CharField(max_length=1024, blank=True, null=True)
#     canonical_smiles = models.CharField(max_length=2048, blank=True, null=True)
#     def __str__(self):
#         return self.chembl_id
#
#     def save(self, force_insert=False, force_update=False, using=None,
#              update_fields=None):
#         self.chembl_url = 'https://www.ebi.ac.uk/chembl/compound/inspect/{}'.format(self.chembl_id)
#         super(ChEMBLSmallMolecule, self).save()


# @python_2_unicode_compatible
# class Assay(models.Model):
#     chembl_id = models.CharField(max_length=50, blank=True, null=True)
#     chembl_id_url = models.URLField(max_length=2048, null=True, blank=True)
#     description = models.CharField(max_length=2048, blank=True, null=True)
#     assay_id = models.IntegerField(blank=True, null=True)
#     chembl_small_molecules = models.ManyToManyField('ChEMBLSmallMolecule',blank=True)
#
#     def save(self, force_insert=False, force_update=False, using=None,
#              update_fields=None):
#         self.chembl_id_url = 'https://www.ebi.ac.uk/chembl/assay/inspect/{}'.format(self.chembl_id)
#         super(Assay, self).save()
#
#     def __str__(self):
#         return self.chembl_id


# class Activity(models.Model):
#     pchembl_value = models.CharField(max_length=50, blank=True, null=True)
#     standard_value = models.CharField(max_length=50, null=True, blank=True)
#     standard_type = models.CharField(max_length=50, blank=True, null=True)
#     standard_units = models.CharField(max_length=50, blank=True, null=True)
#     doc_id = models.IntegerField(null=True, blank=True)
#     assay_id = models.IntegerField(null=True, blank=True)
#
#     small_molecule_chembl = models.ForeignKey('ChEMBLSmallMolecule', blank=True, null=True)

# @python_2_unicode_compatible
# class Doc(models.Model):
#     title = models.CharField(max_length=1024, blank=True, null=True)
#     doi = models.CharField(max_length=1024, blank=True, null=True)
#     pubmed_id = models.IntegerField(null=True, blank=True)
#     pubmed_id_url = models.URLField(max_length=1024, blank=True, null=True)
#     chembl_id = models.CharField(max_length=50, null=True, blank=True)
#     chembl_id_url = models.URLField(max_length=1024, blank=True, null=True)
#
#     small_molecule_chembls = models.ManyToManyField('ChEMBLSmallMolecule', blank=True)
#     activities = models.ManyToManyField('Activity', blank=True)
#
#     def save(self, force_insert=False, force_update=False, using=None,
#              update_fields=None):
#         self.pubmed_id_url = 'https://www.ncbi.nlm.nih.gov/pubmed/?term={}'.format(self.pubmed_id)
#         self.chembl_id_url = 'https://www.ebi.ac.uk/chembl/doc/inspect/{}'.format(self.chembl_id)
#         super(Doc, self).save()
#
#     def __str__(self):
#         return self.title






