# -*- coding: utf-8 -*-
# Generated by Django 1.11 on 2018-06-15 03:05
from __future__ import unicode_literals

from django.db import migrations, models
import django_rdkit.models.fields


class Migration(migrations.Migration):

    dependencies = [
        ('compounds', '0008_auto_20180613_0233'),
    ]

    operations = [
        # migrations.CreateModel(
        #     name='ChEMBL_small_molecule',
        #     fields=[
        #         ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
        #         ('molecule_chembl_id', models.CharField(blank=True, max_length=200, null=True)),
        #         ('molecule_chembl_id_url', models.URLField(blank=True, max_length=1024, null=True)),
        #         ('molecule_smile', models.CharField(blank=True, max_length=2048, null=True)),
        #         ('formula', models.CharField(blank=True, max_length=1024, null=True)),
        #         ('mol', django_rdkit.models.fields.MolField(blank=True, null=True)),
        #         ('mol_block', models.TextField(blank=True)),
        #         ('bfp', django_rdkit.models.fields.BfpField(blank=True, null=True)),
        #         ('mol_weight', models.DecimalField(blank=True, decimal_places=2, max_digits=9, null=True)),
        #         ('alogp', models.DecimalField(decimal_places=2, max_digits=9, null=True, verbose_name='Calculated AlogP')),
        #         ('hba', models.SmallIntegerField(blank=True, null=True, verbose_name='Number of hydrogen bond acceptor')),
        #         ('hbd', models.SmallIntegerField(blank=True, null=True, verbose_name='Number of hydrogen bond donor')),
        #         ('psa', models.DecimalField(blank=True, decimal_places=2, max_digits=9, null=True, verbose_name='Polar surface area')),
        #         ('rtb', models.SmallIntegerField(blank=True, null=True, verbose_name='Number of rotatable bonds')),
        #     ],
        # ),
        migrations.RemoveField(
            model_name='target',
            name='bioactivity_count',
        ),
        migrations.RemoveField(
            model_name='target',
            name='compound_count',
        ),
        migrations.RemoveField(
            model_name='target',
            name='reference_count',
        ),
        # migrations.AddField(
        #     model_name='chembl_small_molecule_all_info',
        #     name='alogp',
        #     field=models.DecimalField(decimal_places=2, max_digits=9, null=True, verbose_name='Calculated AlogP'),
        # ),
        # migrations.AddField(
        #     model_name='chembl_small_molecule_all_info',
        #     name='bfp',
        #     field=django_rdkit.models.fields.BfpField(blank=True, null=True),
        # ),
        # migrations.AddField(
        #     model_name='chembl_small_molecule_all_info',
        #     name='formula',
        #     field=models.CharField(blank=True, max_length=1024, null=True),
        # ),
        # migrations.AddField(
        #     model_name='chembl_small_molecule_all_info',
        #     name='hba',
        #     field=models.SmallIntegerField(blank=True, null=True, verbose_name='Number of hydrogen bond acceptor'),
        # ),
        # migrations.AddField(
        #     model_name='chembl_small_molecule_all_info',
        #     name='hbd',
        #     field=models.SmallIntegerField(blank=True, null=True, verbose_name='Number of hydrogen bond donor'),
        # ),
        # migrations.AddField(
        #     model_name='chembl_small_molecule_all_info',
        #     name='mol',
        #     field=django_rdkit.models.fields.MolField(blank=True, null=True),
        # ),
        # migrations.AddField(
        #     model_name='chembl_small_molecule_all_info',
        #     name='mol_block',
        #     field=models.TextField(blank=True),
        # ),
        # migrations.AddField(
        #     model_name='chembl_small_molecule_all_info',
        #     name='mol_weight',
        #     field=models.DecimalField(blank=True, decimal_places=2, max_digits=9, null=True),
        # ),
        # migrations.AddField(
        #     model_name='chembl_small_molecule_all_info',
        #     name='psa',
        #     field=models.DecimalField(blank=True, decimal_places=2, max_digits=9, null=True, verbose_name='Polar surface area'),
        # ),
        # migrations.AddField(
        #     model_name='chembl_small_molecule_all_info',
        #     name='rtb',
        #     field=models.SmallIntegerField(blank=True, null=True, verbose_name='Number of rotatable bonds'),
        # ),
        # migrations.AddField(
        #     model_name='drug',
        #     name='drug_state_bool',
        #     field=models.IntegerField(blank=True, null=True),
        # ),
        # migrations.AddField(
        #     model_name='target',
        #     name='uniprot_url',
        #     field=models.URLField(blank=True, max_length=1024, null=True),
        # ),
        # migrations.AddField(
        #     model_name='target',
        #     name='chembl_small_molecules_structure_info',
        #     field=models.ManyToManyField(blank=True, to='compounds.ChEMBL_small_molecule'),
        # ),
    ]
